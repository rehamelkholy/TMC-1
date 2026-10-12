"""
Fetch the KIDA species list (Wakelam et al. 2015, ApJS 217, 20) -- Lee et al. 2021 report 578 entries,
converted from the database's own InChI notation to SMILES.

umda reads a per-SPECIES export ("kida-molecules_05_Jul_2020.csv", with an "InChI" column), not the
reaction-network zip. The current KIDA site's equivalent is the browsable species table at
https://kida.astrochem-tools.org/species: a plain, fully server-rendered, unpaginated HTML page
(no JS framework; `?format=json` is ignored). Each row has the form:

    <tr>
        <td><a id="filter_species" href="/species/675/Al">Al</a></td>
        <td>Al</td>
        <td class="hidden-xs">InChI=1S/Al</td>
        <td><a href="/search?...">...</a></td>
    </tr>

i.e. a Bootstrap-templated table whose `<td>` tags may carry a `class` attribute (the InChI cell always
has `class="hidden-xs"`), so the row pattern accepts arbitrary attributes on every tag.

Output: external/pool/raw/kida_smiles.smi
"""
import os
import re

import http_common

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "external", "pool", "raw")
os.makedirs(OUT_DIR, exist_ok=True)
OUT_PATH = os.path.join(OUT_DIR, "kida_smiles.smi")

BASE_URL = "https://kida.astrochem-tools.org/species"

# Matches one <tr>...</tr> table row and pulls out the Species link text, the Formula cell, and the
# Inchi cell. Every tag allows arbitrary attributes ([^>]*) since the real <td> for the InChI cell
# carries class="hidden-xs".
ROW_RE = re.compile(
    r"<tr[^>]*>\s*"
    r"<td[^>]*>\s*<a[^>]*>\s*([^<]+?)\s*</a>\s*</td>\s*"  # species name (col 1)
    r"<td[^>]*>\s*([^<]*?)\s*</td>\s*"                     # formula (col 2)
    r"<td[^>]*>\s*(InChI=\S+)\s*</td>",                    # InChI (col 3) -- rest of row not needed
    re.IGNORECASE,
)


def _extract_species(html):
    """Return {species_name: inchi} parsed out of one page's HTML."""
    found = {}
    for name, _formula, inchi in ROW_RE.findall(html):
        found[name.strip()] = inchi.strip()
    return found


def main():
    from rdkit import Chem

    print(f"Fetching {BASE_URL} ...")
    try:
        html = http_common.read_text(BASE_URL, timeout=60)
    except Exception as e:
        print(f"Failed to fetch {BASE_URL}: {e}\n"
              f"Manual fallback: open that URL in your browser, save the "
              f"page's HTML (Ctrl+S -> 'Webpage, HTML only'), place it at "
              f"{os.path.join(OUT_DIR, 'kida_species.html')}, and tell me "
              f"so I can adjust this script to read the saved file instead "
              f"of fetching live.")
        return

    all_species = _extract_species(html)
    print(f"  base page: {len(all_species)} species with an InChI found")

    # Defensive pagination probe: only kicks in if the base page looks incomplete relative to the paper's
    # reported 578; harmless (a few extra requests that will 404/repeat-empty) if the site is a single page.
    if len(all_species) < 200:
        print("  base page looks smaller than expected -- probing for pagination ...")
        for page in range(2, 21):
            for pattern in (f"{BASE_URL}?page={page}",):
                try:
                    page_html = http_common.read_text(pattern, timeout=30)
                except Exception:
                    continue
                new_species = _extract_species(page_html)
                added = {k: v for k, v in new_species.items() if k not in all_species}
                if not added:
                    continue
                all_species.update(added)
                print(f"  {pattern}: +{len(added)} new (total {len(all_species)})")
        print(f"  after pagination probe: {len(all_species)} species total")

    if not all_species:
        print("No species/InChI rows matched -- the page's HTML structure "
              "must differ from what this script expects. See the manual "
              "fallback note above (save the page's HTML and send it to me).")
        return

    seen = set()
    n_parse_fail = 0
    for name, inchi in all_species.items():
        try:
            mol = Chem.MolFromInchi(inchi)
        except Exception:
            mol = None
        if mol is not None:
            try:
                seen.add(Chem.MolToSmiles(mol, canonical=True))
            except Exception:
                n_parse_fail += 1
        else:
            n_parse_fail += 1

    with open(OUT_PATH, "w") as f:
        for smi in seen:
            f.write(smi + "\n")
    print(f"Wrote {len(seen)} unique canonical SMILES to {OUT_PATH} "
          f"({n_parse_fail} InChI entries didn't convert -- normal for "
          f"bare atoms/ions like 'InChI=1S/Al' or 'e-', which RDKit can't "
          f"represent as a SMILES molecule anyway; paper reports 578 from "
          f"this source).")


if __name__ == "__main__":
    main()
