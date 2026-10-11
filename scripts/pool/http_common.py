"""
Shared HTTP helper for the fetch_*.py scripts.

files.docking.org (and possibly other hosts) reject requests carrying Python's default urllib User-Agent
string ("Python-urllib/3.x") with a 403, an ordinary anti-bot rule. Every fetch script therefore goes
through this module, which sets a browser-like User-Agent.

download_to_file() downloads to a temporary path first and only os.replace()s it onto out_path once the
transfer finishes cleanly; on any exception the temp file is removed and out_path is left untouched
(either absent -- eligible for a retry -- or the last known-good copy). Without this, a connection that
died mid-transfer would leave a truncated body on disk, which a resumability check of the form "does the
file already exist" would treat as complete, silently parsing fewer SMILES than the file should contain.
"""
import os
import urllib.request

# A generic, current desktop-Chrome-shaped User-Agent: the goal is just "looks like an ordinary browser
# request".
_HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
        "(KHTML, like Gecko) Chrome/124.0.0.0 Safari/537.36"
    ),
    "Accept": "*/*",
}


def urlopen(url, timeout=60):
    """Like urllib.request.urlopen, but with browser-like headers set."""
    req = urllib.request.Request(url, headers=_HEADERS)
    return urllib.request.urlopen(req, timeout=timeout)


def download_to_file(url, out_path, timeout=60):
    """Like urllib.request.urlretrieve, but with browser-like headers set (urlretrieve itself cannot pass
    custom headers), and safe against leaving a corrupt/partial file at out_path if the connection dies
    mid-transfer. Downloads to out_path + '.part' (a per-attempt temp name, safe for concurrent callers
    targeting different out_paths) and only renames it onto out_path once the full body has been written.

    Also checks the transferred byte count against the server's Content-Length (when it sent one) and
    raises if they do not match: a connection that closes mid-body does not always surface as a Python
    exception (a graceful short read can look identical to a normal end-of-body), so this turns that
    silent case into an explicit, retryable failure."""
    tmp_path = out_path + ".part"
    try:
        with urlopen(url, timeout=timeout) as resp, open(tmp_path, "wb") as out_f:
            expected_len = resp.headers.get("Content-Length")
            expected_len = int(expected_len) if expected_len is not None else None
            n_written = 0
            while True:
                chunk = resp.read(1 << 20)
                if not chunk:
                    break
                out_f.write(chunk)
                n_written += len(chunk)
            if expected_len is not None and n_written != expected_len:
                raise IOError(
                    f"truncated download from {url}: got {n_written} bytes, "
                    f"server declared Content-Length={expected_len}"
                )
        os.replace(tmp_path, out_path)  # atomic on both POSIX and Windows
    except BaseException:
        try:
            os.remove(tmp_path)
        except OSError:
            pass
        raise


def read_text(url, timeout=60, encoding="utf-8", errors="ignore"):
    """Fetch a URL and return its body decoded as text."""
    with urlopen(url, timeout=timeout) as resp:
        return resp.read().decode(encoding, errors=errors)
