"""Post hoc descriptive analysis for Discussion 4.2.
No label, feature or protocol change; touches no sealed result.
Usage: python 30_descriptive_discussion.py <candidates.csv> <output.json>
"""
import sys, json, hashlib, warnings
import numpy as np, pandas as pd, scipy, sklearn
from scipy.stats import mannwhitneyu, percentileofscore, spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

inp, outp = sys.argv[1], sys.argv[2]
print("sklearn", sklearn.__version__, "| scipy", scipy.__version__,
      "| pandas", pd.__version__, "| numpy", np.__version__)
print("input SHA-256", hashlib.sha256(open(inp, "rb").read()).hexdigest())

df = pd.read_csv(inp)
assert len(df) == 1434
F = ["dist", "sigma", "logcd", "mu", "mw", "charged", "radical"]
df["charged"] = df["charged"].astype(int)
df["radical"] = df["radical"].astype(int)
df["positive"] = df["positive"].astype(bool)
df["date_dt"] = pd.to_datetime(df["date"])
X = StandardScaler().fit(df[F].values).transform(df[F].values)

cut = pd.Timestamp("2022-01-01")
known = df[df.positive & (df.date_dt < cut)]
test = df[df.positive & (df.date_dt >= cut)]
never = df[~df.positive]
n = [len(known), len(test), len(never)]
assert n == [34, 29, 1371], n

# Step 1 -- coefficients
coef, conv_warnings = {}, []
for key, npos_exp in [("2022-01-01", 34), ("2023-01-01", 44), ("2024-01-01", 50),
                      ("2025-01-01", 59), ("final", 63)]:
    y = df.positive.values if key == "final" else \
        (df.positive & (df.date_dt < pd.Timestamp(key))).values
    assert int(y.sum()) == npos_exp, (key, int(y.sum()))
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        m = LogisticRegression(C=1.0, max_iter=5000).fit(X, y.astype(int))
    conv_warnings += [f"{key}: {str(x.message)}" for x in w]
    c = m.coef_[0]
    coef[key] = {"n_pos": int(y.sum()), "n_iter": int(m.n_iter_[0])}
    coef[key].update({f: round(float(v), 4) for f, v in zip(F, c)})
print("convergence warnings:", conv_warnings if conv_warnings else "none")

# Step 2 -- descriptives
desc = {}
for f in ["sigma", "dist", "logcd", "mw", "mu"]:
    nv = never[f].values
    desc[f] = {
        "med_known": round(float(known[f].median()), 4),
        "med_test": round(float(test[f].median()), 4),
        "med_never": round(float(never[f].median()), 4),
        "mwu_p_test_never": float(f"{mannwhitneyu(test[f].values, nv, alternative='two-sided').pvalue:.3g}"),
        "median_pct_test_in_never": round(float(np.median(
            [percentileofscore(nv, v, kind="rank") for v in test[f].values])), 2),
    }

# Step 3 -- shares, sigma counts, Spearman
share = {k: [round(float(g[k].mean()), 4) for g in (known, test, never)]
         for k in ["charged", "radical"]}
q90, q25 = never["sigma"].quantile(0.90), never["sigma"].quantile(0.25)
sigma_counts = {"test_top10pct": int((test["sigma"] >= q90).sum()),
                "test_bottom25pct": int((test["sigma"] <= q25).sum())}
spearman = {k: round(float(spearmanr(df["sigma"], df[k]).statistic), 4)
            for k in ["dist", "mw", "logcd"]}
print("never-sigma q90 =", q90, " q25 =", q25)

res = {"coef": coef, "n": n, "desc": desc, "share": share,
       "sigma_counts": sigma_counts, "spearman": spearman}
with open(outp, "w") as fh:
    json.dump(res, fh, indent=1, sort_keys=True)
print("wrote", outp)
