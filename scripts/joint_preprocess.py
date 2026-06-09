"""
Joint preprocessing script (intersection mode): Year 1 + Year 2 filtered independently, take intersection genera
"""
import pandas as pd
import numpy as np
import re

THRESH_MEAN = 0.0005   # >0.05%
THRESH_PREV = 0.05     # >5%

# ============================================================
# 1. Year 1: raw ASV -> genus level -> filter
# ============================================================
print("===== Processing Year 1 =====")
asv1 = pd.read_csv("data_raw/asv_abundance_year1.xls", sep="\t")
tax_cols = ["domain","kingdom","phylum","class","order","family","genus","species","asv"]
meta_cols = ["Total","Percent","Prevalence"]
s1 = [c for c in asv1.columns if c not in tax_cols and c not in meta_cols]

g1 = asv1.groupby("genus")[s1].sum()
g1.index = g1.index.str.replace(".*g__", "", regex=True)
g1_rel = g1.div(g1.sum(axis=0), axis=1)

keep1 = (g1_rel.mean(axis=1) > THRESH_MEAN) & ((g1_rel > 0).mean(axis=1) > THRESH_PREV)
g1_filt = g1_rel[keep1]
print(f"  {g1.shape[0]} genera → {g1_filt.shape[0]} after filter")
print(f"  Cumulative abundance: {g1_rel.mean(axis=1)[keep1].sum()*100:.2f}%")

# ============================================================
# 2. Year 2: TSV -> genus level -> filter (CP/DP/QP/XP only)
# ============================================================
print("\n===== Processing Year 2 =====")
asv2 = pd.read_csv("data_raw/asv_abundance_year2.tsv", sep="\t")

def parse_genus(lineage):
    parts = lineage.split(";")
    for p in parts:
        if p.startswith("g__"):
            return p[3:]
    last = parts[-1]
    if last.startswith("f__"):
        return "unclassified_" + last[3:]
    return "Unknown"

asv2["genus"] = asv2["Consensus Lineage"].apply(parse_genus)

y2_cols = [c for c in asv2.columns
           if c not in ["#OTU ID", "Consensus Lineage", "genus"]
           and not re.match(r'[A-Z]{2}\d{6}_\d', c)]

g2 = asv2.groupby("genus")[y2_cols].sum()
g2_rel = g2.div(g2.sum(axis=0), axis=1)

keep2 = (g2_rel.mean(axis=1) > THRESH_MEAN) & ((g2_rel > 0).mean(axis=1) > THRESH_PREV)
g2_filt = g2_rel[keep2]
print(f"  {g2.shape[0]} genera → {g2_filt.shape[0]} after filter")
print(f"  Cumulative abundance: {g2_rel.mean(axis=1)[keep2].sum()*100:.2f}%")

# ============================================================
# 3. Genus name normalization -> intersection
# ============================================================
def norm(name):
    name = name.replace("_", "-").replace(" ", "")
    name = re.sub(r'^unclassified-[a-z]--', 'unclassified-', name)
    return name

y1_names = {norm(n) for n in g1_filt.index}
y2_names = {norm(n) for n in g2_filt.index}

shared_names = sorted(y1_names & y2_names)
print(f"\n===== Intersection =====")
print(f"  Year 1 genera: {len(y1_names)}")
print(f"  Year 2 genera: {len(y2_names)}")
print(f"  Intersection:   {len(shared_names)}")
print(f"  Year 1 only:    {sorted(y1_names - y2_names)}")
print(f"  Year 2 only:    {sorted(y2_names - y1_names)}")

# Build mapping: norm_name -> original_name
y1_map = {norm(n): n for n in g1_filt.index}
y2_map = {norm(n): n for n in g2_filt.index}

# ============================================================
# 4. Extract both years with shared genera -> aggregate by month
# ============================================================
panda_month = {"CP": 3, "XP": 8, "QP": 10, "DP": 12}

# --- Year 1 monthly aggregation ---
y1_shared = g1_rel.loc[[y1_map[n] for n in shared_names]].copy()
y1_shared.index = shared_names  # Unified genus names
y1_shared.columns = y1_shared.columns.str.replace(r"(.*)(\d{6})_(\d+)", r"\2_\1\3", regex=True)
y1_months = {c: int(c[4:6]) for c in y1_shared.columns}

y1_monthly = pd.DataFrame(index=y1_shared.index)
for m in sorted(set(y1_months.values())):
    cols = [c for c, mo in y1_months.items() if mo == m]
    y1_monthly[m] = y1_shared[cols].mean(axis=1)
y1_monthly.columns = [f"Month_{m}" for m in sorted(set(y1_months.values()))]
print(f"\n  Year 1 monthly: {y1_monthly.shape[0]} genera × {y1_monthly.shape[1]} months")

# --- Year 2 monthly aggregation ---
y2_shared = g2_rel.loc[[y2_map[n] for n in shared_names]].copy()
y2_shared.index = shared_names

y2_monthly = pd.DataFrame(index=y2_shared.index)
for pfx, mo in sorted(panda_month.items(), key=lambda x: x[1]):
    cols = [c for c in y2_shared.columns if c.startswith(pfx)]
    y2_monthly[mo] = y2_shared[cols].mean(axis=1)
y2_monthly.columns = [f"Month_{m}" for m in sorted(panda_month.values())]
print(f"  Year 2 monthly: {y2_monthly.shape[0]} genera × {y2_monthly.shape[1]} months")

# ============================================================
# 5. Environmental data
# ============================================================
# Year 1 metadata
meta1 = pd.read_csv("data_raw/metadata_year1.csv", encoding="utf-8-sig")
meta1 = meta1[~meta1.iloc[:, 2].isna() & (meta1.iloc[:, 2] != "")]
meta1.columns = list(range(len(meta1.columns)))
meta1 = meta1.rename(columns={6:"Temperature",7:"Crude_Fat",8:"Crude_Protein",
                               9:"Crude_Fiber",10:"Carbohydrate"})
meta1.iloc[:, 0] = meta1.iloc[:, 0].astype(int)
m1 = meta1.groupby(meta1.columns[0])[
    ["Temperature","Crude_Fat","Crude_Protein","Crude_Fiber","Carbohydrate"]].mean().T
m1.columns = [f"Month_{int(c)}" for c in m1.columns]
m1.index = ["Temperature","Crude_Fat","Crude_Protein","Crude_Fiber","Carbohydrate"]

# Year 2 metadata
meta2 = pd.read_csv("data_raw/metadata_year2.csv", encoding="gbk")
meta2.columns = ["Sample","Age","Sex","Weight","Temperature",
                 "Crude_Fat","Crude_Protein","Crude_Fiber","Carbohydrate"]
meta2 = meta2[meta2["Sample"].isin(y2_cols)]
meta2["Panda"] = meta2["Sample"].str.extract(r'([A-Z]+)')[0]
meta2["Month"] = meta2["Panda"].map(panda_month)
m2 = meta2.groupby("Month")[
    ["Temperature","Crude_Fat","Crude_Protein","Crude_Fiber","Carbohydrate"]].mean().T
m2.columns = [f"Month_{int(c)}" for c in m2.columns]
m2.index = ["Temperature","Crude_Fat","Crude_Protein","Crude_Fiber","Carbohydrate"]

# ============================================================
# 6. Merge & output
# ============================================================
all_y1 = pd.concat([y1_monthly, m1])
all_y2 = pd.concat([y2_monthly, m2])

p_v = y1_monthly.shape[0]
p_e = m1.shape[0]
print(f"\nFinal: {p_v} bacteria + {p_e} env = {p_v + p_e} variables")
print(f"Genera: {list(y1_monthly.index)}")

y1_monthly.to_csv("data_processed/bacteria_month.csv", encoding="utf-8-sig")
all_y1.to_csv("data_processed/all_month.csv", encoding="utf-8-sig")
y2_monthly.to_csv("data_processed/bacteria_month_year2.csv", encoding="utf-8-sig")
all_y2.to_csv("data_processed/all_month_year2.csv", encoding="utf-8-sig")

# Top genera comparison
print("\nYear 1 top 10:")
for g, v in y1_monthly.mean(axis=1).sort_values(ascending=False).head(10).items():
    print(f"  {g:30s} {v:.4f}")
print("\nYear 2 top 10:")
for g, v in y2_monthly.mean(axis=1).sort_values(ascending=False).head(10).items():
    print(f"  {g:30s} {v:.4f}")

print("\n===== Done =====")
