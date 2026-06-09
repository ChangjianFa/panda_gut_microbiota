"""
Year 2 microbiome data preprocessing script
Extract Year 2 samples from asv_abundance_year2.tsv -> genus level -> filter -> monthly average
DP=Dec, QP=Oct, XP=Aug, CP=Mar
"""
import pandas as pd
import numpy as np
import re

# ============================================
# 1. Load Year 2 ASV data
# ============================================
print("===== Loading Year 2 ASV data =====")
asv = pd.read_csv("data_raw/asv_abundance_year2.tsv", sep="\t")
# ASV ID in column 1, Consensus Lineage in last column
asv_cols = [c for c in asv.columns if c not in ["#OTU ID", "Consensus Lineage"]]

# Year 1 samples in this file (DD/FX/GG/JN/MD/ME/ML format)
y1_pattern = re.compile(r'[A-Z]{2}\d{6}_\d')
y1_cols = [c for c in asv_cols if y1_pattern.match(c)]
# Year 2 samples (DP/QP/XP/CP)
y2_cols = [c for c in asv_cols if c not in y1_cols]

print(f"Year 2 samples: {len(y2_cols)} -> {sorted(y2_cols)}")

# ============================================
# 2. Parse genus names (Consensus Lineage)
# ============================================
def parse_genus(lineage):
    """Extract genus name g__xxx from Consensus Lineage"""
    parts = lineage.split(";")
    for p in parts:
        if p.startswith("g__"):
            name = p[3:]
            # Keep underscore format (consistent with Year 1)
            return name
    # If no genus-level annotation, use the last taxonomic level
    last = parts[-1]
    if last.startswith("f__"):  # family level
        return "unclassified_" + last[3:]
    return "Unknown"

asv["genus"] = asv["Consensus Lineage"].apply(parse_genus)

# ============================================
# 3. Genus-level aggregation (read counts)
# ============================================
print(f"Unique genera: {asv['genus'].nunique()}")

genus_counts = asv.groupby("genus")[y2_cols].sum()
print(f"Genera after aggregation: {genus_counts.shape[0]}")
print(f"Total reads: {genus_counts.values.sum():,.0f}")

# ============================================
# 4. Relative abundance + filter
# ============================================
genus_relab = genus_counts.div(genus_counts.sum(axis=0), axis=1)

mean_abund = genus_relab.mean(axis=1)
prevalence = (genus_relab > 0).mean(axis=1)
keep = (mean_abund > 0.0005) & (prevalence > 0.05)

bacteria = genus_relab[keep]
print(f"Genera after filter (mean>0.05% & prev>5%): {bacteria.shape[0]}")
print(f"Cumulative abundance: {mean_abund[keep].sum()*100:.2f}%")

# ============================================
# 5. Aggregate by month (by panda individual)
# Panda -> Month: DP=12, QP=10, XP=8, CP=3
# ============================================
panda_month = {"DP": 12, "QP": 10, "XP": 8, "CP": 3}
ordered_months = [3, 8, 10, 12]  # Spring Summer Autumn Winter
ordered_pandas = ["CP", "XP", "QP", "DP"]

# Avoid duplicate pandas issue in ordered_pandas
col_month = {}
for col in bacteria.columns:
    pfx = re.match(r'([A-Z]+)', col).group(1)
    if pfx in panda_month:
        col_month[col] = panda_month[pfx]

print(f"\nSample -> Month mapping:")
for pfx in ordered_pandas:
    cols = [c for c in bacteria.columns if c.startswith(pfx)]
    print(f"  {pfx} ({panda_month[pfx]}): {cols}")

# Monthly average
bacteria_month = pd.DataFrame(index=bacteria.index)
for m in ordered_months:
    cols = [c for c, mo in col_month.items() if mo == m]
    bacteria_month[m] = bacteria[cols].mean(axis=1)

bacteria_month.columns = [f"Month_{m}" for m in ordered_months]
print(f"\nBacteria month matrix: {bacteria_month.shape[0]} genera x {bacteria_month.shape[1]} months")

# ============================================
# 6. Metadata processing
# ============================================
print("\n===== Loading Year 2 metadata =====")
meta = pd.read_csv("data_raw/metadata_year2.csv", encoding="gbk")
meta.columns = ["Sample", "Age", "Sex", "Weight_kg", "Temperature",
                "Crude_Fat", "Crude_Protein", "Crude_Fiber", "Carbohydrate"]
# Keep only Year 2 samples present in ASV data
meta = meta[meta["Sample"].isin(y2_cols)]
meta["Panda"] = meta["Sample"].str.extract(r'([A-Z]+)')[0]
meta["Month"] = meta["Panda"].map(panda_month)
meta.to_csv("data_processed/metadata_year2_clean.csv", index=False, encoding="utf-8-sig")

# Monthly environmental averages
env_cols = ["Temperature", "Crude_Fat", "Crude_Protein", "Crude_Fiber", "Carbohydrate"]
meta_month = meta.groupby("Month")[env_cols].mean().T
meta_month.columns = [f"Month_{int(c)}" for c in meta_month.columns]
meta_month.index = ["Temperature", "Crude_Fat", "Crude_Protein", "Crude_Fiber", "Carbohydrate"]
print("Environmental monthly averages:")
print(meta_month.to_string())

# ============================================
# 7. Merge & output
# ============================================
all_month = pd.concat([bacteria_month, meta_month], axis=0)
p_v = bacteria_month.shape[0]
p_e = meta_month.shape[0]
p = p_v + p_e

print(f"\nVariables: total={p}, bacteria={p_v}, env={p_e}")

bacteria_month.to_csv("data_processed/bacteria_month_year2.csv", encoding="utf-8-sig")
all_month.to_csv("data_processed/all_month_year2.csv", encoding="utf-8-sig")

# ============================================
# 8. Shared genus analysis with Year 1
# ============================================
y1_bacteria = pd.read_csv("data_processed/bacteria_month.csv", index_col=0, encoding="utf-8-sig")
# Normalize genus names: replace - with _ (Year 1 uses dash, Year 2 uses underscore)
y1_names = set(y1_bacteria.index.str.replace("-", "_"))
y2_names = set(bacteria_month.index)

shared = y1_names & y2_names
print(f"\nYear 1 genera: {len(y1_names)}")
print(f"Year 2 genera: {len(y2_names)}")
print(f"Shared genera:  {len(shared)}")
print(f"Year 1 only:    {len(y1_names - y2_names)}")
print(f"Year 2 only:    {len(y2_names - y1_names)}")
print(f"\nShared: {sorted(shared)}")
print(f"Year 1 only: {sorted(y1_names - y2_names)}")
print(f"Year 2 only: {sorted(y2_names - y1_names)}")

print("\n===== Done =====")
