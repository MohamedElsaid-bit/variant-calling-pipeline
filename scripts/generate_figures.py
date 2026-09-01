"""Generate summary figures for the variant calling pipeline from real pipeline output."""
import gzip
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def read_vcf_records(path):
    records = []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            records.append(line.rstrip("\n").split("\t"))
    return records

def get_info_field(info_str, key):
    for kv in info_str.split(";"):
        if kv.startswith(key + "="):
            return float(kv.split("=")[1])
    return None

raw = read_vcf_records("results/vcf/raw_variants.vcf.gz")
final = read_vcf_records("results/vcf/final_filtered.vcf.gz")

# --- Figure 1: variant counts by stage and type ---
n_raw_snp = sum(1 for r in raw if len(r[3]) == 1 and len(r[4]) == 1)
n_raw_indel = len(raw) - n_raw_snp
n_pass = sum(1 for r in final if r[6] == "PASS")
n_filtered_out = len(final) - n_pass
n_pass_snp = sum(1 for r in final if r[6] == "PASS" and len(r[3]) == 1 and len(r[4]) == 1)
n_pass_indel = n_pass - n_pass_snp

fig, ax = plt.subplots(figsize=(7, 5))
stages = ["Raw candidates\n(HaplotypeCaller)", "PASS after\nhard filtering"]
snp_counts = [n_raw_snp, n_pass_snp]
indel_counts = [n_raw_indel, n_pass_indel]
ax.bar(stages, snp_counts, label="SNPs", color="#2c5f8a")
ax.bar(stages, indel_counts, bottom=snp_counts, label="Indels", color="#c0622f")
for i, (s, ind) in enumerate(zip(snp_counts, indel_counts)):
    ax.text(i, s + ind + 1, str(s + ind), ha="center", fontweight="bold")
ax.set_ylabel("Variant count")
ax.set_title("Variant counts before and after hard filtering (chr21 test data)")
ax.legend()
plt.tight_layout()
plt.savefig("results/figures/variant_counts.png", dpi=150)
plt.close()

# --- Figure 2: quality score distribution ---
raw_quals = [float(r[5]) for r in raw if r[5] != "."]
fig, ax = plt.subplots(figsize=(7, 5))
ax.hist(raw_quals, bins=20, color="#2c5f8a", edgecolor="white")
ax.axvline(30, color="#c0622f", linestyle="--", label="QUAL = 30 reference line")
ax.set_xlabel("QUAL score")
ax.set_ylabel("Number of raw candidate variants")
ax.set_title("Distribution of variant call quality scores (pre-filter)")
ax.legend()
plt.tight_layout()
plt.savefig("results/figures/quality_distribution.png", dpi=150)
plt.close()

# --- Figure 3: PASS vs filtered breakdown ---
fig, ax = plt.subplots(figsize=(6, 6))
ax.pie(
    [n_pass, n_filtered_out],
    labels=[f"PASS ({n_pass})", f"Filtered out ({n_filtered_out})"],
    colors=["#2c5f8a", "#c0622f"],
    autopct="%1.0f%%",
    startangle=90,
)
ax.set_title("Hard-filtering outcome (GATK best-practices thresholds)")
plt.tight_layout()
plt.savefig("results/figures/filter_outcome.png", dpi=150)
plt.close()

# --- Figure 4: depth of coverage per PASS variant ---
pass_records = [r for r in final if r[6] == "PASS"]
depths = []
for r in pass_records:
    dp = get_info_field(r[7], "DP")
    if dp is not None:
        depths.append(dp)

fig, ax = plt.subplots(figsize=(7, 5))
ax.hist(depths, bins=15, color="#2c5f8a", edgecolor="white")
ax.set_xlabel("Read depth (DP) at variant site")
ax.set_ylabel("Number of PASS variants")
ax.set_title("Read depth at PASS variant sites")
plt.tight_layout()
plt.savefig("results/figures/depth_distribution.png", dpi=150)
plt.close()

print(f"Raw: {len(raw)} ({n_raw_snp} SNPs, {n_raw_indel} indels)")
print(f"Final PASS: {n_pass} ({n_pass_snp} SNPs, {n_pass_indel} indels)")
print(f"Filtered out: {n_filtered_out}")
print("Figures written to results/figures/")
