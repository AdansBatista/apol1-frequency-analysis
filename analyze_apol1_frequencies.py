"""Calculate APOL1 G1/G2 frequencies and exact confidence intervals."""

import json
import re
from pathlib import Path

import pandas as pd
from scipy.stats import beta

THIS = Path(__file__).resolve()
ROOT = THIS.parent
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
FIGURES = ROOT / "figures"
PROCESSED.mkdir(parents=True, exist_ok=True)
FIGURES.mkdir(parents=True, exist_ok=True)

VARIANTS = [
    {"rsid": "rs73885319", "haplotype": "G1", "aa": "S342G",
     "gnomad_file": "rs73885319_G1.json", "ensembl_file": "rs73885319_G1.json",
    "risk_allele": "G", "abraom_rsids": ["rs73885319"]},
    {"rsid": "rs60910145", "haplotype": "G1", "aa": "I384M",
     "gnomad_file": "rs60910145_G1.json", "ensembl_file": "rs60910145_G1.json",
    "risk_allele": "G", "abraom_rsids": ["rs60910145"]},
    {"rsid": "rs71785313", "haplotype": "G2", "aa": "N388_Y389del",
     "gnomad_file": "rs71785313_G2.json", "ensembl_file": "rs71785313_G2.json",
    "risk_allele": "ATAA", "abraom_rsids": ["rs71785313", "rs143830837"]},
]

# Exact group IDs exclude sex-stratified and compound populations.
GNOMAD_GROUPS = {
    "afr": "African / African American",
    "amr": "Admixed American (Latino)",
    "nfe": "European (non-Finnish)",
    "eas": "East Asian",
    "sas": "South Asian",
    "asj": "Ashkenazi Jewish",
    "fin": "European (Finnish)",
    "mid": "Middle Eastern",
    "ami": "Amish",
    "remaining": "Remaining",
}

KGP_SUPERPOPS = {
    "AFR": "African",
    "AMR": "Admixed American",
    "EAS": "East Asian",
    "EUR": "European",
    "SAS": "South Asian",
    "ALL": "All 1000 Genomes",
}

ALPHA = 0.05  # 95 percent CI


def clopper_pearson(ac, an, alpha=ALPHA):
    """Exact (Clopper-Pearson) binomial confidence interval as percentages."""
    if an == 0:
        return (float("nan"), float("nan"))
    lo = 0.0 if ac == 0 else beta.ppf(alpha / 2, ac, an - ac + 1)
    hi = 1.0 if ac == an else beta.ppf(1 - alpha / 2, ac + 1, an - ac)
    return (lo * 100, hi * 100)


def pct(x):
    return round(x * 100, 4)


rows = []


def load_gnomad():
    print("\n=== gnomAD v4 (exomes + genomes) ===")
    for v in VARIANTS:
        path = RAW / "gnomad" / v["gnomad_file"]
        data = json.loads(path.read_text(encoding="utf-8"))
        variant = data["_data"]["data"]["variant"]
        combined = {}  # group -> [ac, an, hom]
        for src in ("exome", "genome"):
            blk = variant.get(src)
            if not blk:
                continue
            pops = {p["id"]: p for p in (blk.get("populations") or [])}
            for gid, glabel in GNOMAD_GROUPS.items():
                p = pops.get(gid)
                if not p:
                    continue
                ac, an = p["ac"], p["an"]
                hom = p.get("homozygote_count")
                if an > 0:
                    lo, hi = clopper_pearson(ac, an)
                    rows.append({
                        "rsid": v["rsid"], "haplotype": v["haplotype"], "aa_change": v["aa"],
                        "source": f"gnomAD v4 {src}s", "population": gid,
                        "population_label": glabel, "ac": ac, "an": an,
                        "af_pct": pct(ac / an), "ci_low_pct": round(lo, 4),
                        "ci_high_pct": round(hi, 4), "homozygote_count": hom,
                    })
                c = combined.setdefault(gid, [0, 0, 0])
                c[0] += ac
                c[1] += an
                c[2] += hom or 0
        for gid, (ac, an, hom) in combined.items():
            if an > 0:
                lo, hi = clopper_pearson(ac, an)
                rows.append({
                    "rsid": v["rsid"], "haplotype": v["haplotype"], "aa_change": v["aa"],
                    "source": "gnomAD v4 combined", "population": gid,
                    "population_label": GNOMAD_GROUPS[gid], "ac": ac, "an": an,
                    "af_pct": pct(ac / an), "ci_low_pct": round(lo, 4),
                    "ci_high_pct": round(hi, 4), "homozygote_count": hom,
                })
        afr = combined.get("afr")
        amr = combined.get("amr")
        if afr and amr:
            print(f"  {v['rsid']} ({v['haplotype']} {v['aa']}): "
                  f"AFR {pct(afr[0]/afr[1])}%  AMR {pct(amr[0]/amr[1])}%  "
                  f"(AFR hom={afr[2]}, AMR hom={amr[2]})")


def load_ensembl():
    print("\n=== 1000 Genomes phase 3 (Ensembl REST) ===")
    for v in VARIANTS:
        path = RAW / "ensembl" / v["ensembl_file"]
        data = json.loads(path.read_text(encoding="utf-8"))
        risk = v["risk_allele"]
        pops = data.get("populations", [])
        for sp, label in KGP_SUPERPOPS.items():
            tag = f"1000GENOMES:phase_3:{sp}"
            entries = [p for p in pops if p.get("population") == tag]
            if not entries:
                continue
            risk_entry = next((p for p in entries if p.get("allele") == risk), None)
            total = sum(e.get("allele_count", 0) for e in entries)
            ac = risk_entry.get("allele_count", 0) if risk_entry else 0
            freq = risk_entry.get("frequency", 0.0) if risk_entry else 0.0
            lo, hi = clopper_pearson(ac, total) if total else (float("nan"), float("nan"))
            rows.append({
                "rsid": v["rsid"], "haplotype": v["haplotype"], "aa_change": v["aa"],
                "source": "1000 Genomes (Ensembl)", "population": sp,
                "population_label": label, "ac": ac, "an": total,
                "af_pct": round(freq * 100, 4) if freq is not None else None,
                "ci_low_pct": round(lo, 4), "ci_high_pct": round(hi, 4),
                "homozygote_count": None,
            })
        afr = next((r for r in rows if r["rsid"] == v["rsid"]
                    and r["source"].startswith("1000") and r["population"] == "AFR"), None)
        amr = next((r for r in rows if r["rsid"] == v["rsid"]
                    and r["source"].startswith("1000") and r["population"] == "AMR"), None)
        if afr and amr:
            print(f"  {v['rsid']} ({v['haplotype']}): risk allele '{risk}'  "
                  f"AFR {afr['af_pct']}%  AMR {amr['af_pct']}%")


def load_abraom():
    print("\n=== ABraOM (Brazilian SABE609 reference, hg19 WES) ===")
    path = RAW / "abraom" / "apol1_variants_wes_hg19.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    table = data.get("data", [])
    rsid_re = re.compile(r"rs\d+")
    by_rsid = {}
    for row in table:
        snp_field = row.get("avsnp147", "")
        found = rsid_re.findall(snp_field)
        for rs in found:
            by_rsid.setdefault(rs, row)
    for v in VARIANTS:
        row = next((by_rsid[rsid] for rsid in v["abraom_rsids"] if rsid in by_rsid), None)
        if not row:
            print(f"  {v['rsid']} ({v['haplotype']}): not found in ABraOM WES")
            continue
        ac = int(row.get("Allele_ALT_count", 0))
        an = int(row.get("Allele_number", 0))
        cohort = row.get("Cohort", "ABraOM")
        if an > 0:
            lo, hi = clopper_pearson(ac, an)
            rows.append({
                "rsid": v["rsid"], "haplotype": v["haplotype"], "aa_change": v["aa"],
                "source": "ABraOM (Brazil)", "population": cohort,
                "population_label": "Brazilian (SABE elderly cohort)", "ac": ac, "an": an,
                "af_pct": pct(ac / an), "ci_low_pct": round(lo, 4),
                "ci_high_pct": round(hi, 4),
                "homozygote_count": int(row.get("HomozygousALT_count", 0)),
            })
            print(f"  {v['rsid']} ({v['haplotype']} {v['aa']}): "
                  f"{pct(ac/an)}%  (AC={ac}, AN={an}, cohort={cohort})")


def make_figure(df):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    g1_rsid, g2_rsid = "rs73885319", "rs71785313"

    def get(source, rsid, pop):
        sub = df[(df["source"] == source) & (df["rsid"] == rsid) & (df["population"] == pop)]
        if sub.empty:
            return None
        r = sub.iloc[0]
        return (r["af_pct"], r["ci_low_pct"], r["ci_high_pct"])

    cols = [
        ("African\n(gnomAD AFR)", "gnomAD v4 combined", "afr"),
        ("Admixed American\n(gnomAD AMR)", "gnomAD v4 combined", "amr"),
        ("Brazilian\n(ABraOM)", "ABraOM (Brazil)", "SABE609"),
        ("European\n(gnomAD NFE)", "gnomAD v4 combined", "nfe"),
    ]
    labels, g1_vals, g1_err, g2_vals, g2_err = [], [], [], [], []
    g1_missing, g2_missing = [], []
    for label, source, pop in cols:
        labels.append(label)
        for rsid, vals, errs, miss in ((g1_rsid, g1_vals, g1_err, g1_missing),
                                       (g2_rsid, g2_vals, g2_err, g2_missing)):
            got = get(source, rsid, pop)
            if got is None:
                vals.append(0.0)
                errs.append([0.0, 0.0])
                miss.append(True)
            else:
                af, lo, hi = got
                vals.append(af)
                errs.append([max(af - lo, 0), max(hi - af, 0)])
                miss.append(False)

    x = np.arange(len(labels))
    width = 0.38
    fig, ax = plt.subplots(figsize=(9, 5.5))
    g1_err_t = np.array(g1_err).T
    g2_err_t = np.array(g2_err).T
    b1 = ax.bar(x - width / 2, g1_vals, width, yerr=g1_err_t, capsize=4,
                label="G1 (rs73885319, S342G)", color="#2c6fbb")
    b2 = ax.bar(x + width / 2, g2_vals, width, yerr=g2_err_t, capsize=4,
                label="G2 (rs71785313, N388_Y389del)", color="#c8542b")
    ax.set_ylabel("Risk-allele frequency (%)")
    ax.set_title("APOL1 G1 and G2 risk-allele frequencies across populations\n"
                 "gnomAD v4 and ABraOM, error bars are 95% Clopper-Pearson CIs")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()
    for bars, miss in ((b1, g1_missing), (b2, g2_missing)):
        for rect, is_missing in zip(bars, miss):
            h = rect.get_height()
            if is_missing:
                ax.annotate("n/a", xy=(rect.get_x() + rect.get_width() / 2, 0),
                            xytext=(0, 3), textcoords="offset points",
                            ha="center", va="bottom", fontsize=8, color="grey")
            elif h > 0:
                ax.annotate(f"{h:.1f}", xy=(rect.get_x() + rect.get_width() / 2, h),
                            xytext=(0, 3), textcoords="offset points",
                            ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    out = FIGURES / "apol1_g1_g2_frequencies.png"
    fig.savefig(out, dpi=200)
    print(f"\n[figure] {out}")


def main():
    if not RAW.exists():
        raise SystemExit("Raw data not found. Run download_data.py first.")

    load_gnomad()
    load_ensembl()
    load_abraom()

    df = pd.DataFrame(rows)
    order = ["rsid", "haplotype", "aa_change", "source", "population",
             "population_label", "ac", "an", "af_pct", "ci_low_pct",
             "ci_high_pct", "homozygote_count"]
    df = df[order].sort_values(["haplotype", "rsid", "source", "population"]).reset_index(drop=True)

    out_csv = PROCESSED / "apol1_allele_frequencies.csv"
    df.to_csv(out_csv, index=False)
    print(f"\n[table] {out_csv}  ({len(df)} rows)")

    print("\n" + "=" * 70)
    print("HEADLINE: G1 (rs73885319) and G2 (rs71785313) risk-allele frequency")
    print("=" * 70)
    headline = df[df["rsid"].isin(["rs73885319", "rs71785313"])
                  & df["source"].isin(["gnomAD v4 combined", "ABraOM (Brazil)",
                                        "1000 Genomes (Ensembl)"])
                  & df["population"].isin(["afr", "amr", "nfe", "SABE609", "AFR", "AMR"])]
    cols_show = ["haplotype", "rsid", "source", "population", "af_pct",
                 "ci_low_pct", "ci_high_pct", "homozygote_count"]
    with pd.option_context("display.max_rows", None, "display.width", 120):
        print(headline[cols_show].to_string(index=False))

    make_figure(df)
    print("\nDone.")


if __name__ == "__main__":
    main()
