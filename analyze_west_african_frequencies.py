"""Calculate regional APOL1 G1/G2 frequencies from 1000 Genomes data."""

import json
import re
from pathlib import Path

import matplotlib
import pandas as pd
from scipy.stats import beta

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


THIS = Path(__file__).resolve()
ROOT = THIS.parent
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
FIGURES = ROOT / "figures"
PROCESSED.mkdir(parents=True, exist_ok=True)
FIGURES.mkdir(parents=True, exist_ok=True)


VARIANTS = [
    {
        "rsid": "rs73885319",
        "haplotype": "G1",
        "change": "S342G",
        "risk_allele": "G",
        "ensembl_file": "rs73885319_G1.json",
        "abraom_rsids": ("rs73885319",),
    },
    {
        "rsid": "rs71785313",
        "haplotype": "G2",
        "change": "N388_Y389del",
        "risk_allele": "ATAA",
        "ensembl_file": "rs71785313_G2.json",
        "abraom_rsids": ("rs71785313", "rs143830837"),
    },
]

POPULATIONS = {
    "ESN": ("Esan in Nigeria", "Nigeria", "West Africa"),
    "GWD": ("Gambian in Western Divisions", "The Gambia", "West Africa"),
    "MSL": ("Mende in Sierra Leone", "Sierra Leone", "West Africa"),
    "YRI": ("Yoruba in Ibadan", "Nigeria", "West Africa"),
    "ACB": ("African Caribbean in Barbados", "Barbados", "African diaspora"),
    "ASW": ("African ancestry in southwest USA", "United States", "African diaspora"),
}


def clopper_pearson(allele_count, allele_number, alpha=0.05):
    """Return the exact binomial confidence interval as percentages."""
    if allele_number == 0:
        return float("nan"), float("nan")
    low = (
        0.0
        if allele_count == 0
        else beta.ppf(alpha / 2, allele_count, allele_number - allele_count + 1)
    )
    high = (
        1.0
        if allele_count == allele_number
        else beta.ppf(
            1 - alpha / 2,
            allele_count + 1,
            allele_number - allele_count,
        )
    )
    return low * 100, high * 100


def load_1000_genomes_rows():
    """Extract risk-allele counts for named 1000 Genomes populations."""
    rows = []
    for variant in VARIANTS:
        path = RAW / "ensembl" / variant["ensembl_file"]
        data = json.loads(path.read_text(encoding="utf-8"))
        populations = data.get("populations", [])
        for code, (label, country, group) in POPULATIONS.items():
            population_id = f"1000GENOMES:phase_3:{code}"
            entries = [
                entry
                for entry in populations
                if entry.get("population") == population_id
            ]
            if not entries:
                continue
            risk_entry = next(
                (
                    entry
                    for entry in entries
                    if entry.get("allele") == variant["risk_allele"]
                ),
                None,
            )
            allele_number = sum(entry.get("allele_count", 0) for entry in entries)
            allele_count = risk_entry.get("allele_count", 0) if risk_entry else 0
            low, high = clopper_pearson(allele_count, allele_number)
            rows.append(
                {
                    "group": group,
                    "population": code,
                    "population_label": label,
                    "country": country,
                    "source": "1000 Genomes phase 3",
                    "haplotype": variant["haplotype"],
                    "rsid": variant["rsid"],
                    "protein_change": variant["change"],
                    "allele_count": allele_count,
                    "allele_number": allele_number,
                    "frequency_pct": allele_count / allele_number * 100,
                    "ci_low_pct": low,
                    "ci_high_pct": high,
                }
            )
    return rows


def load_abraom_rows():
    """Extract G1 and merged-ID G2 counts from the ABraOM SABE609 exome data."""
    path = RAW / "abraom" / "apol1_variants_wes_hg19.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    rsid_pattern = re.compile(r"rs\d+")
    by_rsid = {}
    for record in data.get("data", []):
        for rsid in rsid_pattern.findall(record.get("avsnp147", "")):
            by_rsid.setdefault(rsid, record)

    rows = []
    for variant in VARIANTS:
        record = next(
            (
                by_rsid[rsid]
                for rsid in variant["abraom_rsids"]
                if rsid in by_rsid
            ),
            None,
        )
        if record is None:
            continue
        allele_count = int(record["Allele_ALT_count"])
        allele_number = int(record["Allele_number"])
        low, high = clopper_pearson(allele_count, allele_number)
        rows.append(
            {
                "group": "Brazilian comparison",
                "population": "SABE609",
                "population_label": "ABraOM SABE609",
                "country": "Brazil",
                "source": "ABraOM",
                "haplotype": variant["haplotype"],
                "rsid": variant["rsid"],
                "protein_change": variant["change"],
                "allele_count": allele_count,
                "allele_number": allele_number,
                "frequency_pct": allele_count / allele_number * 100,
                "ci_low_pct": low,
                "ci_high_pct": high,
            }
        )
    return rows


def make_figure(frame):
    """Plot regional G1 and G2 frequencies with exact confidence intervals."""
    order = ["ESN", "GWD", "MSL", "YRI", "ACB", "ASW", "SABE609"]
    labels = [
        "ESN\nNigeria",
        "GWD\nThe Gambia",
        "MSL\nSierra Leone",
        "YRI\nNigeria",
        "ACB\nBarbados",
        "ASW\nUSA",
        "ABraOM\nSao Paulo",
    ]
    colors = {"G1": "#176b87", "G2": "#c64b32"}
    x_positions = np.arange(len(order))
    width = 0.36

    fig, axis = plt.subplots(figsize=(10.5, 5.8))
    for offset, haplotype in ((-width / 2, "G1"), (width / 2, "G2")):
        values = []
        errors = [[], []]
        for population in order:
            row = frame[
                (frame["population"] == population)
                & (frame["haplotype"] == haplotype)
            ].iloc[0]
            frequency = row["frequency_pct"]
            values.append(frequency)
            errors[0].append(frequency - row["ci_low_pct"])
            errors[1].append(row["ci_high_pct"] - frequency)
        bars = axis.bar(
            x_positions + offset,
            values,
            width,
            yerr=np.array(errors),
            capsize=3,
            color=colors[haplotype],
            label=f"{haplotype} risk allele",
        )
        for bar, value in zip(bars, values):
            axis.annotate(
                f"{value:.1f}",
                (bar.get_x() + bar.get_width() / 2, bar.get_height()),
                xytext=(0, 4),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=8,
            )

    axis.axvline(3.5, color="#9a9a9a", linewidth=1)
    axis.axvline(5.5, color="#9a9a9a", linewidth=1)
    axis.text(1.5, 60, "West African reference populations", ha="center", fontsize=9)
    axis.text(4.5, 60, "African diaspora", ha="center", fontsize=9)
    axis.text(6, 60, "Brazil", ha="center", fontsize=9)
    axis.set_ylabel("Risk-allele frequency (%)")
    axis.set_xticks(x_positions)
    axis.set_xticklabels(labels)
    axis.set_ylim(0, 64)
    axis.set_title("APOL1 G1 and G2 frequencies within West Africa and the diaspora")
    axis.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.01, 0.9))
    axis.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    output = FIGURES / "apol1_west_african_frequencies.png"
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output


def main():
    if not RAW.exists():
        raise SystemExit("Raw data not found. Run download_data.py first.")

    rows = load_1000_genomes_rows() + load_abraom_rows()
    frame = pd.DataFrame(rows)
    output_csv = PROCESSED / "apol1_regional_frequencies.csv"
    frame.to_csv(output_csv, index=False, float_format="%.6f")
    output_figure = make_figure(frame)

    display_columns = [
        "population",
        "haplotype",
        "allele_count",
        "allele_number",
        "frequency_pct",
        "ci_low_pct",
        "ci_high_pct",
    ]
    with pd.option_context("display.width", 120):
        print(frame[display_columns].to_string(index=False, float_format="%.3f"))
    print(f"\n[table] {output_csv}")
    print(f"[figure] {output_figure}")


if __name__ == "__main__":
    main()