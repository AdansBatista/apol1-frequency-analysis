"""Calculate APOL1 frequencies from a Galaxy-extracted 1000 Genomes VCF."""

import argparse
from pathlib import Path

THIS = Path(__file__).resolve()
ROOT = THIS.parent
PANEL = ROOT / "data" / "raw" / "kgp" / "sample_metadata_phase3.txt"
DEFAULT_VCF = ROOT / "data" / "galaxy" / "apol1_region_1kg.vcf"

TARGETS = [
    {"rsid": "rs73885319", "hap": "G1", "aa": "S342G", "pos": 36265860, "ref": "A", "alt": "G"},
    {"rsid": "rs60910145", "hap": "G1", "aa": "I384M", "pos": 36265988, "ref": "T", "alt": "G"},
        {"rsid": "rs71785313", "hap": "G2", "aa": "N388_Y389del", "pos": 36265995, "ref": "AATAATT", "alt": "A"},
]


def load_panel():
    sample2pop = {}
    if not PANEL.exists():
        raise SystemExit("Sample panel not found. Run download_data.py first.")
    for i, line in enumerate(PANEL.read_text(encoding="utf-8").splitlines()):
        parts = line.split()
        if i == 0 and parts and parts[0].lower() == "sample":
            continue
        if len(parts) >= 3:
            sample2pop[parts[0]] = parts[2]  # super_pop
    return sample2pop


def tally(record_fields, sample_cols, sample2pop, alt_index=1):
    """Count alt-allele copies and called alleles overall and per super-population."""
    fmt = record_fields[8].split(":")
    gt_i = fmt.index("GT") if "GT" in fmt else 0
    overall = [0, 0]  # ac, an
    perpop = {}
    for col, sid in sample_cols:
        gt = record_fields[col].split(":")[gt_i]
        for a in gt.replace("|", "/").split("/"):
            if a == ".":
                continue
            overall[1] += 1
            if a == str(alt_index):
                overall[0] += 1
            pop = sample2pop.get(sid)
            if pop:
                pp = perpop.setdefault(pop, [0, 0])
                pp[1] += 1
                if a == str(alt_index):
                    pp[0] += 1
    return overall, perpop


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vcf", nargs="?", type=Path, default=DEFAULT_VCF)
    args = parser.parse_args()
    vcf_path = args.vcf
    if not vcf_path.exists():
        raise SystemExit(f"Region VCF not found: {vcf_path}")
    sample2pop = load_panel()
    print(f"Panel samples with super-pop label: {len(sample2pop)}")

    sample_cols = []  # (column_index, sample_id)
    records = []
    for line in vcf_path.read_text(encoding="utf-8").splitlines():
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            header = line.split("\t")
            for idx in range(9, len(header)):
                sample_cols.append((idx, header[idx]))
            continue
        if line.strip():
            records.append(line.split("\t"))
    print(f"Region variants: {len(records)} | sample columns: {len(sample_cols)}")

    print("\n" + "=" * 78)
    print("APOL1 G1/G2 frequencies from 1000G high-coverage (Galaxy-extracted region)")
    print("=" * 78)
    superpops = ["AFR", "AMR", "EAS", "EUR", "SAS"]
    for t in TARGETS:
        rec = next(
            (
                r
                for r in records
                if int(r[1]) == t["pos"] and r[3] == t["ref"] and t["alt"] in r[4].split(",")
            ),
            None,
        )
        print(f"\n{t['hap']} {t['rsid']} ({t['aa']}) at chr22:{t['pos']}")
        if rec is None:
            print("  not found in region VCF")
            continue
        overall, perpop = tally(rec, sample_cols, sample2pop)
        ac, an = overall
        print(f"  ID in VCF: {rec[2]}  REF>ALT: {rec[3]}>{rec[4]}")
        print(f"  Overall (all samples): AC={ac} AN={an} AF={ac/an*100:.3f}%" if an else "  no calls")
        for sp in superpops:
            if sp in perpop:
                a, n = perpop[sp]
                print(f"    {sp}: AC={a} AN={n} AF={a/n*100:.3f}%")


if __name__ == "__main__":
    main()
