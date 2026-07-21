"""Download the public inputs used by the APOL1 frequency analyses."""

import argparse
import datetime
import json
import sys
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data" / "raw"

# Three APOL1 kidney-risk variants (GRCh38 coordinates)
VARIANTS = [
    {"rsid": "rs73885319", "variantId": "22-36265860-A-G",
     "haplotype": "G1", "label": "G1_SNP1_S342G",
     "note": "S342G missense, part of G1 risk haplotype"},
    {"rsid": "rs60910145", "variantId": "22-36265988-T-G",
     "haplotype": "G1", "label": "G1_SNP2_I384M",
     "note": "I384M missense, part of G1 risk haplotype"},
    {"rsid": "rs71785313", "variantId": None,
     "haplotype": "G2", "label": "G2_DEL_N388_Y389del",
    "note": "6 bp deletion p.N388_Y389del, the G2 risk allele"},
]

GNOMAD_API = "https://gnomad.broadinstitute.org/api"
ENSEMBL_API = "https://rest.ensembl.org"

# Canonical phase 3 panel used for population labels.
KGP_SAMPLES_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
    "integrated_call_samples_v3.20130502.ALL.panel"
)
KGP_CHR22_TBI_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "1000G_2504_high_coverage/working/20201028_3202_phased/"
    "CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz.tbi"
)

ABRAOM_SEARCH_URL = "https://abraom.ib.usp.br/search.php"
ABRAOM_AJAX_URL = "https://abraom.ib.usp.br/script.php"


def utc_now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat().replace("+00:00", "Z")

manifest = {
    "downloaded_at": utc_now(),
    "files": [],
    "errors": [],
}

UA = "apol1-frequency-analysis/1.0 (academic)"

def log_file(path, url, status, size_bytes=None, note=""):
    rel = str(path.relative_to(DATA.parent))
    manifest["files"].append({
        "path": rel,
        "url": url,
        "status": status,
        "size_bytes": size_bytes,
        "note": note,
        "fetched_at": utc_now(),
    })

def log_error(component, url, error):
    manifest["errors"].append({
        "component": component,
        "url": url,
        "error": str(error),
        "at": utc_now(),
    })

def http_get(url, timeout=60, accept=None):
    headers = {"User-Agent": UA}
    if accept:
        headers["Accept"] = accept
    req = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read(), resp.status

def http_post_json(url, payload, timeout=60):
    body = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(
        url, data=body,
        headers={"Content-Type": "application/json", "User-Agent": UA},
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read()), resp.status

def save_if_missing(path, url, fetcher, note="", force=False):
    """Save a download unless a cached copy is available."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not force:
        size = path.stat().st_size
        print(f"  [skip] {path.name}  ({size:,} bytes already on disk)")
        log_file(path, url, "cached", size, note)
        return True
    try:
        data, status = fetcher()
    except urllib.error.HTTPError as e:
        print(f"  [fail] {path.name}  HTTP {e.code} {e.reason}")
        log_error(path.name, url, f"HTTP {e.code}")
        return False
    except Exception as e:
        print(f"  [fail] {path.name}  {e}")
        log_error(path.name, url, str(e))
        return False
    if isinstance(data, str):
        data = data.encode("utf-8")
    path.write_bytes(data)
    print(f"  [ ok ] {path.name}  ({len(data):,} bytes)")
    log_file(path, url, status, len(data), note)
    return True

def fetch_gnomad(force=False):
    print("\n=== gnomAD v4.1: population-stratified frequencies ===")
    out_dir = DATA / "gnomad"
    out_dir.mkdir(parents=True, exist_ok=True)
    # Population AF is derived from AC/AN in the gnomAD v4 schema.
    query_templates = [
        ("rich-v4",
         '{ variant(rsid:"%s",dataset:gnomad_r4){variantId rsids '
         'exome{ac an af populations{id ac an homozygote_count}}'
         'genome{ac an af populations{id ac an homozygote_count}}'
         '}}'),
        ("rich-no-hom",
         '{ variant(rsid:"%s",dataset:gnomad_r4){variantId rsids '
         'exome{ac an af populations{id ac an}}'
         'genome{ac an af populations{id ac an}}'
         '}}'),
        ("medium",
         '{ variant(rsid:"%s",dataset:gnomad_r4){variantId rsids '
         'exome{ac an af}genome{ac an af}'
         '}}'),
        ("minimal",
         '{ variant(rsid:"%s",dataset:gnomad_r4){variantId rsids}}'),
    ]
    for v in VARIANTS:
        rsid = v["rsid"]
        out = out_dir / f"{rsid}_{v['haplotype']}.json"
        if out.exists() and not force:
            print(f"  [skip] {out.name}  (already on disk)")
            log_file(out, GNOMAD_API + f"#{rsid}", "cached", out.stat().st_size, v["note"])
            continue
        success = False
        for label, tmpl in query_templates:
            query = tmpl % rsid
            try:
                data, status = http_post_json(GNOMAD_API, {"query": query})
                if data.get("errors"):
                    continue
                v_data = (data.get("data") or {}).get("variant")
                if not v_data:
                    continue
                wrapped = {"_query_type": label, "_rsid": rsid, "_data": data}
                out.write_text(json.dumps(wrapped, indent=2), encoding="utf-8")
                size = out.stat().st_size
                print(f"  [ ok ] {out.name}  ({size:,} bytes, query={label})")
                log_file(out, GNOMAD_API + f"#{rsid}", "ok", size,
                         v["note"] + f" [query: {label}]")
                success = True
                break
            except Exception:
                continue
        if not success:
            print(f"  [fail] {out.name}: all gnomAD queries failed")
            log_error("gnomad", GNOMAD_API + f"#{rsid}", "all queries failed")

def fetch_ensembl(force=False):
    print("\n=== Ensembl REST: 1000 Genomes population frequencies ===")
    out_dir = DATA / "ensembl"
    out_dir.mkdir(parents=True, exist_ok=True)
    for v in VARIANTS:
        rsid = v["rsid"]
        out = out_dir / f"{rsid}_{v['haplotype']}.json"
        url = f"{ENSEMBL_API}/variation/human/{rsid}?pops=1;genotypes=0"
        def fetch():
            data, status = http_get(url, accept="application/json")
            return data, status
        save_if_missing(out, url, fetch, note=v["note"], force=force)

def fetch_kgp_meta(force=False):
    print("\n=== 1000 Genomes: sample metadata + chr22 .tbi index ===")
    out_dir = DATA / "kgp"
    out_dir.mkdir(parents=True, exist_ok=True)

    out = out_dir / "sample_metadata_phase3.txt"
    def fetch_samples():
        data, status = http_get(KGP_SAMPLES_URL, timeout=60)
        return data, status
    save_if_missing(out, KGP_SAMPLES_URL, fetch_samples,
                    note="1000G phase 3 sample and population panel", force=force)

    out_tbi = out_dir / "chr22.vcf.gz.tbi"
    def fetch_tbi():
        data, status = http_get(KGP_CHR22_TBI_URL, timeout=60)
        return data, status
    save_if_missing(out_tbi, KGP_CHR22_TBI_URL, fetch_tbi,
                     note="tabix index for 30x NYGC chr22 phased VCF", force=force)

def fetch_abraom(force=False):
    print("\n=== ABraOM: APOL1 WES variants ===")
    out_dir = DATA / "abraom"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_json = out_dir / "apol1_variants_wes_hg19.json"

    # These values reproduce the coding, PASS-quality filters in the web interface.
    body = urllib.parse.urlencode({
        "table": "abraomdb",
        "str": "APOL1",
        "gatk": "PASS",
        "cegh": "'FDP', 'FAB'",
        "exonic": "'intronic', 'downstream', 'upstream', 'UTR3', 'UTR5', 'NA'",
    }).encode("utf-8")

    def fetch():
        request = urllib.request.Request(
            ABRAOM_AJAX_URL,
            data=body,
            method="POST",
            headers={
                "Content-Type": "application/x-www-form-urlencoded",
                "User-Agent": UA,
                "Referer": ABRAOM_SEARCH_URL,
                "X-Requested-With": "XMLHttpRequest",
                "Accept": "application/json",
            },
        )
        with urllib.request.urlopen(request, timeout=90) as response:
            content = response.read()
            json.loads(content)
            return content, response.status

    save_if_missing(
        out_json,
        ABRAOM_AJAX_URL,
        fetch,
        note="APOL1 coding variants from the SABE609 exome cohort",
        force=force,
    )

def write_manifest():
    DATA.mkdir(parents=True, exist_ok=True)
    out = DATA.parent / "manifest.json"
    out.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n[manifest] {out}")

def write_readme():
    out = DATA.parent / "README.md"
    text = f"""# APOL1 analysis data

Generated by `download_data.py` on {manifest['downloaded_at']}.

## Layout

| Folder | Source | What's there |
|---|---|---|
| `raw/gnomad/` | gnomAD v4.1 GraphQL API | Population-stratified frequencies for G1 SNPs (rs73885319, rs60910145) and G2 (rs71785313) |
| `raw/ensembl/` | Ensembl REST API | Same three variants with 1000 Genomes phase 3 population frequencies |
| `raw/kgp/` | 1000 Genomes FTP | Sample-to-population panel and chr22 tabix index |
| `raw/abraom/` | ABraOM | APOL1 coding variants from the SABE609 exome cohort |
| `processed/` | Analysis scripts | Generated frequency tables |

## Provenance

`manifest.json` records each source URL, status, file size, and timestamp.
"""
    out.write_text(text, encoding="utf-8")
    print(f"[readme]   {out}")

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true",
                        help="re-download even if files already exist")
    args = parser.parse_args()
    if args.force:
        print("[force] re-downloading all files")
    fetch_gnomad(force=args.force)
    fetch_ensembl(force=args.force)
    fetch_kgp_meta(force=args.force)
    fetch_abraom(force=args.force)
    write_manifest()
    write_readme()

    print("\n" + "=" * 70)
    print("DOWNLOAD SUMMARY")
    print("=" * 70)
    n_ok = sum(1 for f in manifest["files"] if f["status"] not in ("cached",))
    n_cached = sum(1 for f in manifest["files"] if f["status"] == "cached")
    n_err = len(manifest["errors"])
    total_bytes = sum((f["size_bytes"] or 0) for f in manifest["files"])
    print(f"  Downloaded this run: {n_ok}")
    print(f"  Already on disk    : {n_cached}")
    print(f"  Errors             : {n_err}")
    print(f"  Total disk          : {total_bytes/1024:.1f} KB")
    print("=" * 70)
    if manifest["errors"]:
        print("Errors:")
        for e in manifest["errors"]:
            print(f"  - {e['component']}: {e['error'][:120]}")
    return 0 if n_err == 0 else 1

if __name__ == "__main__":
    sys.exit(main())
