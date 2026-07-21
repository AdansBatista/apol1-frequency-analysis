# APOL1 G1/G2 Frequency Analysis

**Author:** Adans Schmidt Batista (Student ID: 3774163)\
**Course:** BIOL625, Bioinformatics and Genomics (Rev. 1)\
**Assignment:** Assignment 3, Final Project: Revised Submission

This repository reproduces APOL1 G1 and G2 risk-allele frequencies from gnomAD v4,
1000 Genomes, and the Brazilian ABraOM SABE609 exome cohort. It also compares named
West African populations and validates the 1000 Genomes aggregate results from a
Galaxy-extracted regional VCF.

G1 is represented by rs73885319 (S342G) and rs60910145 (I384M). G2 is rs71785313
(N388_Y389del). ABraOM stores G2 under the former identifier rs143830837, which
dbSNP merged into rs71785313.

## Scripts

| Script | Role |
|---|---|
| `download_data.py` | Downloads public gnomAD, Ensembl/1000 Genomes, and ABraOM inputs |
| `analyze_apol1_frequencies.py` | Calculates continental and Brazilian frequencies with exact confidence intervals |
| `galaxy_region_frequencies.py` | Tallies G1 and G2 from a Galaxy-extracted 1000 Genomes VCF |
| `analyze_west_african_frequencies.py` | Calculates detailed West African and diaspora frequencies |
| `make_variant_alignment_figure.py` | Draws the G1 substitutions and G2 protein deletion |
| `make_movement_map.py` | Draws the West Africa-to-Americas historical context map |

## Setup

Python 3.12 or newer is recommended.

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

## Run order

```powershell
python download_data.py
python analyze_apol1_frequencies.py
python analyze_west_african_frequencies.py
python make_variant_alignment_figure.py
python make_movement_map.py
```

Use `--force` with `download_data.py` to replace cached inputs.

The Galaxy validation requires the APOL1 regional VCF and the 1000 Genomes sample
panel downloaded above:

```powershell
python galaxy_region_frequencies.py path/to/apol1_region.vcf
```

## Outputs

Raw inputs and the provenance manifest are written under `data/`. Processed tables are
written to `data/processed/`, and figures are written to `figures/`.

## Main reproduced values

| Group | G1 frequency | G2 frequency |
|---|---:|---:|
| gnomAD African/African American | 22.87% | 13.58% |
| 1000 Genomes African super-population | 25.95% | 12.86% |
| ABraOM SABE609, Sao Paulo | 2.87% | 1.97% |

The named West African analysis shows wider regional variation than the continental
averages. G1 ranges from 12.35% in MSL to 49.49% in ESN, while G2 ranges from
7.87% in YRI to 19.47% in GWD.

## Code authorship

These scripts are original project code prepared with disclosed large language model
assistance. I reviewed, executed, and verified the code and its outputs.
