# Supplementary Table S3 — Microbiota instrument coverage at the 19 colocalization CpGs

Exposure side: 456 independent MiBioGen instrument SNPs (211 taxa, genome-wide significant, p < 1e-5, clumped).
Mediator side: ARIES/mQTLdb significant cis-mQTL SNPs (p < 1e-5) at each CpG, retrieved per CpG.
A pair is testable only if the same SNP appears on both sides.

| Region | CpG | mQTL SNPs (ARIES) | Microbiota instrument SNPs in the same set | Taxa represented |
|---|---|---:|---:|---:|
| CFH CFHR 1q32 | cg11348583 | 793 | 0 | 0 |
| CFH CFHR 1q32 | cg13682187 | 752 | 0 | 0 |
| CFH CFHR 1q32 | cg07209298 | 556 | 0 | 0 |
| CFH CFHR 1q32 | cg13441112 | 0 | 0 | 0 |
| DEFA1 8p23 | cg07297068 | 85 | 0 | 0 |
| DEFA1 8p23 | cg05969298 | 37 | 0 | 0 |
| DEFA1 8p23 | cg19304352 | 22 | 0 | 0 |
| DEFA1 8p23 | cg00705004 | 0 | 0 | 0 |
| HORMAD2 22q12 | cg00539347 | 329 | 0 | 0 |
| HORMAD2 22q12 | cg27665648 | 0 | 0 | 0 |
| ITGAM ITGAX 16p11 | cg10451425 | 643 | 0 | 0 |
| ITGAM ITGAX 16p11 | cg04018474 | 610 | 0 | 0 |
| ITGAM ITGAX 16p11 | cg02846316 | 531 | 0 | 0 |
| ITGAM ITGAX 16p11 | cg03100639 | 476 | 0 | 0 |
| MADCAM1 19p13 | cg23386659 | 24 | 0 | 0 |
| MADCAM1 19p13 | cg27535410 | 0 | 0 | 0 |
| ST6GAL1 3q27 | cg17448192 | 127 | 0 | 0 |
| ST6GAL1 3q27 | cg25599673 | 1 | 0 | 0 |
| TNFRSF13B 17p11 | cg02671941 | 1 | 0 | 0 |
| **Total (19 CpGs)** | | **4987** | **0** | **0** |

## Region-level instrument availability (plus-minus 500 kb of the lead IgAN SNP)

| Region | Taxa with >=1 instrument | Instrument SNPs | Max SNPs per taxon | Smallest p |
|---|---:|---:|---:|---:|
| CFH CFHR 1q32 | 0 | 0 | 0 | NA |
| DEFA1 8p23 | 0 | 0 | 0 | NA |
| HORMAD2 22q12 | 1 | 6 | 6 | 3.331e-06 |
| ITGAM ITGAX 16p11 | 0 | 0 | 0 | NA |
| MADCAM1 19p13 | 7 | 19 | 13 | 3.74e-09 |
| ST6GAL1 3q27 | 3 | 10 | 8 | 1.675e-07 |
| TNFRSF13B 17p11 | 3 | 8 | 4 | 1.261e-06 |

Four of the seven regions carry at least one microbiota instrument, but no instrument reaches any of the 19 target CpGs.
The three remaining regions (CFH/CFHR 1q32, DEFA1 8p23, ITGAM/ITGAX 16p11) have no microbiota instrument within the window at all.

Source data: `results/g1_tool_cpg_overlap.tsv`, `results/g1_microbiota_region_summary.tsv`.
Machine-readable version of the region-level table: `tables/TableS3_microbiota_coverage.tsv` (cited in Supplementary Table S3).