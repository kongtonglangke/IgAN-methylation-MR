# External data sources

No raw third-party data are redistributed in this repository. Every input is publicly available
from the sources below. Access dates were taken from the timestamps of the locally retrieved
files, which span **30 October to 1 November 2026**.

| Dataset | Role in the study | Source | Accession / endpoint | Accessed |
|---|---|---|---|---|
| GoDMC blood mQTL | methylation QTL instruments | mQTL database (mQTLdb) http://www.mqtldb.org/ | GoDMC release | 2026-10-30 |
| ARIES mQTL | developmental-window mQTL, replication arm | mQTLdb http://www.mqtldb.org/ | ARIES, five windows | 2026-10-30 |
| FinnGen R12 | IgAN outcome GWAS | https://www.finngen.fi/ | endpoint N14_IGA_NEPHROPATHY (822 cases, 499,526 controls) | 2026-10-30 |
| NHGRI-EBI GWAS Catalog | IgAN association records for instrument-first pairing | https://www.ebi.ac.uk/gwas/ | MONDO_0005342 (183 records, 182 variants) | 2026-10-30 |
| MiBioGen | gut microbiota instruments | https://mibiogen.gcc.rug.nl/ | 211 taxa, pruned at r2 < 0.001 within 1 Mb (456 variants) | 2026-10-30 |
| Gd-IgA1 GWAS effects | O-glycosylation pathway arm | Kiryluk et al., PLoS Genet 2017;13(2):e1006609 (PMID 28187132) | sentinel variants only; no genome-wide summary statistics are public | 2026-10-30 |
| GSE35069 | sorted blood cell methylation | NCBI GEO https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35069 | 60 samples, 10 fractions | 2026-11-01 |
| ENCODE cCRE | candidate cis-regulatory elements | UCSC track `encodeCcreCombined` (hg38) via https://api.genome.ucsc.edu/ | accessed through the UCSC API because the ENCODE Portal returned HTTP 405 to programmatic access | 2026-11-01 |
| ChIP-Atlas ATAC | chromatin accessibility | https://chip-atlas.org/ | `ATC.ALL.05.AllAg.AllCell` (hg38), tabix-indexed, read remotely | 2026-11-01 |
| 1000 Genomes phase 3 | structural variants and LD panel | https://www.internationalgenome.org/ | phase 3 merged SV call set; EUR panel | 2026-10-30 |
| hg19 to hg38 chain | coordinate lift-over | UCSC https://hgdownload.soe.ucsc.edu/ | `hg19ToHg38.over.chain.gz` | 2026-10-30 |

## Notes

- mQTLdb and ARIES release only CpG–SNP pairs passing a per-cohort significance threshold.
  This is a structural property of the resource and is discussed in the paper (Sections 2.1, 3.8, 4.4, 4.5).
- The ENCODE Portal HTTP 405 workaround and the ChIP-Atlas remote-tabix workflow are documented in
  `scripts/g4_ccre.py` and `scripts/g4_atac_direct.py`.
