# Targeted funding of understudied genes: the NIH IDG R03 program

## Contents
- scripts/: analysis scripts, numbered in run order
- data/raw/: downloaded source files (versions and dates below)
- data/processed/: tables made by the scripts and used for the figures

## Data sources
- IDG FOAs (RFA-RM-18-021, 19-011, 20-019, 21-012, 22-024): grants.nih.gov
- IDG R03 award list: https://druggablegenome.net/node/78
- JensenLab protein and PubTator counts, mention files (Oct 2026): https://download.jensenlab.org
- STRING v12.0 aliases: https://string-db.org
- HCOP (2026-10-01): https://www.genenames.org/tools/hcop/
- BioGRID 5.0.263: https://thebiogrid.org
- ChEMBL 37: https://www.ebi.ac.uk/chembl/
- GO goa_human.gaf (2026-05-21): https://geneontology.org
- RCSB PDB (accessed 2026-11-03): https://www.rcsb.org
- GenCC (2026-11-03): https://thegencc.org
- NIH RePORTER and ExPORTER: https://reporter.nih.gov

## Steps
1. Install Python 3.12 and matplotlib (3.6.3 or newer).
2. Run the scripts in order, 01 to 17 (lettered scripts such as 12b
   come after their number). Scripts that call NIH or NCBI services
   take your email: python <script>.py --email you@example.org

The data/raw/ files are the exact versions used in the paper. Rerunning
the download steps fetches newer versions and may change counts slightly.
