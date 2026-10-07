# GoMapMan 2 (MapMan4) annotation of Arabidopsis genes

Last updated: 06.10.2026

## Description

MapMan4 bins of the Arabidopsis genes from GoMapMan 2 (https://gomapman2.nib.si), for CKN v2.0.1 (`mapman` node
column, `scripts/ckn_v2.0.1.py` at the repository root) and for the PSS functional clusters (`mapman` property: the
bins of a cluster's Arabidopsis genes, see `skm-pss/pss-updates/2026.09.28-data-fixes.cypher` block 18 and
`regular-updates.cypher`).

## ETL Methods

### Raw data access

From the GoMapMan 2 exports page (downloaded 2026-10-06), in `raw_data/` (not in git):

- `MapMan4_gene-binCode_pairs.zip`: Identifier, plantName, binCode (13 plants; Arabidopsis: 32168 pairs)
- `MapMan4_ontology.zip`: binCode, binName, ... (8751 bins)

`raw_data/archive/` has the files of an earlier, classic MapMan exploration (`scripts/gomapman.ipynb`).

### Transformation

`python scripts/ath_mapman4.py` writes `transformed_data/ath_mapman4.tsv`: per gene (AGI, upper case), its bins,
`|`-joined (as in CKN and the PSS downloads), each `<bin code>_<full bin name>` (the names of the bin and its parents
joined with `.`; bin names have no `.`, `;` or `|`). Bin 35 (not assigned / not annotated) is left out. 16170 genes,
20683 gene-bin pairs. Copy it into the Neo4j import folder (`skm-neo4j/import`) for the PSS updates.

## Additional notes

Compared with CKN v2's GMM column (classic MapMan): 15379 nodes have MapMan4 bins, against 17607 with real classic
bins (and 5810 more with only bin 35): 3962 protein-coding genes are not assigned in MapMan4, and the MapMan4 export has
no metabolites (57 had classic bins).
