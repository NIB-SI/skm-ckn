# O'Mally (2016) DAP-seq data

Last updated: 28.12.2024

## Description



## Citations: 

Ronan C. O’Malley, Shao-shan Carol Huang, Liang Song, Mathew G. Lewsey, Anna Bartlett, Joseph R. Nery, Mary Galli, Andrea Gallavotti, and Joseph R. Ecker (2016). Cistrome and Epicistrome Features Shape the Regulatory DNA Landscape, Cell. 2016 May 19;165(5):1280-92. doi: 10.1016/j.cell.2016.04.038. 

Lavrekha, V.V., Levitsky, V.G., Tsukanov, A.V., Bogomolov, A.G., Grigorovich, D.A., Omelyanchuk, N., Ubogoeva, E.V., Zemlyanskaya, E.V. and Mironova, V. (2022). CisCross: A gene list enrichment analysis to predict upstream regulators in Arabidopsis thaliana. Frontiers in Plant Science, 13, p.942710.

## ETL Methods

### Raw data access

__Original data deposit:__  https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE60141

__CisCross peak calling:__  https://plamorph.sysbio.ru/ciscross/AboutCisCross.html

__Arabidopsis gene models:__  https://www.arabidopsis.org/download/list?dir=Genes%2FAraport11_genome_release

__EPD The Eukaryotic Promoter Database:__ https://epd.expasy.org/ftp/epdnew/A_thaliana/current/At_EPDnew.sga

### Transformation

#### 1. Fetch sample metadata

    01_fetch_sample_metadata_and_data.ipynb

#### 2. Prepare dapseq data

    02_prepare-dapseq-data.ipynb

#### 3. Prepare target data

    03_prepare-target-data.ipynb


Input: reference_files/prepare_araport11_gff.ipynb

#### 4. Extract tf-target pairs

    04_extract-tf-targets.ipynb

## Additional notes

-