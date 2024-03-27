# Dataset name

Last updated: dd.mm.yyyy

## Description

https://dianalab.e-ce.uth.gr/tarbasev9

## Citation: 

Giorgos Skoufos, Panos Kakoulidis, Spyros Tastsoglou, Elissavet Zacharopoulou, Vasiliki Kotsira, Marios Miliotis, Galatea Mavromati, Dimitris Grigoriadis, Maria Zioga, Angeliki Velli, Ioanna Koutou, Dimitra Karagkouni, Steve Stavropoulos, Filippos S Kardaras, Anna Lifousi, Eustathia Vavalou, Armen Ovsepian, Anargyros Skoulakis, Sotiris K Tasoulis, Spiros V Georgakopoulos, Vassilis P Plagianakos, Artemis G Hatzigeorgiou, TarBase-v9.0 extends experimentally supported miRNA–gene interactions to cell-types and virally encoded miRNAs, Nucleic Acids Research, 2023, DOI:  https://doi.org/10.1093/nar/gkad1071

## ETL Methods

### Raw data access

https://dianalab.e-ce.uth.gr/tarbasev9/data/Other_species_TarBase-v9.tsv.gz

The headers and the descriptions of the fields in the interactions files are provided below :

species: The name of the species.
mirna_name: The miRNA name according to miRBase-v22.
mirna_id: The miRNA id according to miRBase-v22.
gene_name: The gene name according to Ensembl 109.
gene_id: The gene id according to Ensembl 109.
gene_location: The region inside the gene body in which the miRNA binds (e.g., 3UTR or CDS).
transcript_name: The transcript name according to Ensembl 109.
transcript_id: The transcript id according to Ensembl 109.
chromosome: The chromosome from which the gene originates.
start: The start (x) coordinate of the MRE/Binding site of the miRNA in the gene (Ensembl 109 – we are using an absolute coordinate system relative to the entire genome).
end: The end (y) coordinate of the MRE/Binding site of the miRNA in the gene (Ensembl 109 – we are using an absolute coordinate system relative to the entire genome).
strand: The strand from which the gene originates (+ indicates forward strand while – the opposite).
experimental_method: The experimental method from which the interaction was derived.
regulation: The regulation type of the interaction (i.e., whether the miRNA down- or up-regulates the target gene – Negative means downregulation while Positive the opposite).
tissue: The tissue that was utilized to conduct the experiment.
cell_line: The cell line that was utilized to conduct the experiment.
article_pubmed_id: PubMed id of the article from which the interaction was derived from.
confidence: In case of biological replicates, a value of 1 means that the same interaction appeared in at least two of the replicates.
interaction_group: This field is a way to distinguish the best (in terms of microCLIP score) interaction from the rest of the interactions that appear in the same MRE.
cell_type: The cell-type that was utilized for the experiment (e.g., Fibroblasts, Macrophages etc.).
microt_score: The predicted interaction score of microT-CDS-2023 for the same, experimentally supported interaction. The score of microT-CDS-2023 is always between 0 and 1.
comment: Free text that further annotates the interaction.

Note: The string “NA” (i.e., Not Available) is used to denote missing values.


### Transformation


## Additional notes

