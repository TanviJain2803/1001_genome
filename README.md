# 1001_genome
Ecker Lab Salk's 1001 genomes project pipelines

## XGB ML Classifier for Feature Selection Pipeline
The pipeline attempts to select the genes associated with different phenotypes based on the varying transcriptomes of different accessions via Pseudobulking.   Here is the diagram representing the workflow:

<img width="563" alt="Screenshot 2023-08-16 at 4 16 45 PM" src="https://github.com/TanviJain2803/1001_genome/assets/73011639/9dccbe7f-2c7b-4539-a997-418f9aae73ce">

Some associated results from this pipeline and linear modeling can be found here: https://docs.google.com/presentation/d/1PbmYhv_nslxEZFPF7JLtSIL2ANF63-WbFdWRNgYipP4/edit?usp=sharing

<img width="594" alt="Screenshot 2023-08-16 at 4 15 29 PM" src="https://github.com/TanviJain2803/1001_genome/assets/73011639/f75b9d2f-d6c2-4136-b322-55813ba9e99b">

Pipeline also includes:
- calculating % of cell specific genes in each cluster (can print out lists of the cell specific genes and further analyse in virtual plant)

Thoughts on overall progress of workflow:
Right now we are using XGBoost for feature selection, howere there is also potential to use LASSO and maybe other models to further refine the workflow.
When uploaded to virtual plant we see the following results:

<img width="686" alt="Screenshot 2023-08-16 at 4 21 59 PM" src="https://github.com/TanviJain2803/1001_genome/assets/73011639/aded90a8-7e68-4330-96d5-72f4ed0989be">


We hope to see more of an overlap between Linear Modelling genes and more quantity in ML classifier pipeline after working on the LMM normalization and getting more data. The pipeline has shown good progress because of higher % in cell type-specific genes, which people would be more interested in. We have tried logit regression but it didn't make sense because we converted the gene count into binary numbers rather than output, and Arthur recommended against it.

## Linear Regression Pipeline
Using pseudobulk on Seurat Cluster to form gene x accession matrices for each cluster. We also import phenotype metadata (accessions and their respective phenotype val). Then using linear pearson correlation R values and its respective adj p values to filter it out.

Since scatterplots look very 0 inflated, and most heatmaps show positve r-vals, we are suspecting we need to try out a different kind of normalisation (LMM) that might improve the results. There are not a lot of cell specific genes showing up through this pipeline compared to the XGB ML model. We tried Logit regression as well, however
Flowering Time, Plant Dry Weight, and Lattitude shows good results

<img width="851" alt="Screenshot 2023-08-16 at 4 20 02 PM" src="https://github.com/TanviJain2803/1001_genome/assets/73011639/51413250-2222-4f8c-9b48-3470de6def07">

Pipeline includes:
- calculating number of genes and ecotypes in each cluster
- creating scatterplots with linear regression visualisation for any gene in any cluster (each point is an ecotype)
  notes: currently looks 0 inflated, which I suspect is why there are more positive correlations than negative correlations (maybe LMM normalisation will fix this)
- calculating % of cell specific genes in each cluster (can print out lists of the cell specific genes and further analyse in virtual plant)



