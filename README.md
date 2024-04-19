# 1001_genome
Ecker Lab Salk's 1001 genomes project pipelines 


The 1001 Genomes project is a collaborative database that presents a comprehensive collection of whole-genome sequences from over 1001 accessions or ecotypes of the reference plant Arabidopsis thaliana. This data offers a spectrum of genetic and phenotypic diversity of arabidopsis, opening avenues for extensive exploration. \
At the Ecker Lab, we are building a 1001 Genome dataset to a single cell resolution. My primary research inquiry aims to explore the connections between gene expression patterns and specific phenotypes within distinct cell types of Arabidopsis. Addressing this question holds significance as it will not only enhance our comprehension of the relationship between cell types and phenotype but also has the potential to uncover crucial regulatory mechanisms within the context of plant development, adaptation, and responses to environmental stimuli.
More information on the 1001 Genomes Collaborative project can be found [here](https://1001genomes.org/)

## Lab Presentations
- Project Debrief + Clustering analysis + Linear Regression Pipeline presentation can be found [here](https://docs.google.com/presentation/d/15_TUG_j7n7F7XbB6NzNo_Mnry9L66Qx0Rx_94nymkOg/edit?usp=sharing)
- Project Debrief 2 + Wet Lab Workflow + Progress Updates + QC Optimisation Pipeline presentation can be found [here](https://docs.google.com/presentation/d/1R8Yn-lpiQqHgJAPi4vQClwj7IpXDpCbJKrpJmTsqXsc/edit?usp=sharing)
## Run Integration Pipeline 
The pipeline integrates different run count matrices and performs pre-processing steps (QC Filtering, Normalisation, Scaling, Variable Feature Selection, Linear Dimension Reduction (PCA)) and then integrates the runs.

<img width="563" alt="Screenshot 2023-08-16 at 4 16 45 PM" src="https://github.com/TanviJain2803/1001_genome/blob/main/Screenshot%202024-04-18%20at%204.45.30%20PM.png">

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
<img width="851" alt="Screenshot 2023-08-16 at 4 20 02 PM" src="https://github.com/TanviJain2803/1001_genome/blob/main/Screenshot%202024-04-18%20at%204.46.32%20PM.png">
Using pseudobulk on Seurat Cluster to form gene x accession matrices for each cluster. We also import phenotype metadata (accessions and their respective phenotype val). Then using linear pearson correlation R values and its respective adj p values to filter it out.

Since scatterplots look very 0 inflated, and most heatmaps show positve r-vals, we are suspecting we need to try out a different kind of normalisation (LMM) that might improve the results. There are not a lot of cell specific genes showing up through this pipeline compared to the XGB ML model. We tried Logit regression as well, however
Flowering Time, Plant Dry Weight, and Lattitude shows good results

<img width="851" alt="Screenshot 2023-08-16 at 4 20 02 PM" src="https://github.com/TanviJain2803/1001_genome/assets/73011639/51413250-2222-4f8c-9b48-3470de6def07">

Pipeline includes:
- calculating number of genes and ecotypes in each cluster
- creating scatterplots with linear regression visualisation for any gene in any cluster (each point is an ecotype)
  notes: currently looks 0 inflated, which I suspect is why there are more positive correlations than negative correlations (maybe LMM normalisation will fix this)
- calculating % of cell specific genes in each cluster (can print out lists of the cell specific genes and further analyse in virtual plant)

## ClusTree Pipeline
<img width="851" alt="Screenshot 2023-08-16 at 4 20 02 PM" src="https://github.com/TanviJain2803/1001_genome/blob/main/Screenshot%202024-04-18%20at%204.45.51%20PM.png">
Creating a clustering hierarchy tree using [MRTree package](https://htmlpreview.github.io/?https://github.com/pengminshi/MRtree/blob/master/vignettes/MRtree-tutorial.html). The tree uses multiple resolutions where the first resolution is 0 (all nuclei), and then a manually forced split into root and shoot (explained further below), and other findCluster() resolutions with the final resolution being 1 (leaf nodes).
The pipeline will output the following:
- New .RData seurat obj with manual included resolution saved as 'integrated_snn_res.0.02'
- UMAP .jpeg of the manual split
- ClusTree .jpeg showing hierarchy of clusters with their respective pie charts

Note about manual resolution split: We used our run2.2 nuclei and some 10x root nuclei (overlapped it). We found out what nuclei in run2.2 is root, then we imported the colnames of run2.2 nuclei into the `Genomes_Cluster` obj. We used "integrated_snn_res.1" to see what clusters have run2.2 root nuclei. If a cluster has run2.2 root nuclei, it is classified as a root cluster. The first split in the clustree is then manually forced to be split into root and shoot using the biological truth. We then verified if the latter cluster is shoot using photosynthetic markers on a heatmap:
<img width="863" alt="Screenshot 2023-10-16 at 1 58 31 PM" src="https://github.com/TanviJain2803/1001_genome/assets/73011639/556e285a-32c6-4bb4-8cb9-d8eeebec8451">
Since the heatmap gives us promising results, we are confident in the manual splitting of the initial resolution

