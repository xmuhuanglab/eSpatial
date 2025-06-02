# **eSpatial: a computational framework to decode the spatially divergent combinations of enhancers**
## **Introduction**

eSpatial is a framework to decipher the spatial regulation of enhancer clusters controlling the same gene based on spatial chromatin accessibility profiling.

### Workflow:

 Briefly, eSpatial analysis comprises the following seven key steps:

**Step 1. Prepare input matrix (Input)**

**Step 2. Reduce dimension and construct spatial neighbor network (Preprocessing)**

**Step 3. Define cell types**

**Step 4. Detect spatial domains**

**Step 5. Define enhancer clusters of genes**

**Step 6. Identify spatial enhancer units based on spatial patterns of enhancers**

**Step 7. Decode the combinations of divergent enhancer units**

![image](https://github.com/user-attachments/assets/bf4d5bf3-e6f2-4301-b490-fd24609b598a)


## **Installation**

#### Docker setup
First install [Docker](https://docs.docker.com/get-docker/) and then pull eSpatial inside container as follows.
```
docker pull jliu5212/espatial:v1
```

## **Tutorial**
Tutorial [notebook](https://github.com/xmuhuanglab/eSpatial/tree/main/notebooks) describes how to run eSpatial on a small example dataset.

## **Contact:**
For any inquiries or assistance, please feel free to open an issue.

## **Reference**
Hong, D., Shu, M., Liu, J. et al. Divergent combinations of enhancers encode spatial gene expression. Nat Commun 16, 5091 (2025). [https://doi.org/10.1038/s41467-025-60482-1](https://doi.org/10.1038/s41467-025-60482-1)
