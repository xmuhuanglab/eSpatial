# **eSpatial: a computational framework to decode the spatially divergent combinations of enhancers**
## **Introduction**

eSpatial is a computational framework to decode gene expression regulation by the spatially divergent combinations of enhancers through integrating spatial profiling of gene expression and chromatin accessibility. 

### Workflow:

Briefly, it contains the following seven steps: 

**Step 1. Prepare input matrix (Input)**

**Step 2. Reduce dimension and construct spatial neighbor network (Preprocessing)**

**Step 3. Detect spatial domains**

**Step 4. Identify spatial-specific genes and cis-regulatory elements (sGEs/sCREs)**

**Step 5. Define enhancer clusters of genes**

**Step 6. Depict spatial patterns of enhancers for** **genes expressed in multiple spatial domains**

**Step 7. Decode the spatial enhancer code**

![image](https://github.com/xmuhuanglab/eSpatial/assets/95668602/25a2a385-636a-41a1-b3e8-5137f0be3106)

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
