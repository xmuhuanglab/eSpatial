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

**Step 7. Decode the combinations of spatially divergent enhancers**

![image](https://github.com/xmuhuanglab/eSpatial/assets/95668602/6a5c6f2b-2940-4b15-bd68-9949fcfe1cb8)

## **Installation**

#### Docker setup
First install [Docker](https://docs.docker.com/get-docker/) and then pull eSpatial inside container as follows.

```
docker pull jliu5212/espatial:v1
```

## **Tutorial**
Tutorial [notebook]() describes how to run SCARlink on a small example dataset.

## **Contact:**
For any inquiries or assistance, please feel free to open an issue.
