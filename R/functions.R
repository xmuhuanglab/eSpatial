# Functions
# ------------------ The function GPCor is for calculating the correlation of gene-peak pairs
############## INPUT: 1. peak-cell matrix; 2. scRNA-seq matrix or gene activity matrix; 3. reference genome in Granges form
############## OUTPUT: A dataframe of peak-gene correlation
GPCor <- function(cre.mat=cre.mat,  # peak-cell matrix
                  exp.mat=rna.mat,  # scRNA-seq matrix or gene activity matrix
                  normalizeRNAMat=F, # if the rna matrix need normalization
                  normalizeATACMat = F,
                  genome = "hg19", # reference genome, must be one of "hg19", "mm10", or "hg38"
                  windowPadSize = 100000, # base pairs padded on either side of gene TSS
                  proPadSize = 2000, # base pairs padded on either side of gene TSS for enhancer
                  nCores=8 # How many registerCores to use
){
  # check genome
  if(!genome %in% c("hg19", "hg38", "mm10")){
    stop("You must specify one of hg19, hg38 or mm10 as a genome build for currently supported TSS annotations..\n")
  }
  # Normalize input matrix
  if (normalizeATACMat == T) {
    cre.mat.cc <- centerCounts(cre.mat)
  }
  cre.mat.cc <- cre.mat
  if(normalizeRNAMat==T){
    rna.mat <- NormalizeData(exp.mat, normalization.method = "LogNormalize", scale.factor = 10000) 
  }
  rna.mat = exp.mat
  # Create a se object for normalized scATAC data
  peaks <- StringToGRanges(rownames(cre.mat))
  peaks$Peak <- rownames(cre.mat)
  names(peaks) <- rownames(cre.mat)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = cre.mat.cc),
    rowRanges = peaks
  )
  
  # genePeak overlap
  Ov <- genePeakOv(ATAC.se = se,
                   RNAmat = rna.mat,
                   genome = genome,
                   windowPadSize = windowPadSize,
                   proPadSize = proPadSize)
  
  #---
  genes <- as.character(unique(Ov$gene_name))
  library(doParallel)
  getDoParRegistered()
  registerDoParallel(nCores) # registerCores
  
  GPTab <- foreach(g=genes,.combine = 'rbind',.inorder=TRUE,
                   .errorhandling = 'remove') %dopar% {
                     cat("Running gene: ",g,which(genes == g),"\n")
                     Ovd <- Ov %>% filter(gene_name == g)
                     ObsCor   <- PeakGeneCor(ATAC = cre.mat,
                                             RNA = rna.mat,
                                             peakRanges = peaks,
                                             OV = Ovd)
                     
                   }
  
  closeAllConnections() # closed cores
  return(GPTab)
}

# Nomalized data
centerCounts <- function(obj,
                         doInChunks=TRUE,
                         chunkSize=1000){
  if(!class(obj) %in% c("SummarizedExperiment","RangedSummarizedExperiment","dgCMatrix","dgeMatrix","Matrix"))
    stop("Supplied object must be either of class SummarizedExperiment or sparse Matrix ..\n")
  
  if(ncol(obj) > 10000)
    doInChunks <- TRUE
  
  if(doInChunks){
    cat("Centering counts for cells sequentially in groups of size ",
        chunkSize, " ..\n\n")
    starts <- seq(1,ncol(obj),chunkSize)
  } else{
    starts <- 1
  }
  
  counts.l <- list()
  
  for(i in 1:length(starts)){
    beginning <- starts[i]
    if(i==length(starts)) # If it's the last one
    {
      ending <- ncol(obj)
    } else {
      ending <- starts[i]+chunkSize-1
    }
    
    cat("Computing centered counts for cells: ",beginning," to ", ending,"..\n")
    
    if(class(obj) == "RangedSummarizedExperiment" | class(obj)=="SummarizedExperiment"){
      m <- SummarizedExperiment::assay(obj[, beginning:ending])} else {
        m <- obj[,beginning:ending] # Assumes Matrix format
      }
    cellMeans <- Matrix::colMeans(m)
    cat("Computing centered counts per cell using mean reads in features ..\n\n")
    # Center cell counts based on its mean RIP count
    cCounts <- Matrix::t(Matrix::t(m)/cellMeans)
    
    counts.l[[i]] <- cCounts
    
    gc()
  }
  
  cat("Merging results..\n")
  centered.counts <- do.call("cbind",counts.l)
  cat("Done!\n")
  
  if(class(obj) == "RangedSummarizedExperiment" | class(obj)=="SummarizedExperiment"){
    SummarizedExperiment::assay(obj) <- centered.counts
    return(obj)
  } else {
    return(centered.counts)
  }
}

#Gene peak overlaps
genePeakOv <- function(ATAC.se, # SummarizedExperiment object of scATAC data
                       RNAmat, # Paired normalized scRNA-seq data, with gene names as rownames
                       genome, # Must be one of "hg19", "mm10", or "hg38"
                       geneList=NULL, # 2 or more valid gene symbols (if only running on subset of genes)
                       windowPadSize=100000, # base pairs padded on either side of gene TSS
                       proPadSize = 2000 # base pairs padded on either side of gene TSS for enhancer
                       
) {
  stopifnot(inherits(ATAC.se,"RangedSummarizedExperiment"))
  stopifnot(inherits(RNAmat,c("Matrix","matrix")))
  
  if(!all.equal(ncol(ATAC.se),ncol(RNAmat)))
    stop("Input ATAC and RNA objects must have same number of cells")
  
  message("Assuming paired scATAC/scRNA-seq data ..")
  
  ATACmat <- assay(ATAC.se) # Rownames preserved
  
  if(is.null(rownames(RNAmat)))
    stop("RNA matrix must have gene names as rownames")
  
  # Check for peaks/genes with 0 accessibility/expression
  
  if(any(Matrix::rowSums(assay(ATAC.se))==0)){
    message("Peaks with 0 accessibility across cells exist ..")
    message("Removing these peaks prior to running correlations ..")
    peaksToKeep <- Matrix::rowSums(assay(ATAC.se))!=0
    ATAC.se <- ATAC.se[peaksToKeep,] # Subset ranges
    ATACmat <- ATACmat[peaksToKeep,]
    message("Important: peak indices in returned gene-peak maps are relative to original input SE")
  }
  
  
  peakRanges <- granges(ATAC.se) # Peak ranges
  names(peakRanges) <- rownames(ATAC.se)
  if(any(Matrix::rowSums(RNAmat)==0)){
    message("Genes with 0 expression across cells exist ..")
    message("Removing these genes prior to running correlations ..")
    genesToKeep <- Matrix::rowSums(RNAmat)!=0
    RNAmat <- RNAmat[genesToKeep,]
  }
  
  cat("Number of peaks in ATAC data:",nrow(ATACmat),"\n")
  cat("Number of genes in RNA data:",nrow(RNAmat),"\n")
  
  
  if (!genome %in% c("hg19", "hg38", "mm10"))
    stop("You must specify one of hg19, hg38 or mm10 as a genome build for currently supported TSS annotations..\n")
  # switch(genome, hg19 = {
  #   TSSg <- hg19TSSRanges
  # }, hg38 = {
  #   TSSg <- hg38TSSRanges
  # }, mm10 = {
  #   TSSg <- mm10TSSRanges
  # })
  if(genome %in% c("hg19", "hg38", "mm10")){
    #load(paste0('../data/', genome, '_refSeq.Rdata'))
    TSSg <- get(paste0(genome, 'TSSRanges'))
  }
  
  # Keep genes that have annotation and are in RNA matrix
  names(TSSg) <- as.character(TSSg$gene_name)
  
  if(!is.null(geneList)){
    if(length(geneList)==1)
      stop("Please specify more than 1 valid gene symbol")
    
    if(any(!geneList %in% names(TSSg))){
      cat("One or more of the gene names supplied is not present in the TSS annotation specified: \n")
      cat(geneList[!geneList %in% names(TSSg)], sep = ", ")
      cat("\n")
      stop()
    }
    
    TSSg <- TSSg[geneList]
  }
  
  # Get peak summit
  cat("\nTaking peak summits from peak windows ..\n")
  peakSummits <- resize(peakRanges,width = 1,fix = "center")
  
  # Pad promoter by this much *either side*
  proflank <- GenomicRanges::flank(TSSg,
                                   width = proPadSize,
                                   both = TRUE)
  # Find overlap of all peaks to all genes given window
  # Subject is Peaks, query is Gene
  cat("Finding overlapping promoter-gene pairs ..\n")
  geneProOv <- findOverlaps(query = proflank,subject = peakSummits)
  numpgPairs <- length(geneProOv)
  
  cat("Found ",numpgPairs,"total promoter-peak pairs for given TSS window ..\n")
  
  cat("Number of peak summits that overlap any gene promoter window: ",length(unique(subjectHits(geneProOv))),"\n")
  cat("Number of gene promoter windows that overlap any peak summit: ",length(unique(queryHits(geneProOv))),"\n\n")
  peakSummits <- peakSummits[-unique(subjectHits(geneProOv))]
  
  # Checking in case some genes in RNA don't overlap our TSS annotations
  genesToKeep <- intersect(names(TSSg),rownames(RNAmat))
  cat("\nNum genes overlapping TSS annotation and RNA matrix being considered: ",length(genesToKeep),"\n")
  
  # Match gene order in RNA matrix and TSS ranges
  RNAmat <- RNAmat[genesToKeep,]
  TSSg <- TSSg[genesToKeep]
  
  
  # Pad TSS by this much *either side*
  TSSflank <- GenomicRanges::flank(TSSg,
                                   width = windowPadSize,
                                   both = TRUE)
  # Find overlap of all peaks to all genes given window
  # Subject is Peaks, query is Gene
  cat("Finding overlapping peak-gene pairs ..\n")
  genePeakOv <- findOverlapPairs(query = TSSflank,subject = peakSummits)
  genePeakOv <- data.frame(genePeakOv@first,genePeakOv@second)
  numPairs <- nrow(genePeakOv)
  
  cat("Found ",numPairs,"total gene-peak pairs for given TSS window ..\n")
  
  cat("Number of peak summits that overlap any gene TSS window: ",length(unique(genePeakOv$Peak)),"\n")
  cat("Number of gene TSS windows that overlap any peak summit: ",length(unique(genePeakOv$gene_name)),"\n\n")
  
  return(genePeakOv)
  
}

# Gene peak correlation
PeakGeneCor <- function(ATAC, # Normalized reads in peaks counts (rownames should  be "Peak1","Peak2" etc.)
                        RNA, # Normalized gene expression counts
                        OV, # Gene TSS - Peak overlap pairs object (Genes: query, Peaks: subject)
                        peakRanges,
                        seed = 2022,
                        mtd="spearman") {
  
  
  geneIndices <- as.character(OV$gene_name)
  peakIndices <- as.character(OV$Peak)
  
  uniquegenes <- unique(geneIndices)
  uniquepeaks <- unique(peakIndices)
  
  
  idy = uniquepeaks
  ovlpn = length(idy)
  idyexcl = which(!1:length(peakRanges) %in% idy)
  gene.name = uniquegenes
  
  gene.val = RNA[uniquegenes,];
  
  data.raw = ATAC
  peak.raw = peakRanges
  # use set
  data.use = ATAC[idy,];
  peak.use = peak.raw[uniquepeaks];
  
  # shuf set
  data.shuf <- data.use
  peak.shuf <- peak.use
  colheader <- colnames(data.shuf)
  #Shuffle row-wise:
  set.seed(seed)
  data.shuf <- data.shuf[sample(nrow(data.shuf)),]
  #Shuffle column-wise:
  set.seed(seed)
  data.shuf <- data.shuf[,sample(ncol(data.shuf))]
  colnames(data.shuf) <- colheader
  
  # rdm set
  set.seed(seed)
  
  #    rdm <- sample(idyexcl, 1000)
  rdm <- sample(idyexcl, ovlpn)
  data.rdm = data.raw[rdm,];
  peak.rdm = peak.raw[rdm];
  
  
  # rdm shuf set
  data.rdmshuf <- data.rdm
  peak.rdmshuf <- peak.rdm
  colheader <- colnames(data.rdmshuf)
  #Shuffle row-wise:
  set.seed(seed)
  data.rdmshuf <- data.rdmshuf[sample(nrow(data.rdmshuf)),]
  #Shuffle column-wise:
  set.seed(seed)
  data.rdmshuf <- data.rdmshuf[,sample(ncol(data.rdmshuf))]
  colnames(data.rdmshuf) <- colheader
  
  corrFunc <- function(var1, var2, method) {
    result = cor.test(var1, var2, method = method)
    data.frame(result[c("estimate","p.value","statistic")],
               stringsAsFactors=FALSE)
  }
  
  #-----------------------------------
  #cat("Real set: performing statitical test on... \n", file = stderr())
  peaks.id = seq(nrow(data.use));
  corr = lapply(peaks.id, function(t){
    corrFunc(as.numeric(gene.val), as.numeric(data.use[t,]), mtd)
  })
  corr <- do.call(rbind, corr)
  
  peak.use$estimate <- corr[, "estimate"]
  peak.use$statistic <- corr[, "statistic"]
  peak.use$method <- mtd
  peak.use$Pval <- corr[, "p.value"]
  peak.use$FDR <- p.adjust(peak.use$Pval, method = "BH")
  peak.use$class <- "corr"
  peak.use$Gene <- gene.name
  
  #cat("Done ... \n", file = stderr())
  
  #------------------------------
  #cat("Shuffle set: performing statitical test on... \n", file = stderr())
  peaks.id = seq(nrow(data.shuf));
  corr = lapply(peaks.id, function(t){
    corrFunc(as.numeric(gene.val), as.numeric(data.shuf[t,]), mtd)
  })
  corr <- do.call(rbind, corr)
  
  peak.shuf$estimate <- corr[, "estimate"]
  peak.shuf$statistic <- corr[, "statistic"]
  peak.shuf$method <- mtd
  peak.shuf$Pval <- corr[, "p.value"]
  peak.shuf$FDR <- p.adjust(peak.shuf$Pval, method = "BH")
  peak.shuf$class <- "shuf"
  peak.shuf$Gene <- gene.name
  
  #cat("Done ... \n", file = stderr())
  
  #--------------------------------
  #cat("Random set: performing statitical test on... \n", file = stderr())
  peaks.id = seq(nrow(data.rdm));
  corr = lapply(peaks.id, function(t){
    corrFunc(as.numeric(gene.val), as.numeric(data.rdm[t,]), mtd)
  })
  corr <- do.call(rbind, corr)
  
  peak.rdm$estimate <- corr[, "estimate"]
  peak.rdm$statistic <- corr[, "statistic"]
  peak.rdm$method <- mtd
  peak.rdm$Pval <- corr[, "p.value"]
  peak.rdm$FDR <- p.adjust(peak.rdm$Pval, method = "BH")
  peak.rdm$class <- "random"
  peak.rdm$Gene <- gene.name
  #cat("Done ... \n", file = stderr())
  
  
  #----------------------
  #cat("Rdmshuf set: performing statitical test on... \n", file = stderr())
  peaks.id = seq(nrow(data.rdmshuf));
  corr = lapply(peaks.id, function(t){
    corrFunc(as.numeric(gene.val), as.numeric(data.rdmshuf[t,]), mtd)
  })
  corr <- do.call(rbind, corr)
  peak.rdmshuf$estimate <- corr[, "estimate"]
  peak.rdmshuf$statistic <- corr[, "statistic"]
  peak.rdmshuf$method <- mtd
  peak.rdmshuf$Pval <- corr[, "p.value"]
  peak.rdmshuf$FDR <- p.adjust(peak.rdmshuf$Pval, method = "BH")
  peak.rdmshuf$class <- "rdmShuf"
  peak.rdmshuf$Gene <- gene.name
  #cat("Done ... \n", file = stderr())
  
  #############################
  # OUTPUT
  peak.use.df <- as.data.frame(peak.use)
  peak.shuf.df <- as.data.frame(peak.shuf)
  peak.rdm.df <- as.data.frame(peak.rdm)
  peak.rdmshuf.df <- as.data.frame(peak.rdmshuf)
  peak.df <- rbind(peak.use.df, peak.shuf.df)
  peak.df <- rbind(peak.df, peak.rdm.df)
  peak.df <- rbind(peak.df, peak.rdmshuf.df)
  
  return(peak.df);
  
}



SpecificityScore <- function(object = obj, # Seurat object with cluster information in the metadata
                             assays = "RNA", # Assays to calculate specificity score
                             clusters = NA # Clusters to calculate specificity score
                             ) {
  DefaultAssay(object) <- assays
  calculateFoldChange <- function(cluster) {
    avg_fc <- FoldChange.Seurat(obj, ident.1 = cluster)
    return(avg_fc$avg_log2FC)
  }
  
  if (assays %in% c("RNA", "peaks")) {
    object <- switch(
      assays,
      RNA = NormalizeData(object),
      peaks = RunTFIDF(object)
    )
    
    features <- rownames(object)
    mat <- data.frame(matrix(nrow = length(features), ncol = length(clusters)), row.names = features)
    colnames(mat) <- clusters
    
    mat[, clusters] <- t(lapply(clusters, calculateFoldChange))
    
    return(mat)
  } else {
    stop("Invalid assay type. Supported types are 'RNA' and 'peaks'.")
  }
}


geneModules = function(exp.mat = exp.mat, # Normalized gene expression matrix
                       k = 6, # Number of k-means
                       clusters = colnames(exp.mat) 
                       ){
  exp.Z <- t(scale(t(exp.mat), center = TRUE))
  keep.genes <- names(which(rowMeans(exp.Z) != 0))
  exp.km <- kmeans(exp.Z[keep.genes, ], k, nstart = 25)
  cluster <- exp.km$cluster
  kmeansInfo <- data.frame(
    row.names = names(cluster),
    genes = names(cluster),
    kmeans = paste0("k", cluster)
  )
  # order kmeans
  kDomain <- lapply(seq_along(1:k), function(i) {
    tmp <- cluster[cluster == i]
    genes <- names(tmp)
    means <- colMeans(exp.mat[genes, ])
    df <- matrix(data = means, nrow = 1)
  }) %>% do.call(what = rbind)
  
  rownames(kDomain) <- paste0("k", 1:k)
  colnames(kDomain) <- colnames(exp.mat)
  kDomain <- kDomain[,clusters]
  
  rank <- apply(kDomain, 1, function(row) {
    paste(order(row, decreasing = TRUE), collapse = "")
  })
  sorted_clusters <- rownames(kDomain)[order(rank)]
  
  kDF <- lapply(seq_along(1:k), function(i) {
    tmp <- kmeansInfo %>% filter(kmeans == sorted_clusters[i])
    tmp$kmeans = paste0("k",i)
    return(tmp)
  }) %>% do.call(what = rbind)
  kDF$kmeans = factor(kDF$kmeans,levels = paste0("k",1:k))
  return(kDF)
}

enhancerOrder = function(exp.mat = exp.mat, # Normalized gene expression matrix
                         cre.mat = cre.mat, # Normalized chromatin accessibility matrix
                         genemodule = genemodule,# Gene module information
                         gppairs = GPTabFilt, # Gene enhancer pairs 
                         clusters = clusters,
                         cutoff = 0.02 # A threshold to identify the clusters each gene module expressed
                         ){
  # Calculate gene module expression in each domain 
  kmeans <- as.character(unique(genemodule$kmeans))
  kDomain <- lapply(kmeans, function(i) {
    genes <- genemodule$genes[which(genemodule$kmeans == i)]
    genes <- genemodule %>% filter(kmeans == i) %>% pull(genes)
    means <- colMeans(exp.mat[genes, ])
    df <- matrix(data = means, nrow = 1)
  }) %>% do.call(what = rbind)
  
  rownames(kDomain) <- kmeans
  colnames(kDomain) <- colnames(exp.mat)
  kDomain <- kDomain[,clusters]
  
  peakod <- data.frame()
  for (module in kmeans) {
    domain <- kDomain[module,]
    domain <- names(domain[domain > cutoff])
    

    kgene <- genemodule %>% filter(kmeans == module) %>% pull(genes)
    kgene <- intersect(kgene, rownames(exp.mat))
    kpeak <- GPTabFilt %>% filter(Gene %in% kgene) %>% pull(Peak)
    kpeak <- intersect(rownames(cre.mat), kpeak)
    

    cre.sub <- cre.mat[kpeak, domain]
    
    if (length(domain) > 1) {
   
      peakmodule <- data.frame(peak = kpeak, module = NA, score = NA, row.names = kpeak, genemodule = module)
      peakmodule$module <- colnames(cre.sub)[max.col(cre.sub[kpeak,])]
      peakmodule$score <- unlist(lapply(kpeak, function(x) cre.sub[x, which.max(cre.sub[x, ])]))
      
      pod <- lapply(domain, function(m1) {
        peakmodule %>% filter(module == m1) %>% arrange(desc(score))
      }) %>% bind_rows()
      peakod <- bind_rows(peakod, pod)
    }
    
    if (length(domain) == 1) {
      pod <- data.frame(peak = names(cre.sub), module = domain, score = cre.sub, genemodule = module)
      pod <- pod %>% arrange(desc(score))
      peakod <- bind_rows(peakod, pod)
    }
  }
  return(peakod)
}

enhancerPattern = function(exp.mat = exp.mat, # Normalized gene expression matrix
                           cre.mat = cre.mat, # Normalized chromatin accessibility matrix
                           genemodule = gene.module, # Gene module information
                           gppairs = GPTabFilt, # Gene enhancer pairs 
                           clusters = clusters,
                           cutoff.cluster = 0.02, # A threshold to identify the clusters each gene module expressed
                           cutoff.bin = 0 # A threshold to binary the specificity score
                           ) {
  
  # Calculate gene module expression in each domain 
  kmeans <- as.character(unique(genemodule$kmeans))
  kDomain <- lapply(kmeans, function(i) {
    genes <- genemodule$genes[genemodule$kmeans == i]
    means <- colMeans(exp.mat[genes, ])
    matrix(data = means, nrow = 1)
  }) %>% do.call(what = rbind)
  
  rownames(kDomain) <- kmeans
  colnames(kDomain) <- colnames(exp.mat)
  kDomain <- kDomain[, clusters]
  
  peakod <- data.frame()
  cutoff.bin <- mean(cre.mat) + sd(cre.mat)
  for (m in seq_along(kmeans)) {
    module <- kmeans[m]
    domain <- kDomain[module, ]
    domain <- names(domain[domain > cutoff.cluster])
    
    kgene <- genemodule %>% filter(kmeans == module) %>% pull(genes)
    kgene <- intersect(kgene, rownames(exp.mat))
    
    #---peak module
    kpeak <- GPTabFilt %>% filter(Gene %in% kgene) %>% pull(Peak)
    kpeak <- intersect(rownames(cre.mat), kpeak)
    cre.sub <- cre.mat[kpeak, domain]
    
    cre.bin <- ifelse(cre.sub >= cutoff.bin, 1, 0)
    
    if (length(domain) > 1) {
      #---gene module
      comb <- lapply(kpeak, function(x) {
        com <- paste0(cre.bin[x, ], collapse = "")
        data.frame(
          peaks = x,
          expDomainNum = length(domain),
          geneModule = module,
          type = com,
          enhDomainNum = sum(cre.bin[x, ])
        )
      }) %>% do.call(what = rbind)
      
      #---filtered random combination
      ncom <- round(nrow(comb) * 0.005)
      comInfo <- data.frame(table(comb$type))
      radom <- comInfo$Var1[comInfo$Freq < ncom]
      comb$type[comb$type %in% radom] <- "0-radom"
      combType <- unique(comb$type)
      combType <- combType[order(combType, decreasing = TRUE)]
      
      pod <- lapply(seq_along(combType), function(x) {
        pp <- comb %>% filter(type == combType[x])
        pp$type <- paste0("M", x)
        return(pp)
      }) %>% do.call(what = rbind)
    }
    
    if (length(domain) == 1) {
      comb <- unique(cre.bin)
      pod <- data.frame(
        peaks = names(cre.bin),
        expDomainNum = length(domain),
        geneModule = paste0("m", m),
        type = cre.bin,
        enhDomainNum = cre.bin
      )
      pod <- lapply(comb, function(x) {
        pp <- pod %>% filter(type == x)
      }) %>% do.call(what = rbind)
    }
    peakod <- rbind(peakod, pod)
  }
  #colnames(peakod) <- c("Enhancer","exp.Domain","gene.Module","enh.Pattern","enh.Domain")
  return(peakod)
}
enhancerCombination <- function(enh.pattern = enh.pattern, # Enhancer pattern information
                                gppairs = GPTabFilt # Gene enhancer pairs 
                                ) {
  genes <- gppairs %>% 
    filter(Peak %in% enh.pattern$peaks) %>% 
    pull(Gene) %>% 
    unique()
  
  mylapply <- ifelse(test = TRUE, yes = pblapply, no = lapply)
  
  god <- mylapply(X = genes, FUN = function(x) {
    enh <- gppairs %>% filter(Gene %in% x) %>% pull(Peak)
    tmp <- enh.pattern %>% filter(peaks %in% enh)
    
    geneInfo <- data.frame(
      gene = x,
      enhNum = length(enh),
      expNum = unique(tmp$expDomainNum),
      comNum = length(unique(tmp$type))
    )
  }) %>% do.call(what = rbind)
  colnames(god) <- c("Gene","enh.Num","exp.Domain","enh.Pattern.Num")
  return(god)
}

plotHeatmap = function(matrix = mat,
                       scale = TRUE,
                       row_split = genemodule$kmeans,
                       col_split = NULL,
                       col = c("#352A86", "#343DAE", "#0262E0", "#1389D2", 
                               "#2DB7A3", "#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D"),
                       name = "expression"){
  if(scale){
    matrix_hm <- t(apply(matrix, 1, scale))
    rownames(matrix_hm) <- rownames(matrix)
    colnames(matrix_hm) <- colnames(matrix)
  }else{
    matrix_hm <- matrix
  }
  
  p <- Heatmap(
    matrix_hm,
    name = name,
    show_column_names = TRUE,
    show_row_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_names_rot = 90,
    col = col,
    row_split = row_split,
    row_title_rot = 0,
    height = 10,
    border = TRUE,
    border_gp = gpar(lex = 1)
  )
  return(p)
}


FoldChange.Seurat <- function(
  object,
  ident.1 = NULL,
  ident.2 = NULL,
  group.by = NULL,
  subset.ident = NULL,
  assay = NULL,
  slot = 'data',
  reduction = NULL,
  features = NULL,
  pseudocount.use = 1,
  mean.fxn = NULL,
  base = 2,
  fc.name = NULL,
  ...
) {
  if (!is.null(x = group.by)) {
    if (!is.null(x = subset.ident)) {
      object <- subset(x = object, idents = subset.ident)
    }
    Idents(object = object) <- group.by
  }
  if (!is.null(x = assay) && !is.null(x = reduction)) {
    stop("Please only specify either assay or reduction.")
  }
  # select which data to use
  if (is.null(x = reduction)) {
    assay <- assay %||% DefaultAssay(object = object)
    data.use <- object[[assay]]
    cellnames.use <-  colnames(x = data.use)
  } else {
    data.use <- object[[reduction]]
    cellnames.use <- rownames(data.use)
  }
  cells <- IdentsToCells(
    object = object,
    ident.1 = ident.1,
    ident.2 = ident.2,
    cellnames.use = cellnames.use
  )
  # check normalization method
  norm.command <- paste0("NormalizeData.", assay)
  norm.method <- if (norm.command %in% Command(object = object) && is.null(x = reduction)) {
    Command(
      object = object,
      command = norm.command,
      value = "normalization.method"
    )
  } else if (length(x = intersect(x = c("FindIntegrationAnchors", "FindTransferAnchors"), y = Command(object = object)))) {
    command <- intersect(x = c("FindIntegrationAnchors", "FindTransferAnchors"), y = Command(object = object))[1]
    Command(
      object = object,
      command = command,
      value = "normalization.method"
    )
  } else {
    NULL
  }
  fc.results <- FoldChange(
    object = data.use,
    cells.1 = cells$cells.1,
    cells.2 = cells$cells.2,
    features = features,
    slot = slot,
    pseudocount.use = pseudocount.use,
    mean.fxn = mean.fxn,
    base = base,
    fc.name = fc.name,
    norm.method = norm.method
  )
  return(fc.results)
}

IdentsToCells <- function(
  object,
  ident.1,
  ident.2,
  cellnames.use
) {
  #
  if (is.null(x = ident.1)) {
    stop("Please provide ident.1")
  } else if ((length(x = ident.1) == 1 && ident.1[1] == 'clustertree') || is(object = ident.1, class2 = 'phylo')) {
    if (is.null(x = ident.2)) {
      stop("Please pass a node to 'ident.2' to run FindMarkers on a tree")
    }
    tree <- if (is(object = ident.1, class2 = 'phylo')) {
      ident.1
    } else {
      Tool(object = object, slot = 'BuildClusterTree')
    }
    if (is.null(x = tree)) {
      stop("Please run 'BuildClusterTree' or pass an object of class 'phylo' as 'ident.1'")
    }
    ident.1 <- tree$tip.label[GetLeftDescendants(tree = tree, node = ident.2)]
    ident.2 <- tree$tip.label[GetRightDescendants(tree = tree, node = ident.2)]
  }
  if (length(x = as.vector(x = ident.1)) > 1 &&
      any(as.character(x = ident.1) %in% cellnames.use)) {
    bad.cells <- cellnames.use[which(x = !as.character(x = ident.1) %in% cellnames.use)]
    if (length(x = bad.cells) > 0) {
      stop(paste0("The following cell names provided to ident.1 are not present in the object: ", paste(bad.cells, collapse = ", ")))
    }
  } else {
    ident.1 <- WhichCells(object = object, idents = ident.1)
  }
  # if NULL for ident.2, use all other cells
  if (length(x = as.vector(x = ident.2)) > 1 &&
      any(as.character(x = ident.2) %in% cellnames.use)) {
    bad.cells <- cellnames.use[which(!as.character(x = ident.2) %in% cellnames.use)]
    if (length(x = bad.cells) > 0) {
      stop(paste0("The following cell names provided to ident.2 are not present in the object: ", paste(bad.cells, collapse = ", ")))
    }
  } else {
    if (is.null(x = ident.2)) {
      ident.2 <- setdiff(x = cellnames.use, y = ident.1)
    } else {
      ident.2 <- WhichCells(object = object, idents = ident.2)
    }
  }
  return(list(cells.1 = ident.1, cells.2 = ident.2))
}
