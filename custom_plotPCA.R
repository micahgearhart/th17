plotPCA.DESeqTransform = function(object, intgroup="condition", ntop=500, dims=c(1,2),returnData=FALSE)
{
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    colData(object)[[intgroup]]
  }
  
  # assembly the data for the plot
  d <- data.frame(a=pca$x[,dims[1]], b=pca$x[,dims[2]], group=group, intgroup.df, name=colnames(object))
  colnames(d)[1:2]<-paste0("PC",dims)
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[dims]
    attr(d, "rotation")<-pca$rotation
    return(d)
  }
  
  ggplot(data=d, aes_string(x=paste0("PC",dims[1]), y=paste0("PC",dims[2]), color=intgroup[1],shape=intgroup[2])) + geom_point(size=3,stroke = 2) + 
    xlab(paste0("PC",dims[1],": ",round(percentVar[dims[1]] * 100),"% variance")) +
    ylab(paste0("PC",dims[2],": ",round(percentVar[dims[2]] * 100),"% variance")) 
    # removed coord_fixed()
}

g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 
