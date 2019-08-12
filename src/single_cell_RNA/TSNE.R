function(
  object
  cells.use=NULL
  dims.use=1:10
  k.seed=1
  do.fast=T
  add.iter=0
  genes.use=NULL
  reduction.use="pca"
  dim_embed=2
  
  cells.use=Seurat:::set.ifnull(cells.use,colnames(object@data))
  if (is.null(genes.use)) {
    dim.code=Seurat:::translate.dim.code(reduction.use); dim.codes=paste(dim.code,dims.use,sep="")
    data.use=FetchData(object,dim.codes)
  }
  if (!is.null(genes.use)) {
    genes.use=ainb(genes.use,rownames(object@data))
    data.use=t(object@data[genes.use,cells.use])
  }
  
  
  #data.dist=as.dist(mahalanobis.dist(data.use))
  if (do.fast) {
    set.seed(k.seed); data.tsne=Rtsne::Rtsne(as.matrix(data.use),dims=dim_embed,...)
    data.tsne=data.frame(data.tsne$Y)
  }
  if (!(do.fast)) {
    set.seed(k.seed); data.tsne=data.frame(tsne(data.use,k=dim_embed,...))
  }
  if (add.iter > 0) {
    data.tsne=data.frame(tsne(data.use,initial_config = as.matrix(data.tsne),max_iter = add.iter,...))
  }
  colnames(data.tsne)=paste("tSNE_",1:ncol(data.tsne),sep="")
  rownames(data.tsne)=cells.use
  object@tsne.rot=data.tsne
  return(object)
}