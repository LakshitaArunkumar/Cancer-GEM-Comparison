#!/usr/bin/env Rscript
# Public reconstruction of comparative gene-expression module workflow.
# Source expression matrix and GEM definitions must be locally authorized.
# No original TCGA matrices or collaborator GEM tables are distributed here.
suppressPackageStartupMessages({library(GSVA);library(ggplot2)})
args <- commandArgs(trailingOnly=TRUE)
if(length(args)==0 || '--help' %in% args){
 cat('Usage: Rscript scripts/compare_gems.R --expression FILE --modules FILE --outdir DIR\n')
 quit(status=if(length(args))0 else 2)
}
if(length(args)%%2) stop('Use --key value pairs')
opts <- setNames(as.list(args[seq(2,length(args),2)]),sub('^--','',args[seq(1,length(args),2)]))
opt <- function(k){if(is.null(opts[[k]]))stop('Missing --',k);opts[[k]]}
exprfile <- opt('expression'); modfile <- opt('modules'); outdir <- opt('outdir')
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
# Input: expression CSV has column gene and sample-expression columns; values
# must already be appropriately normalized and gene names consistent.
e <- read.csv(exprfile,check.names=FALSE,stringsAsFactors=FALSE)
if (!'gene' %in% names(e))stop('expression CSV needs gene column')
if (anyDuplicated(e$gene)||anyNA(e$gene)||any(!nzchar(e$gene)))stop('Nonunique/missing gene IDs')
expr <- as.matrix(e[,-1,drop=FALSE]);rownames(expr)<-e$gene
suppressWarnings(storage.mode(expr)<-'double')
if(any(!is.finite(expr)))stop('Expression values must be finite numeric values')
if(ncol(expr)<4)stop('At least 4 samples recommended')
# Input: modules CSV with comparison,cell_type,origin,gem,gene (one gene per row).
mod <- read.csv(modfile,stringsAsFactors=FALSE)
needed <- c('comparison','cell_type','origin','gem','gene')
if(!all(needed %in% names(mod)))stop('modules CSV needs ',paste(needed,collapse=','))
mod <- mod[!is.na(mod$gene)&nzchar(mod$gene),]
if(!nrow(mod))stop('No module genes')
for (cmp in unique(mod$comparison)) for (ct in unique(mod$cell_type[mod$comparison==cmp])) {
 sub <- mod[mod$comparison==cmp & mod$cell_type==ct,]
 if(length(unique(sub$origin))!=2)next
 # Prefix GEM names by origin so module names cannot collide.
 sub$key <- paste(sub$origin,sub$gem,sep='__')
 sets <- split(sub$gene,sub$key); sets <- lapply(sets,function(z)intersect(unique(z),rownames(expr)))
 sets <- sets[lengths(sets)>=5]
 if(length(sets)<2 || length(unique(sub$origin[sub$key %in% names(sets)]))<2){
  warning('Insufficient intersected GEMs for ',cmp,' ',ct);next
 }
 # Demonstrative GSVA implementation; original analysis was exploratory.
 par <- GSVA::gsvaParam(expr,sets,minSize=5,maxSize=Inf)
 scores <- GSVA::gsva(par,verbose=FALSE)
 # Matrix is GEM x sample; similarity compares GEM score profiles ACROSS samples.
 cosine <- function(x,y){sum(x*y)/sqrt(sum(x*x)*sum(y*y))}
 origins <- unique(sub$origin)
 left <- names(sets)[startsWith(names(sets),paste0(origins[1],'__'))]
 right <- names(sets)[startsWith(names(sets),paste0(origins[2],'__'))]
 cross <- expand.grid(GEM_1=left,GEM_2=right,stringsAsFactors=FALSE)
 cross$cosine <- mapply(function(a,b)cosine(scores[a,],scores[b,]),cross$GEM_1,cross$GEM_2)
 cross$pearson <- mapply(function(a,b)cor(scores[a,],scores[b,]),cross$GEM_1,cross$GEM_2)
 cross$shared_genes <- mapply(function(a,b)length(intersect(sets[[a]],sets[[b]])),cross$GEM_1,cross$GEM_2)
 dest <- file.path(outdir,gsub('[^A-Za-z0-9_-]','_',cmp),gsub('[^A-Za-z0-9_-]','_',ct))
 dir.create(dest,recursive=TRUE,showWarnings=FALSE)
 write.csv(cross[order(-cross$cosine),],file.path(dest,'cross_cancer_similarity.csv'),row.names=FALSE)
 cross$label <- paste(cross$GEM_1,cross$GEM_2,sep=' : ')
 p <- ggplot(cross,aes(GEM_1,GEM_2,fill=cosine))+geom_tile()+
  scale_fill_gradient2(limits=c(-1,1))+labs(title=paste(cmp,ct,'GSVA cosine'),x=NULL,y=NULL)+
  theme_minimal()+theme(axis.text.x=element_text(angle=45,hjust=1))
 ggsave(file.path(dest,'cosine_heatmap.png'),p,width=9,height=6,dpi=160)
 message(cmp,' / ',ct,': ',nrow(cross),' cross-origin GEM pairs')
}
