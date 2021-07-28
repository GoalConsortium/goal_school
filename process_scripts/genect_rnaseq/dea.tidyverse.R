#!/cm/shared/apps/R/intel/3.2.1/bin/Rscript
library(edgeR)
library(DESeq2)
library(tidyverse)
library("RColorBrewer")
library(qusage)
library("heatmaply")

gradient_col <- ggplot2::scale_fill_gradient2(
  low = "blue", high = "red", 
  midpoint = 0.5, limits = c(0, 1)
)

rowMax <- function(x) apply(x,1,max)

#################### Read in Data ################################
tbl <- read_delim('countTable.txt',delim="\t")
tbl2 <- read_delim('countTable.logCPM.txt',delim="\t")
dtbl <- read_delim('design.txt',delim="\t")
genenames <- select(tbl,ENSEMBL,SYMBOL,TYPE)

ct <- as.data.frame(select(tbl,!(ENSEMBL:TYPE)))
row.names(ct) <- tbl$ENSEMBL
libSizes <- as.vector(colSums(ct))
logcpm <- as.data.frame(select(tbl2,!(ENSEMBL:TYPE)))
row.names(logcpm) <- tbl2$ENSEMBL

colData <- dtbl %>% mutate(SampleID= factor(SampleID,levels=colnames(ct))) %>% arrange(SampleID)
grpnames <- levels(factor(as.character(colData$SampleGroup)))
colData$SampleGroup <- factor(colData$SampleGroup, levels=grpnames)

if (length(libSizes > 1000000) < 1) {
  print(paste("Samples are filtered if there is < 1M reads.  There are less than no remaining sample(s) after this filter",sep=' '))
  q()
}

dds <- DESeqDataSetFromMatrix(countData=ct,colData= colData,design= ~ SampleGroup)
dds <- dds[ rowMax(counts(dds)) > 30, ]
dds <- dds[ colSums(counts(dds)) > 1000000]

tmp.tab <- aggregate(t(logcpm),by=list(as.character(colData$SampleGroup)),FUN=mean)
row.names(tmp.tab) <- tmp.tab$Group.1
mean.by.group <- as.data.frame(round(log2(t(tmp.tab[,c(2:ncol(tmp.tab))])+1), digits = 2))
mean.by.group$ENSEMBL <- row.names(mean.by.group)
#################### Run Sample Comparison ################################
dds <- DESeq(dds)
rld <- rlogTransformation(dds, blind=TRUE)
sampleDists <- as.data.frame(as.matrix(dist(t(assay(rld)))))
sdata <- colData %>% mutate(SampleID= factor(SampleID,levels=colnames(sampleDists))) %>% arrange(SampleID)

png(file="samples_heatmap.png",bg ="transparent",height=768,width=1024)
ggheatmap(as.matrix(sampleDists),col_side_colors = colData$SampleGroup,seriate = "mean", show_dendrogram = c(FALSE, TRUE))
dev.off()

#Compare Samples using PCA
png(file="pca.png",bg ="transparent",height=768,width=1024)
print(plotPCA(rld, intgroup="SampleGroup"),col.hab=col.blocks)
dev.off()
#################### Run DESEQ2 ################################
#Do all pairwise comparisons
contrast <- resultsNames(dds)
cond<- levels(colData(dds)$SampleGroup)
a <- length(cond)-1
for (i in 1:a) {
  for (j in 2:length(cond)) {
    if (i == j) {
      next
    } else {
      res <- as.data.frame(results(dds,contrast=c("SampleGroup",cond[i],cond[j])))
      res$ENSEMBL <- row.names(res)
      output <- inner_join(inner_join(genenames,res),mean.by.group)
      output <- mutate(output, rawP = pvalue, logFC = log2FoldChange, fdr = padj, bonf = p.adjust(rawP, method ='bonferroni'))
      write_delim(output,file=paste(cond[i],'_',cond[j],'.deseq2.txt',sep=""),quote=FALSE,delim='\t')
      filt.out <- filter(output,fdr < 0.05)    
      if (nrow(filt.out) > 2) {
        s <- filter(tbl2, ENSEMBL %in% filt.out$ENSEMBL)
        subset <- as.data.frame(select(s,!(ENSEMBL:TYPE)))
        row.names(subset) <- s$SYMBOL
        ngenes <- nrow(subset)
        textscale <- (1/(ngenes/30))
        if (textscale > 1) {
          textscale <-1
        }
        if (textscale < 0.1) {
          textscale <- 0.1
        }
        png(file=paste(cond[i],'_',cond[j],'.heatmap.deseq2.png',sep=""),height=768,width=1024)
        ggheatmap(subset[rowMax(subset) > 1,],col_side_colors = colData$SampleGroup,row_dend_left = TRUE,scale='row',cexRow=textscale)
        dev.off()
      }
    }
  }
}

###################### Run EdgeR ########################
###################### Run QuSage ######################
run.qusage <- 0
if (file.exists('geneset.gmt')) {
  MSIG.geneSets <- read.gmt('geneset.gmt')
  run.qusage <- 1
}

design <- model.matrix(~colData$SampleGroup)
d <- DGEList(counts=ct,group=colData$SampleGroup,lib.size=libSizes)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
png(file="mds.png",bg ="transparent",height=768,width=1024)
col.blocks <- as.numeric(colData$SampleGroup)
plotMDS(d, labels=colData$SampleGroup,col=col.blocks)
dev.off()

cond <-levels(d$samples$group)
colnames(design) <- levels(d$samples$group)
a <- length(cond)-1
for (i in 1:a) {
  for (j in 2:length(cond)) {
    if (i == j) {
      next
    } else {
      c <- exactTest(d, pair=c(cond[j],cond[i]))
      res <- c$table
      res$ENSEMBL <- row.names(res)
      output <- inner_join(inner_join(genenames,res),mean.by.group)
      output <- mutate(output, rawP = PValue, fdr = p.adjust(output$PValue, method ='fdr'), bonf = p.adjust(PValue, method ='bonferroni'))
      write_delim(output,file=paste(cond[i],'_',cond[j],'.deseq2.txt',sep=""),quote=FALSE,delim='\t')
      filt.out <- filter(output,fdr < 0.05)    
      if (nrow(filt.out) > 2) {
        s <- filter(tbl2, ENSEMBL %in% filt.out$ENSEMBL)
        subset <- as.data.frame(select(s,!(ENSEMBL:TYPE)))
        row.names(subset) <- s$SYMBOL
        ngenes <- nrow(subset)
        textscale <- (1/(ngenes/30))
        if (textscale > 1) {
          textscale <-1
        }
        if (textscale < 0.1) {
          textscale <- 0.1
        }
        png(file=paste(cond[i],'_',cond[j],'.heatmap.deseq2.png',sep=""),height=768,width=1024)
        ggheatmap(subset[rowMax(subset) > 1,],col_side_colors = colData$SampleGroup,row_dend_left = TRUE,scale='row',cexRow=textscale)
        dev.off()
        if (run.qusage > 0) {
          gcont <- paste(cond[j],cond[i],sep='-')
          qs.results = qusage(logcpm, grps,gcont,MSIG.geneSets)
          save(qs.results,file=paste(cond[i],'_',cond[j],'.qusage.rda',sep=""))
        }
      }
    }
  }
}

###################### Run Limma VOOM ########################

# design <- model.matrix(~0+grps)
# colnames(design) <- grpnames
# a <- length(cond)-1
# design.pairs <- c()
# k <- 0
# for (i in 1:a) {
#   for (j in 2:length(cond)) {
#     if (i == j) {
#       next
#     } else {
#     k <- k+1
#     design.pairs[k] <- paste(cond[i],'-',cond[j],sep='')
#     }
#     }
#     }

#contrast.matrix <- makeContrasts(design.pairs,levels=design)

#d <- DGEList(counts=countTable,group=grps,lib.size=libSizes)
#d <- calcNormFactors(d)
#d <- estimateCommonDisp(d)

# v <- voom(d,design,plot=TRUE)
# fit <- lmFit(v,design)
# fit2 <- contrasts.fit(fit, contrast.matrix)
# fit2 <- eBayes(fit2)

# comps <- colnames(fit2$coefficients)
# for (i in 1:length(comps)) {
#   res <- topTable(fit2,coef=i,number=Inf,sort.by="P")
#   res2 <- merge(genenames,res,by.x='ensembl',by.y='row.names',all.y=TRUE,all.x=FALSE)
#   output <- merge(res2,mean.by.group,by.y="row.names",by.x='symbol')
#   output$rawP <- output$P.Value
#   output$logFC <- output$adj.P.Val
#   output$fdr <- p.adjust(output$rawP, method ='fdr')
#   output$bonf <- p.adjust(output$rawP, method ='bonferroni')
#   write.table(output,file=paste(cond[i],'_',cond[j],'.voom.txt',sep=""),quote=FALSE,row.names=FALSE,sep='\t')
#   filt.out <- na.omit(output[output$fdr < 0.05,])
#   if (nrow(filt.out) > 2) {
#   subset <- logcpm[row.names(logcpm) %in% filt.out$ensembl,]
#   gnames <- filt.out[c('ensembl','symbol')]
#   s <- merge(gnames,subset,by.x="ensembl",by.y="row.names",all.x=FALSE,all.y=TRUE,sort=FALSE)
#   STREE <- hclust(dist(t(subset)))
#   zscores <- scale(t(subset))
#   ngenes <- length(colnames(zscores))
#   textscale <- (1/(ngenes/30))
#   if (textscale > 1) {
#      textscale <-1
#   }
#   if (textscale < 0.1) {
#     textscale <- 0.1
#   }
#   png(file=paste(cond[i],'_',cond[j],'.heatmap.voom.png',sep=""),height=768,width=1024)
#   heatmap.2(zscores, col = bluered(100),Rowv = as.dendrogram(STREE), RowSideColors = col.blocks,dendrogram='row', cexCol=textscale,labCol=s$symbol,srtRow=45,srtCol=45,trace="none", margins=c(5, 5))
#   legend("topright",legend=grpnames,col=rainbow(length(grpnames)),pch=20,cex=0.5)
#   dev.off()
#   }
# }

