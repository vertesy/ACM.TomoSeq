######################################################################################################
# 01.TomoHeartQC.R
######################################################################################################
# source("~/GitHub/Ideas_and_New_Projects/HRT.Scripts/2017.04.state/Functions.HRT.R")

# Functions ------------------------



gene.barplot.HRT <- function(gene="PLIN1__chr15", Expr.Matrix= HeartSlices_ERCC_normalized, color=m$s$cluster.colors$vec) {
  Labels = c(
    "HeartSlices_HQ" = "Transcripts (raw)",
    "HeartSlices_ERCC_normalized" = "Transcripts (ERCC normalized)",
    "HeartSlices_norm" = "Transcripts (sum normalized)",
    "HeartSlices_zscore" = "Transcripts (z-score)"
  )
  YaxisLabel = Labels[as.character(substitute(Expr.Matrix))];   print(YaxisLabel)
  wbarplot(Expr.Matrix[gene,], plotname = paste(gene, substitute(Expr.Matrix), collapse = '', sep = '.'),
           col = color, ylab=YaxisLabel, main = id2titlecaseitalic(gene, prefix = "Expression of"))

}



# ls_gene_IDs$""
similar.genes.HRT <- function(gene = "PLIN1__chr15", Expr.Matrix=HeartSlices_ERCC_normalized, number = 10, corThr =.9, dist_ = "pearson", plotResults=T){
  stopifnot(gene %in% rownames(Expr.Matrix))
  stopifnot(dist_ %in% c("pearson", "spearman", "euclidean","manhattan")  )
  target = as.named.vector(Expr.Matrix[gene, ], WhichDimNames = 2)

  if (dist_ == "pearson" | dist_ == "spearman") {
    corWt = as.named.vector(  cor(x = t(Expr.Matrix), y =target, method =  dist_,use = "na.or.complete" )      );
    if (number) {simGenes = tail (sort(corWt), n = number)
    } else if (corThr) {    simGenes = corWt[corWt >= corThr]; number=length(simGenes)}
    mth = "correlation"
  }
  if (dist_ == "euclidean" | dist_ == "manhattan") {
    distances = signif(apply(t(Expr.Matrix), 2, function(x) dist(rbind(target, x ), method = dist_)))
    if (number) {simGenes = head(sort(distances), n = number)
    } else if (corThr) {    simGenes = distances[distances <= corThr]; number=length(simGenes)}
    mth = "distance in gene Expr.Matrix"
  }
  print(paste("The ",number,"most similar genes based on", dist_, mth,"are: ", paste(names(simGenes), collapse = ", ")))

  if (plotResults) {
    try(dev.off, silent = T)
    DAT = Expr.Matrix[names(simGenes), ]
    # DAT = log10(DAT+1)
    rownames(DAT) = id2name(rownames(DAT))
    MinCor = if(!number) corThr else min(iround(simGenes))
    plotnameLastPlot = p0("The ", number," most similar genes to ", gene, " (",dist_,"; thr. ",MinCor,")")
    pheatmap::pheatmap(DAT, cluster_cols = F, main = plotnameLastPlot, show_colnames = F)
    assign("plotnameLastPlot", value = plotnameLastPlot, envir = .GlobalEnv)
  }
  return(simGenes)
}



# # source: http://www.gettinggeneticsdone.com/2014/05/r-volcano-plots-to-visualize-rnaseq-microarray.html
wplot_Volcano.HRT <- function(DEseqResults, thr_log2fc_ = thr_log2fc, thr_padj_ =thr_padj,saveit =T,
                               pname_ =F, highlight =T, Comp.Colors = m$s$cluster.colors$categ[CPL],
                               showNames = 0 ,CexLB=.5, Symmetric = F, genes=F, ...) {
  if (pname_ == FALSE) { pname_ = substitute(DEseqResults) }
  Columns =c("Gene", "log2FoldChange", "padj")
  DE = as.data.frame(DEseqResults)
  if (!is.null(rownames(DE)) & !"Gene" %in% colnames(DE)) {    DE = cbind( "Gene" = rownames(DE), DE)  }
  if (sum(! (Columns %in% colnames(DE)))) { any_print("A dataframe with 3 columns needed:", Columns )}
  DEseqResults = DEseqResults[!is.na(rowSums(DEseqResults)), ]

  # Make a basic volcano plot
  # subb = p0("Colors highlight genes enriched in respective clusters,\n with: p.adj< ",thr_padj_, ", and log2(FC)> ",thr_log2fc_,".")
  subb = p0("Cluster colored genes: p.adj< ",thr_padj_, " & log2(FC)> ",thr_log2fc_,".")
  logFC.4plot = DEseqResults$"log2FoldChange"
  padj.4plot = DEseqResults$"padj"
  XLM = range(logFC.4plot)
  XLM = if (Symmetric)  c(-max(abs(XLM)), max(abs(XLM)) ); print(XLM)
  with(DE, plot(x = logFC.4plot, y = -log10(padj.4plot), main=pname_, sub=subb, pch=20, cex=.5, col = rgb(0,0,0,.25), xlim=XLM, xlab="log2(fold-change)", ylab= "-log10(adjusted p-value)"))


  if (highlight) { # Add colored points:
    with(subset(DE, (padj < thr_padj_) & (log2FoldChange) < -thr_log2fc_), points(log2FoldChange, -log10(padj), pch=20, cex=.75, col=Comp.Colors[2]))
    with(subset(DE, (padj < thr_padj_) & (log2FoldChange) > thr_log2fc_), points(log2FoldChange, -log10(padj), pch=20, cex=.75, col=Comp.Colors[1]))
  }
  # Label points with the textxy function from the calibrate plot
  if (showNames) {
    if (genes == F) {
      with(subset(DE, padj<thr_padj_ & abs(log2FoldChange)>thr_log2fc_), calibrate::textxy(log2FoldChange, -log10(padj), labs=id2name(Gene), cex=CexLB))
    } else {
      with(DE[genes, ], calibrate::textxy(log2FoldChange, -log10(padj), labs=id2name(Gene), cex=CexLB))
    }
  }
  wlegend(Comp.Colors, title = "Region", OverwritePrevPDF = F)
  if (saveit) { wplot_save_this(plotname = paste0(pname_,".volcano") )  }
}



# wplot_Volcano.HRT <- function(DEseqResults, thr_log2fc_ = thr_log2fc, thr_padj_ =thr_padj,saveit =T, pname_ =F, highlight =T, Comp.Colors = m$s$cluster.colors$categ[CPL]
#                                , showNames = 0 ,CexLB=.5, Symmetric = F, ...) {
#   if (pname_ == FALSE) { pname_ = substitute(DEseqResults) }
#   Columns =c("Gene", "log2FoldChange", "padj")
#   DE = as.data.frame(DEseqResults)
#   if (!is.null(rownames(DE)) & !"Gene" %in% colnames(DE)) {    DE = cbind( "Gene" = rownames(DE), DE)  }
#   if (sum(! (Columns %in% colnames(DE)))) { any_print("A dataframe with 3 columns needed:", Columns )}
#   DEseqResults = DEseqResults[!is.na(rowSums(DEseqResults)), ]
#
#   # Make a basic volcano plot
#   # subb = p0("Colors highlight genes enriched in respective clusters,\n with: p.adj< ",thr_padj_, ", and log2(FC)> ",thr_log2fc_,".")
#   subb = p0("Cluster colored genes: p.adj< ",thr_padj_, " & log2(FC)> ",thr_log2fc_,".")
#   logFC.4plot = DEseqResults$"log2FoldChange"
#   padj.4plot = DEseqResults$"padj"
#   XLM = range(logFC.4plot)
#   XLM = if (Symmetric)  c(-max(abs(XLM)), max(abs(XLM)) ); print(XLM)
#   with(DE, plot(x = logFC.4plot, y = -log10(padj.4plot), main=pname_, sub=subb, pch=20, cex=.5, col = rgb(0,0,0,.25), xlim=XLM, xlab="log2(fold-change)", ylab= "-log10(adjusted p-value)"))
#
#
#   if (highlight) { # Add colored points:
#     with(subset(DE, (padj < thr_padj_) & (log2FoldChange) < -thr_log2fc_), points(log2FoldChange, -log10(padj), pch=20, cex=.75, col=Comp.Colors[2]))
#     with(subset(DE, (padj < thr_padj_) & (log2FoldChange) > thr_log2fc_), points(log2FoldChange, -log10(padj), pch=20, cex=.75, col=Comp.Colors[1]))
#   }
#   # Label points with the textxy function from the calibrate plot
#   if (showNames) {
#     # thr_X = rev(sort(abs(logFC.4plot)))[showNames]
#     with(subset(DE, padj<thr_padj_ & abs(log2FoldChange)>thr_log2fc_), calibrate::textxy(log2FoldChange, -log10(padj), labs=id2name(Gene), cex=CexLB))  }
#
#   wlegend(Comp.Colors, title = "Region", OverwritePrevPDF = F)
#   if (saveit) { wplot_save_this(plotname = paste0(pname_,".volcano") )  }
# }



ERCC_normalize <-function(mat = HeartSlices[, HQ_slices], pattern="^ERCC-") { # normalize each column to the sum of ERCC reads
  ERCCs = grepv(pattern = pattern,rownames(mat))
  ERCC.expression = mat[ERCCs, ]
  cs.ERCC = colSums(ERCC.expression, na.rm = T)
  norm_mat = (t(t(mat) / cs.ERCC)) * median(cs.ERCC)
  iprint("colMedians: ", head(iround(colMedians(norm_mat))))
  return(norm_mat)
}


# prepare4plotMA ------------------------------------------------------------------------
prepare4plotMA <- function(DESeq_results, thr_padj_=thr_padj, thr_log2fc_ =F) { # highlight results using 2 thresholds
  DE = as.data.frame(DESeq_results)[, c("baseMean", "log2FoldChange", "padj")]
  index_notNA = !is.na(DE$"padj") & !is.na(DE$"log2FoldChange")
  index_isSign = (DE$"padj" <= thr_padj_)
  if (thr_log2fc_ != F) {
    index_FoldChange = (na.omit.strip(DE$log2FoldChange <= -thr_log2fc_ | DE$log2FoldChange >=  thr_log2fc_))
    DE$"padj" = (index_isSign & index_FoldChange & index_notNA)
  } else { DE$"padj" = index_isSign }
  return(DE)
}
