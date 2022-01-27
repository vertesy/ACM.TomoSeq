######################################################################################################
# 00.HRT.Frame.DSP.replicate.2020.01.10.R
######################################################################################################
# source("~/GitHub/Projects/HRT/Scripts.DSP.replicate/00.HRT.Frame.DSP.replicate.2020.01.10.R")

rm(list = ls(all.names = TRUE));
try(dev.off(),silent = T)

# Functions ---s---------------------
try (source('~/GitHub/Packages/CodeAndRoll/CodeAndRoll.R'),silent= F)
irequire(MarkdownReportsDev)
irequire(pheatmap)
irequire(dendextend)
irequire(openxlsx)
irequire(corrr)
irequire(ggplot2)
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("DESeq2")
irequire(DESeq2)
irequire(calibrate)

source("~/GitHub/Projects/HRT/Scripts/Functions.HRT.R")

# Setup ----------------------------
InputDir = "~/GitHub/_Papers/ACM.TomoSeq/"
InputFile = p0(InputDir, "Data/DSP.heart/DSP_heart.TranscriptCounts.tsv.gz")
OutDirOrig = OutDir = "~/Google_Drive/Avano/HRT/Analysis.DSP.2020.01/"
setup_MarkdownReports(OutDir = OutDir, scriptname = "00.HRT.Frame.DSP.replicate.2020.01.10.R")


# Metadata  ------------------------------------------------------------------------------------------------------------
MetaDataDir =  p0(InputDir, "Metadata/")
goi = read.simple.tsv(MetaDataDir, "GenesOfInterest.ARVC.txt")

# Parameters ----------------------------
{
  p=NULL
  p$'cormethod.s.clust' = "spearman"
  p$'clmethod.s.clust' = "complete"

  p$"SampleEnds" = 173
  p$"min.ERCC" = 100
  p$"k" = 4

  p$'MinGenes' =10
  p$'MinTranscripts' = 1000
  p$'MinExpr' = 5
  p$'MinExpPerSlice' = 5
  p$'MinSlicePerGene' = 3
  p$'cormethod.corrplot' = "pearson"

  # DEseq
  p$'padj' = 0.01
  p$'thr_log2fc' = 1
}

md.LogSettingsFromList(p)

# Go  ------------------------------------------------------------------------------------------------------------
HeartSlices = round(convert.tsv.data(read.simple.tsv(InputFile)), digits = 1);
HeartSlices = HeartSlices[ , 1:p$"SampleEnds"]

dim(HeartSlices)
gene_IDs_th = rownames(HeartSlices)
gene_names_th = id2name(gene_IDs_th)
slice_names = colnames(HeartSlices)


# QC & Filtering ------------------------------------------------------------------------------------------------------------
source("~/GitHub/Projects/HRT/Scripts.DSP.replicate/01.HRT.QC.DSP.replicate.R")

oo()


# Metadata ------------------------------------------------------------------------------------------------------------
CreateMetadataAnnotation = T
if (CreateMetadataAnnotation) {
  # Slices ----
  m=NULL
  m$s$slice_names = colnames(HeartSlices_HQ)
  m$s$Transcripts = colSums(HeartSlices_HQ)

  idx_M = grepv(pattern = "__chrM",x = rownames(HeartSlices_HQ), invert = F); idx_M
  m$s$Transcript.chrM = colSums(HeartSlices_HQ[idx_M, ])
  m$s$chrM.pc = iround( m$s$Transcript.chrM / m$s$Transcripts)

  # Genes ----
  m$g$gene_IDs = rownames(HeartSlices_HQ)
  m$g$Names = id2name(m$g$gene_IDs)
  m$g$Chr = id2chr(m$g$gene_IDs)

  ls_gene_IDs = list.fromNames(m$g$gene_IDs)

  m$s$ERCC.expression = HeartSlices[grepv("^ERCC-",rownames(HeartSlices[, HQ_slices])), HQ_slices]

} #if

## Clustering  ------------------------------------------------------------------------------------------------------------

Section.Clustering = T
if (Section.Clustering) {
  try.dev.off()
  slices_cor = cor(HeartSlices_norm, method = p$'cormethod.s.clust')
  MaxNonSelfCor = max(slices_cor[slices_cor<1])
  slices_cor[slices_cor==1] = MaxNonSelfCor

  PlotName = ppp("Section.Correlation",p$'cormethod.s.clust', p$'clmethod.s.clust')
  FileName =ppp("Fig.2B", PlotName)
  x = pheatmap(slices_cor, cluster_rows  = F, cluster_cols = F, main = PlotName, col = colorRamps::matlab.like(100), show_rownames = F, show_colnames = F)
  wplot_save_this(plotname = ppp(FileName, "MatlabCol"), w=9, h=8.75)

  # Clustering ------------------------------

  x = pheatmap(slices_cor, cluster_cols = F, cutree_rows = p$"k", clustering_method = p$'clmethod.s.clust', silent = T)
  # CLust stat
  m$s$'clusterID' = hclust.getClusterID.row(pheatmapObject = x,k = p$"k"); l(m$s$'clusterID')
  # m$s$'cluster.colors'= wcolorize(m$s$'clusterID', RColorBrewerSet = "Dark2", ReturnCategoriesToo = T, show = F)

  ColorFix =T
  if (ColorFix) {
    xx = translate(m$s$'clusterID', oldvalues = 1, newvalues = 5)
    yy = wcolorize(xx, RColorBrewerSet = "Dark2", ReturnCategoriesToo = T, show = F)

    # colorFix
    yy$categ["5"] = "#A9A9A9" # darkgrey
    yy$vec = translate(yy$vec, oldvalues = "#E7298A", newvalues =  "#A9A9A9")

    # nameFix
    names(yy$categ) =1:4
    yy$vec = translate(yy$vec, oldvalues = 5, newvalues = 1)
    m$s$'cluster.colors' =yy
  }
  names(m$s$'cluster.colors'$vec) = names(m$s$'clusterID')


  annot_col.create.pheatmap.vec(slices_cor, m$s$'clusterID', annot_names = "Cluster")
  annot_col$'Cluster' = m$s$'cluster.colors'$categ

  PlotName = ppp("Section.Clustering",p$'cormethod.s.clust', p$'clmethod.s.clust')
  FileName =ppp("Fig.3A", PlotName)
  x = pheatmap(slices_cor, cluster_cols = F, cutree_rows = p$"k", clustering_method = p$'clmethod.s.clust', main = PlotName, col = colorRamps::matlab.like(100)
               , annotation_col = annot,annotation_row = annot,  annotation_colors = annot_col, show_rownames = F, show_colnames = F)
  wplot_save_this(plotname = ppp(FileName, "MatlabCol"), w=7, h=6)

} #if


try.dev.off()
d1=color_branches(x$tree_row, k=p$"k", col = sort(m$s$'cluster.colors'$categ))
plot(d1,  cex=.33) # selective coloring of branches

# wcolorize()
idx_cl1=which_names(m$s$'clusterID'==1)
idx_cl2=which_names(m$s$'clusterID'==2)
idx_cl3=which_names(m$s$'clusterID'==3)
idx_cl4=which_names(m$s$'clusterID'==4)

TranscriptCounts.log10 = log10(m$s$Transcripts)-1
mRNA.content.and.chrM.expression=cbind(
  "log10(mRNA)" = TranscriptCounts.log10,
  "Expression from chrM (% total)" = 100*m$s$chrM.pc)


Fig.S5 = T
if (Fig.S5) {

  try.dev.off()
  pdfA4plot_on.layout(pname = "Fig.S5.Clustering.Results", layout_mat = rbind(1, c(2, 3),4))

  d1=color_branches(x$tree_row, k=p$"k", col = sort(m$s$'cluster.colors'$categ))
  d1 %>% set("labels_cex",.3) %>% plot()  # selective coloring of branches

  CCC = m$s$'cluster.colors'$vec[HQ_slices]

  #C
  plot.new()

  #B
  wpie(table(m$s$'clusterID'),savefile = F, col = m$s$'cluster.colors'$categ, both_pc_and_value = T, plotname = "Cluster sizes")

  #D
  wbarplot(log10(TotalTranscriptCount[HQ_slices]), col=CCC, savefile = F, ylab="log10(Transcripts)", main = "Cluster Identity and Transcript counts")
  wlegend("topleft", NamedColorVec = CCC,OverwritePrevPDF = F)

  #E
  wbarplot(mRNA.content.and.chrM.expression[HQ_slices,2], col=CCC, savefile = F, ylab="Expression from chrM (% total)", main = "Cluster Identity and chrM")
  wlegend("topleft", NamedColorVec =  CCC,OverwritePrevPDF = F)

  # OLD F and G:
  TotalTranscriptCount_per_Cl_log10 = split(log10(TotalTranscriptCount[HQ_slices]+1), m$s$'clusterID'[HQ_slices])
  wvioplot_list(TotalTranscriptCount_per_Cl_log10, ylab = "log10(Tr+1)", col = m$s$'cluster.colors'$categ, savefile = F)

  Relative.ChrM = split(100*m$s$chrM.pc, f = m$s$'clusterID'[HQ_slices])
  wvioplot_list(Relative.ChrM, ylab = "Expression from chrM (% total)", col = m$s$'cluster.colors'$categ, savefile = F)

  pdfA4plot_off()

} #if



# DiffExp ------------------------------------------------------------------------------------------------------------
DiffExp = T
if (DiffExp) source("~/GitHub/Projects/HRT/Scripts.DSP.replicate/03.Diff.Exp.HRT.DSP.replicate.R")

# Fig.2D Heatmap Similar.Geneprofiles.zScore ------------------------------------------------------------------------------------------
Fig.2D = T
if (Fig.2D) {
  Fig.2D.marker.genes = c("PLIN1__chr15", "ACTA2__chr10", "DES__chr2", "ZBTB11__chr3")
  ls_top150 = list.fromNames(Fig.2D.marker.genes)

  non.ERCCs = grepv("ERCC", HE_genes, invert = T)
  DAT = getRows(HeartSlices_zscore, non.ERCCs)
  DAT = clip.outliers(DAT,probs = c(.01, .99))

  for (i in 1:length(Fig.2D.marker.genes) ) {
    g = Fig.2D.marker.genes[i]; print(g)
    similar.genes.HRT(Expr.Matrix = DAT, gene = g, dist_ = "pearson")
    fname=ppp("Fig.2D.Similar.genes",g,"heatmap")
    wplot_save_this(plotname = fname, w=5,h=3)

    Data.S2 = TRUE
    if (Data.S2) {
      top150 = iround(sort.decreasing(similar.genes.HRT(Expr.Matrix = DAT, gene = g, plotResults = F, number = 150, dist_ = "pearson")))
      tibble::enframe(top150)

      ls_top150[[i]] <-
        tibble::enframe(top150) %>%
        tidyr::separate(col ="name", into = c("Gene.Symbol", "Chromosome"), sep = "__") %>%
        rename(value = 'Pearson CC.' )
    }

  } #for

  if (Data.S2) {
    Data.S2.150.most.similar.genes = cbind(ls_top150[[1]], ls_top150[[2]], ls_top150[[3]], ls_top150[[4]])
    hs <- createStyle(textDecoration = "BOLD", fontSize=12, fgFill = "darkolivegreen3")
    openxlsx::write.xlsx(Data.S2.150.most.similar.genes, file = "Data.S2.150.most.similar.genes.spearman.xlsx",
                         rowNames=T, firstRow=T, firstCol=T, colWidths = "auto", headerStyle = hs, tabColour ="darkgoldenrod1", creator="Vertesy") #
  }

  plot.ZBTB11.like.genes = FALSE
  if (plot.ZBTB11.like.genes) {
    # Pearson correlation based ------------------------------------------------
    # ZBTB11.like = c("DYNC1LI2__chr16", "SYNCRIP__chr6")
    pdfA4plot_on("Fig.REBUTTAL.ZBTB11.like", rows = 4, cols = 1)
    for (i in 1:length(ZBTB11.like) ) gene.barplot.HRT(gene = ZBTB11.like[i], Expr.Matrix =HeartSlices_ERCC_normalized )
    pdfA4plot_off()

  }

}



# Identify and group peak genes ------------------------------------------------------------------------------------------

# PCA.plot ------------------------------------------------------------------------------------------

Fig.4A_E.Genes.on.PCA = T
if (Fig.4A_E.Genes.on.PCA) {
  try.dev.off()

  idx.nonMito = grepv(m$g$gene_IDs,pattern =  "__chrM", invert = T)
  XX = corrr::correlate(HeartSlices_norm[idx.nonMito, ], method = p$'cormethod.corrplot')
  hrt.PCA.plot=network_plot(XX) #+ geom_point(color='darkblue')

  plotnameLastPlot = ppp("NetworkPlot.of.Cardiac.Slices", p$'cormethod.corrplot')
  plot(hrt.PCA.plot)
  wplot_save_this(plotname = plotnameLastPlot, w=20)
  try.dev.off()

  hrt.PCA.plot.OBJECT= ggplot_build(hrt.PCA.plot)
  hrt.PCA.data= hrt.PCA.plot.OBJECT$data
  XY.coord = hrt.PCA.data[[3]][,1:2]

  # Fig.3C & D  -----------------------------------------------------------------

  Fig.3B_D = T
  if (Fig.3B_D) {
    try.dev.off()
    plotLayout = rbind(1:2, 3:4, c(5,5))
    pdfA4plot_on.layout(pname = "Fig.3B_D.Section.Expression.chrM.PCA", layout_mat = plotLayout)

    plot.new()
    CCC= m$s$'cluster.colors'$vec[rownames(mRNA.content.and.chrM.expression)]
    CorCoeff = iround(c("Pearson" = cor(mRNA.content.and.chrM.expression, method = "pearson")[2],  "Spearman" = cor(mRNA.content.and.chrM.expression, method = "spearman")[2]))
    wplot(mRNA.content.and.chrM.expression, pch=21,bg = CCC, panel_first = grid())#, col=CCC, cex=1.25
    wlegend(m$s$'cluster.colors'$categ, 1, title = "Clusters")
    wlegend.label(p0(names(CorCoeff),": ",CorCoeff), title = "Corr.Coeff.:",poz = 3,  cex = .75)

    # Fig.3D  -----------------------------------------------------------------
    plot.new()
    plotname = "Fig.3D.Clusters.on.PCA"
    wplot(XY.coord, bg= m$s$'cluster.colors'$vec, plotname = plotname, panel_first = grid()
          , pch = 21, xlab="PC1", ylab="PC2", xaxt='n', yaxt='n', axes=T)
    wlegend(NamedColorVec = m$s$'cluster.colors'$categ, title = "Clusters", w =6, h = 5 )

    # Fig.3B  -----------------------------------------------------------------
    plotname = "Fig.3B.Transcript.count"
    wbarplot(as.numeric(log10(TotalTranscriptCount[HQ_slices])), col=CCC, savefile = F
             , main = ""# "Cluster Identity and Transcript counts"
             , xlab ="Epi-to-endocardium position", ylab="log10(Transcripts)")
    pdfA4plot_off()
  } #if

  # Fig.4  -----------------------------------------------------------------
  Genes.Fig.4 = c( 'PLIN1__chr15', 'FABP4__chr8', 'PNPLA2__chr11', 'CIDEC__chr3', 'ACTA2__chr10', 'VIM__chr10',
                   'COL1A2__chr7', 'FN1__chr2', 'ACTN2__chr1', 'DES__chr2', 'TTN__chr2', 'ACTA1__chr1')
  ITGB1.like.myoc = c( 'ITGB1__chr10', 'ZNF664__chr12', 'CAPZA2__chr7', 'HADHB__chr2', 'PTMA__chr2', 'NDUFB6__chr9',
                       'TMEM59__chr1', 'FSTL1__chr3', 'LOC100505738__chr5', 'CLIP1__chr12', 'SLC25A4__chr4', 'HSPA8__chr11',
                       'ATP5J__chr21', 'RPL34__chr4', 'STAT1__chr2', 'PHLDB2__chr3', 'RPS8__chr1', 'FGF12__chr3',
                       'SMYD2__chr1', 'AP1S2__chrX', 'ATAD1__chr10', 'C16orf45__chr16', 'EPAS1__chr2', 'ATP2A2__chr12')
  Genes.Fig.4 = c(Genes.Fig.4, ITGB1.like.myoc)

  plotnameLastPlot = "Fig.4A_E.Genes.on.PCA"
  pdfA4plot_on(plotnameLastPlot, rows = 6, cols = 2, w =4  )
  for (i in 1:length(Genes.Fig.4) ) {
    g  = Genes.Fig.4[i]
    ccc= wcolorize(HeartSlices_norm[g,], set = "topo.colors")
    wplot(XY.coord, bg= ccc, pch = 21, ylab="", xlab="", panel_first = F, axes=F, plotname = id2name(g))
    legend.col(col = topo.colors(100), lev = HeartSlices_norm[g,])

  } #for
  pdfA4plot_off()
}



# save.image(file = "HRT.2019.07.08.RData")


if (T) {
  defWidth = options("width")$width
  options(width = 200)
  sink(file = paste0("sessionInfo.", format(Sys.time(), format = "%Y.%m.%d"), ".txt"))
  devtools::session_info()
  sink()
  options(width = defWidth)
  rm(defWidth)
}
