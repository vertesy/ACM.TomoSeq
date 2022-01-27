######################################################################################################
# 01.TomoHeartQC.R
######################################################################################################
# source("~/GitHub/Ideas_and_New_Projects/HRT.Scripts/01.TomoHeartQC.R")
try(dev.off(),silent = T)

# Functions ------------------------
# Setup ----------------------------
# Metadata  ------------------------------------------------------------------------------------------------------------
# Parameters ----------------------------
# Go  ------------------------------------------------------------------------------------------------------------


TotalTranscriptCount = colSums(HeartSlices)
# wbarplot(TotalTranscriptCount, w=14, vline = p$"SampleEnds")

TotalGeneCount = colSums(HeartSlices>=p$'MinExpr')
# wbarplot(log10(TotalGeneCount+1), vline = p$"SampleEnds")
# ComplexityPerSlice = cbind(TotalTranscriptCount,TotalGeneCount)



try.dev.off()
Fig.S2.Sequencing.quality.control = T
if (Fig.S2.Sequencing.quality.control) {


  ## Filter Slices  ------------------------------------------------------------------------------------------------------------

  HQ_slices = which_names( (TotalTranscriptCount>=p$'MinTranscripts') );l(HQ_slices)


  plotLayout = rbind(c(1:2), c(3,3), 4:5)
  pdfA4plot_on.layout(pname = "Fig.S3.Sequencing.quality.control", layout_mat = plotLayout)

  Legend = flip_value2name(c("3" = sum(TotalTranscriptCount>p$'MinTranscripts')))

  sbb = p0("Sections > ",p$'MinTranscripts' ," transcripts are retained.")
  whist(log10(TotalTranscriptCount), vline = log10(p$'MinTranscripts')
        , xlab = "log10(transcript count per section+1)", ylab="Nr. of sections", filtercol = T, breaks = 100, sub = sbb, main="", savefile = F)
  wlegend(poz = 1, title = "High Quality Sections:",NamedColorVec = Legend, OverwritePrevPDF = F)


  ## Gene Filter  ------------------------------------------------------------------------------------------------------------
  ExpressedSlices = rowSums(HeartSlices >= p$'MinExpPerSlice')
  HE_genes = which_names(ExpressedSlices> p$'MinSlicePerGene')

  pname = ""
  sbb = kollapse("Min.:", p$'MinExpPerSlice', " transcripts in Min.:", p$'MinSlicePerGene', " slices.")
  GeneFilter = iround(cbind("log10(mRNA per section+1))" = log10(rowSums(HeartSlices)), "Sections with Expression" = jitter(ExpressedSlices, amount = 5)))
  idim(GeneFilter)
  geneCol = gplots::rich.colors(8)[c(2,7)]
  ccc = geneCol[(ExpressedSlices> p$'MinSlicePerGene')+1]
  plot(GeneFilter, col=ccc, pch=20, cex=.5, # type = 'n',
       main = pname, sub=sbb)
  lll = geneCol[2]; names(lll) = l(HE_genes)
  wlegend(poz = 1, title = "Highly Expressed Genes:",NamedColorVec =   lll, OverwritePrevPDF = F)
  # plot.new()

  wbarplot(as.numeric(log10(TotalTranscriptCount+1)), hline = log10(p$'MinTranscripts'), w=14, h=7,savefile = F, ylab="log10(transcript count per section+1)", main="")
  pdfA4plot_off()

  forIllustrator = F
  if (forIllustrator) {
    pdfA4plot_on.layout(pname = "Fig.S2.gene.filter.rasterized", layout_mat = plotLayout)
    plot(GeneFilter, col=ccc, pch=20, cex=.5,
         main = pname, sub=sbb)
    lll = geneCol[2]; names(lll) = l(HE_genes)
    pdfA4plot_off()
  }

} #if


## Normalize  ------------------------------------------------------------------------------------------------------------
HeartSlices_HQ = HeartSlices[HE_genes, HQ_slices]; dim(HeartSlices_HQ)
HeartSlices_norm = median_normalize(HeartSlices_HQ)
HeartSlices_zscore = t(scale(t(HeartSlices_norm)))
HeartSlices_ERCC_normalized = ERCC_normalize(mat = HeartSlices[, HQ_slices])
Nr_HQ_slices =l(HQ_slices)


Write.out.Data.S1 = T
if (Write.out.Data.S1) {

  "The high quality sections had on average 19951 transcripts (median: 12587), corresponding to an average of 2338(median: 2141) genes per section"
  round(mean(colSums(HeartSlices[, HQ_slices])))
  round(median(colSums(HeartSlices[, HQ_slices])))

  round(mean(colSums(HeartSlices[, HQ_slices]>=1)))
  round(median(colSums(HeartSlices[, HQ_slices]>=1)))

  Data.S1 = list(
    "Raw" = HeartSlices_HQ,
    "Normalized to sum" = iround(HeartSlices_norm),
    "ERCC-Normalized" = iround(HeartSlices_ERCC_normalized),
    "Z-score" = iround(HeartSlices_zscore)
  )

  names(Data.S1)
  hs <- createStyle(textDecoration = "BOLD", fontSize=12, fgFill = "darkolivegreen3")
  openxlsx::write.xlsx(Data.S1, file = "Data.S1.Gene.Expression.Matrices.xlsx",
                       rowNames=T, firstRow=T, firstCol=T, colWidths = "auto", headerStyle = hs, tabColour ="darkgoldenrod1", creator="Vertesy") #

} #if
