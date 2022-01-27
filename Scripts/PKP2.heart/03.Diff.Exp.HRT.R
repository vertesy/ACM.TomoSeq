######################################################################################################
# DifferentialExpression ZH
######################################################################################################
# source("~/GitHub/ACM.TomoSeq/Scripts/03.Diff.Exp.HRT.R")
try(dev.off(), silent = T)

# Setup --------------------------------------------------------------------------------------------------------------------------------------
create_set_SubDir("03.Diff.Exp")
source('~/GitHub/TheCorvinas/R/DatabaseLinke.R')


# Paramters --------------------------------------------------------------------------------------------------------------------------------------
Pairwise.DE =T

# DE calculation and plotting --------------------------------------------------------------------------------------------------------------------------------------
SampleNames = colnames(HeartSlices_HQ)
clusterID =m$s$Clusters
NrCl = l(unique(clusterID))

# Pairwise DEseq comparison --------------------------------------------------------
if (Pairwise.DE) {
  comparisons.numeric = list("c2vs1" = c(2,1),   "c3vs1" = c(3,1),   "c2vs3" = c(2,3) )
  ls_clust = splititsnames_byValues(m$s$clusterID)

  # create a lists of sections to be compared against each other
  ls_idx = list.fromNames(names(comparisons.numeric))
  for (i in 1:3) ls_idx[[i]] <- list.2.replicated.name.vec(ls_clust[comparisons.numeric[[i]]])
  comparisons = names(ls_idx)

  # create a list of 0/1 comparison tables for DESEQ
  set1stvalueTo0_2ndTo1 <- function(x) t(t(translate(vec = x, oldvalues = unique(x), newvalues = c(0,1))))
  ls_coldata = lapply(ls_idx, set1stvalueTo0_2ndTo1 )

  count__data = round(HeartSlices_HQ)

  DE_objects = paste0("DE_",comparisons)
  DE_results = paste0("res_",comparisons)
  DE_objects_ls = list.fromNames(DE_objects);
  DE_results_ls = list.fromNames(DE_results)

  i=1
  for (i in 1:length(comparisons) ) {
    name=DE_objects[i]; print(name)

    DAT = count__data[ , names(ls_idx[[i]]) ]
    COLDAT = ls_coldata[[i]]; colnames(COLDAT) = comparisons[i]

    DDSo =DESeqDataSetFromMatrix("countData" = DAT, "colData" = COLDAT, "design" = eval(parse(text = paste("~",comparisons[i]) ) ) )
    DE_objects_ls[[name]] = DESeq(DDSo, parallel = T, betaPrior = TRUE )

    name2=DE_results[i]
    res = iround(as.data.frame(results(DE_objects_ls[[name]], contrast=c(comparisons[i],"0", "1")) ))

    res =  res[order(res$'log2FoldChange'), ]
    Hugos = id2name(rownames(res))

    # Add links
    DE_results_ls[[name2]] = cbind( res,
                                    "HGNC" = link_HGNC(Hugos, writeOut = F, Open = F),
                                    "uniprot_human" = link_uniprot_human(Hugos, writeOut = F, Open = F),
                                    "ensembl" = link_ensembl.grc37(Hugos, writeOut = F, Open = F),
                                    "String" = link_String(Hugos, organism = "human", writeOut = F, Open = F),
                                    "pubmed" = link_pubmed(Hugos, additional_terms = " human", writeOut = F, Open = F)
    )
    # write.simple.tsv(DE_results_ls[[name2]], ManualName = paste0(name2, ".tsv"))
  } #for



  Data.S3 = DE_results_ls
  names(Data.S3) = paste("Cluster", substr(comparisons,2,7))
  hs <- createStyle(textDecoration = "BOLD", fontSize=12, fgFill = "darkolivegreen3")
  openxlsx::write.xlsx(Data.S3, file = "Data.S3.Pairwise.Differential.Gene.Expression.xlsx",
                       rowNames=T, firstRow=T, firstCol=T, colWidths = "auto", headerStyle = hs, tabColour ="darkgoldenrod1", creator="Vertesy") #



  iprint("Diff. expression calculations ready")

  # Plot results, with accepted hits higlighted ---------------------------------
  COLZ = c( 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj', 'HGNC', 'uniprot_human', 'ensembl', 'String', 'pubmed')
  COLZ.num = c( 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj')

  All.in.1.pdf = T
  if (All.in.1.pdf) { try.dev.off()
    pname = "Fig.5A.S4A-B.Differential.Gene.Expression";
    pdfA4plot_on(pname = pname, cols = 2, rows = 3)

    for (i in 1:length(comparisons) ) {
      name2=DE_results[i]
      rs = DE_results_ls[[name2]][ , COLZ.num ]
      pname = comparisons[i]
      plotMA(prepare4plotMA(rs, thr_padj_ = p$'padj', thr_log2fc_ = p$'thr_log2fc'), main=pname)
      ccc = m$s$cluster.colors$categ[comparisons.numeric[[i]]]
      wplot_Volcano.HRT(rs, thr_padj_ = p$'padj', thr_log2fc_ = p$'thr_log2fc', pname_ = pname, Comp.Colors = ccc,saveit = F, Symmetric = T)
    }

    for (i in 1:length(comparisons) ) {
      name2=DE_results[i]
      rs = DE_results_ls[[name2]][ , COLZ.num ]
      pname = comparisons[i]
      ccc= m$s$cluster.colors$categ[comparisons.numeric[[i]]]
      wplot_Volcano.HRT(rs, thr_padj_ = p$'padj', thr_log2fc_ = 0, pname_ = pname, Comp.Colors = ccc,saveit = F, Symmetric = T)
      wplot_Volcano.HRT(rs, thr_padj_ = p$'padj', thr_log2fc_ = p$'thr_log2fc', pname_ = pname, Comp.Colors = ccc,saveit = F, Symmetric = T, showNames = T)
    }
    pdfA4plot_off()
  }


  Fig.5A.and.S4AB.DE = T
  if (Fig.5A.and.S4AB.DE) { try.dev.off()
    pname = "Fig.S5.D-E.Pairwise.Differerntial.Expression.Analysis";
    pdfA4plot_on(pname = pname, cols = 2, rows = 3)
i=1
    for (i in 1:length(comparisons) ) {
      name2=DE_results[i]
      rs = DE_results_ls[[name2]][ , COLZ.num ]
      pname = comparisons[i]
      # plot.new()
      ccc= m$s$cluster.colors$categ[comparisons.numeric[[i]]]
      highlight.ls = list(
        c('FABP4', 'ACTA2','PLIN1', 'PLIN4','VIM', 'TTN', 'MYL2', 'TNNI3', 'MYH7', 'ACTN2', 'CKM', 'ACTN2', 'NPPA'),
        c('FABP4', 'ACTA2','PLIN1', 'PLIN4','VIM', 'TTN', 'MYL2', 'TNNI3', 'MYH7', 'ACTN2', 'CKM', 'ACTN2', 'NPPA'),
        c('ZBTB11', 'PLIN1',  'PLIN4', 'RNF11',  'DES',  'CKM',  'AK1',  'SRL')
      )
      setdiff( highlight.ls[[i]], id2name(rownames(rs)))
      idx.plt = which(id2name(rownames(rs)) %in% highlight.ls[[i]])
      wplot_Volcano.HRT(rs, thr_padj_ = p$'padj', thr_log2fc_ = p$'thr_log2fc', pname_ = pname, Comp.Colors = ccc,saveit = F, Symmetric = T, showNames = T, genes = rownames(rs)[idx.plt] )
      wplot_Volcano.HRT(rs, thr_padj_ = p$'padj', thr_log2fc_ = p$'thr_log2fc', pname_ = pname, Comp.Colors = ccc,saveit = F, Symmetric = T)
    }
    pdfA4plot_off()
  }
  iprint("Plots are ready")
} # Pairwise.DE



# Fig.S5.Expression.profiles -----------------------------------------------------------------------
Fig.5.and.S5.Expression.profiles = T
if (Fig.5.and.S5.Expression.profiles) {
  try.dev.off()
  calculate.CDF = F
  if (calculate.CDF) {
    ls_exp_per_region = colsplit(HeartSlices_HQ, m$s$'clusterID'); lapply(ls_exp_per_region, idim)
    ls_MedianExpr_per_region = lapply(ls_exp_per_region, rowSums)
    ls_MedianExpr_per_region = lapply(ls_MedianExpr_per_region, sort, decreasing=T)
    total.expr.FFR = ls_MedianExpr_per_region[[1]]
  } #if

  # ls_gene_IDs$DES__chr2
  Genes.Fig.5 = c("DES__chr2","PLIN1__chr15", "ZBTB11__chr3", "ADH1B__chr4", "RNF11__chr1")
  Genes.Fig.5.Suppl = c("MYH7__chr14", "CKM__chr19","ADH1B__chr4", "PLIN4__chr19")

  pdfA4plot_on("Fig.5B_E.Expression.profiles", rows = 4, cols = 1)
  for (i in 1:length(Genes.Fig.5) ) gene.barplot.HRT(gene = Genes.Fig.5[i], Expr.Matrix =HeartSlices_ERCC_normalized )
  pdfA4plot_off()

  pdfA4plot_on("Fig.S5.Additional.expression.profiles", rows = 4, cols = 1)
  for (i in 1:length(Genes.Fig.5.Suppl) ) gene.barplot.HRT(gene = Genes.Fig.5.Suppl[i], Expr.Matrix =HeartSlices_ERCC_normalized )
  pdfA4plot_off()

} #if


create_set_Original_OutDir(NewOutDir = ParentDir)
