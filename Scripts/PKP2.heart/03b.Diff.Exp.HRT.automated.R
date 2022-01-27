




# DEseq comparison --------------------------------------------------------


# comparisons.numeric = list("c1vs2" = c(1,2),   "c1vs3" = c(1,3),   "c2vs3" = c(2,3) )
comparisons.numeric = list("c2vs1" = c(2,1),   "c3vs1" = c(3,1),   "c2vs3" = c(2,3) )
ls_clust = splititsnames_byValues(m$s$clusterID)

# create a lists of sections to be compared against each other
ls_idx = list.fromNames(names(comparisons.numeric))
for (i in 1:3) ls_idx[[i]] <- list.2.replicated.name.vec(ls_clust[comparisons.numeric[[i]]])
comparisons = names(ls_idx)

# create a list of 0/1 comparison tables for DESEQ
make.01.char.vec <- function(x) as.character(as.factor.numeric(x)-1)
minusmin <- function(x) t(t(make.01.char.vec(x))) # creates a 0/1 vector, fits it into a dataframe column. This makes it 0/1: as.character(as.numeric(as.logical(X)))
ls_coldata = lapply(ls_idx, minusmin )


count__data = round(HeartSlices_HQ)


DE_objects = paste0("DE_",comparisons)
DE_results = paste0("res_",comparisons)
DE_objects_ls = list.fromNames(DE_objects);
DE_results_ls = list.fromNames(DE_results)


i=3
for (i in 1:length(comparisons) ) {
  name=DE_objects[i]; print(name)
  if (F) {
    count__data2 = round(HeartSlices[ , names(ls_idx[[i]] )]); idim(count__data2) # subset dataset
    idx_row_EXP = rowSums(count__data2)>1 # Some genes are not expressed in these slices, remove them
    table(idx_row_EXP)
    DAT = count__data2[idx_row_EXP, ]

  } else {
    DAT = count__data[,names(ls_idx[[i]] )]

  }
  COLDAT = ls_coldata[[i]]; colnames(COLDAT) = comparisons[i]

  idim(COLDAT)
  idim(DAT)

  DDSo =DESeqDataSetFromMatrix("countData" = DAT, "colData" = COLDAT, "design" = eval(parse(text = paste("~",comparisons[i]) ) ) )
  DE_objects_ls[[name]] = DESeq(DDSo, parallel = T, betaPrior = TRUE )

  name2=DE_results[i]
  res = iround(as.data.frame(results(DE_objects_ls[[name]], contrast=c(comparisons[i],"0", "1")) ))
  res =  res[order(res$log2FoldChange), ]
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
    plotMA(prepare4plotMA(rs, thr_padj_ = thr_padj, thr_log2fc_ = thr_log2fc), main=pname)
    ccc= m$s$cluster.colors$categ[comparisons.numeric[[i]]]
    wplot_Volcano.HRT(rs, thr_padj_ = thr_padj, thr_log2fc_ = thr_log2fc, pname_ = pname, Comp.Colors = ccc,saveit = F, Symmetric = T)
  }

  for (i in 1:length(comparisons) ) {
    name2=DE_results[i]
    rs = DE_results_ls[[name2]][ , COLZ.num ]
    pname = comparisons[i]
    ccc= m$s$cluster.colors$categ[comparisons.numeric[[i]]]
    wplot_Volcano.HRT(rs, thr_padj_ = thr_padj, thr_log2fc_ = 0, pname_ = pname, Comp.Colors = ccc,saveit = F, Symmetric = T)
    wplot_Volcano.HRT(rs, thr_padj_ = thr_padj, thr_log2fc_ = thr_log2fc, pname_ = pname, Comp.Colors = ccc,saveit = F, Symmetric = T, showNames = T)
  }

  pdfA4plot_off()
}


Fig.5A.and.S4AB.DE = T
if (Fig.5A.and.S4AB.DE) { try.dev.off()
  pname = "Fig.S5.D-E.Pairwise.Differerntial.Expression.Analysis";
  pdfA4plot_on(pname = pname, cols = 2, rows = 3)

  for (i in 1:length(comparisons) ) {
    name2=DE_results[i]
    rs = DE_results_ls[[name2]][ , COLZ.num ]
    pname = comparisons[i]
    plot.new()
    ccc= m$s$cluster.colors$categ[comparisons.numeric[[i]]]
    wplot_Volcano.HRT(rs, thr_padj_ = thr_padj, thr_log2fc_ = thr_log2fc, pname_ = pname, Comp.Colors = ccc,saveit = F, Symmetric = T)
  }
  pdfA4plot_off()
}





iprint("Plots are ready")
