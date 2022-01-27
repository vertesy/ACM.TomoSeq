
require(MarkdownReportsDev)

# GO_terms = read_clip_tbl()
# GO_termsFDR = as.named.vector(FirstCol2RowNames(GO_terms))
# View(GO_termsFDR)
# dput(GO_termsFDR)
GO_termsFDR = c(`vesicle-mediated transport` = 3.443697499, `intracellular transport` = 2.657577319,
                `immune system process` = 2.397940009, `chaperone-mediated autophagy` = 2.301029996,
                `G1/S transition of mitotic cell cycle` = 1.452225295, autophagy = 1.452225295,
                `central nervous system development` = 1.452225295, `antigen processing and presentation \nof exogenous peptide \nantigen via MHC class II` = 1.452225295,
                `mitotic cell cycle phase transition` = 1.452225295, `regulation of cell cycle` = 1.452225295,
                `regulation of protein localization \nto plasma membrane` = 1.452225295,
                `axon extension` = 1.330683119, `mitotic cell cycle` = 1.32330639,
                `protein-containing complex \nsubunit organization` = 1.311580178
)

try.dev.off()
wbarplot(GO_termsFDR, hline = -log10(c(0.05, 0.01)), ylab ="-log10(FDR)", w=5, lcol = c("red4", "red2")
         , filtercol = F, incrBottMarginBy = 8)
