#load in object
object <- readRDS(file = "output/diaz_astrooligo/diaz11612_oligo.rds")

# make dotplot per celltype with known marker genes
neuro <- DotPlot(object, features = c("RBFOX3", "MYT1", "CACNA1B", "WBSCR17", "FAM153B", "FRMPD4", "GRIN1",
                                                         "GABRG3", "MTUS2", "TENM2")) + labs(title = "neuron")
olig <- DotPlot(object, features = c("TMEM144", "DOCK5", "MBP", "FRMD4B", "ST18", "SHTN1", "PIP4K2A", "ST6GALNAC3",
                                     "MAN2A1", "OPALIN", "KANK4")) + labs(title = "oligo")
endo <- DotPlot(object, features = c("ABCB1", "CD34", "FLT4", "TIE1", "EGFL7", "ADGRL4", "ERG", "PTPRB", "CYYR1", "NOTCH4",
                                     "MECOM", "ANO2", "HSPG2", "BTNL9", "CDA")) + labs(title = "endo")
per <- DotPlot(object, features = c("CD248", "RGS5", "SLC38A11", "CARMN", "EDNRA", "PLXDC1", "GRM8", "PDGFRB",
                                     "TRPC6", "EBF2", "PRR16", "PRKG1", "C1S", "DCN")) + labs(title = "per")
tam <- DotPlot(object, features = c("CD163", "P2RY12", "CD14", "ITGB2", "C1QC", "SLC11A1")) + labs(title = "tam")
opc <- DotPlot(object, features = c("COL9A1", "PDGFRA", "SMOC1", "CA10", "PCDH15", "LINC02283", "AL512308.1",
                                     "ADARB2", "SCN9A")) + labs(title = "opc")
dividing <- DotPlot(object, features = c("TOP2A", "AURKB", "CDC20", "PLK1")) + labs(title = "dividing")
tumor <- DotPlot(object, features = c("EGFR", "GFAP", "SOX2", "VIM", "PTPRZ1", "TENM1")) + labs(title = "tumor")
astro <- DotPlot(object, features = c("LRRC3B", "GABRG1", "ETNPPL", "LINC00499", "TPD52L1", "AC002429.2", "SLC39A12",
                                     "EMX2", "SPON1", "WIF1", "LINC00943", "GJB6", "ZNF98", "FOXG1", "SLC7A10")) + labs(title = "astro")
tcell <- DotPlot(object, features = c("CD2", "CD3D", "TRBC2", "TRAC", "ICOS", "GZMA", "SKAP1", "CD96",
                                      "THEMIS", "SLFN12L")) + labs(title = "tcell")
bcell <- DotPlot(object, features = c("IGLC3", "CD19", "CD79B")) + labs(title = "bcell")
fibro <- DotPlot(object, features = c("CYP1B1", "OGN", "CLU", "CXCL12")) + labs(title = "fibro")

# # plot 2 or more dotplots next to each other -- used for thesis
# ggpubr::ggarrange(tumor, dividing, ncol = 2, nrow = 1, align = 'h')
# ggpubr::ggarrange(astro, neuro, ncol = 2, nrow = 1, align = 'h')
# ggpubr::ggarrange(per, tam, ncol = 2, nrow = 1, align = 'h')
# ggpubr::ggarrange(olig, endo, ncol = 2, nrow = 1, align = 'h')
  
# show the dotplots -- to analyze
print(neuro)
print(olig)
print(endo)
print(per)
print(tam)
print(opc)
print(dividing)
print(tumor)
print(astro)
print(tcell)
print(fibro)


# # tried to subclassify tumor clusters
# mesen <- DotPlot(object, idents = c(0, 4, 5, 6, 7, 8, 10), features = c("YKL40", "MET", "CD44", "RELB", "F13A1", "RNF149", "PLAUR", "CASP4", "ILR4", "TRADD")) + labs(title = "mesen")
# print(mesen)
# pron <- DotPlot(object, idents = c(0, 4, 5, 6, 7, 8, 10), features = c("PDGFRA", "OLIG2", "DDL3", "SOX2", "NKX2", "UBE2E2", "NKAIN", "DLL3", "EBRB3")) + labs(title = "proneural")
# print(pron)
# # neur <- DotPlot(object, idents = c(0, 1, 2, 3, 8, 10, 12), features = c("NEFL", "SLC12A5", "SYT1", "GABRA1")) + labs(title = "neural")
# # print(neur)
# clas <- DotPlot(object, idents = c(0, 4, 5, 6, 7, 8, 10), features = c("EGFR", "AKT2", "SMO", "NOTCH3", "JAG1", "FGFR3", "PDGFA", "NES")) + labs(title = "classical")
# print(clas)
# # 
# 
