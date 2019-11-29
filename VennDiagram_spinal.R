

s_115 = read.csv("Mus_s_115_ls and spinalcord/to ENTREZ_118to115.csv")
s_115_ls = s_115[2]

s_299 = read.csv("Mus_s_115_ls and spinalcord/to ENTREZ_304to299.csv")
s_299_ls = s_299[2]

library(VennDiagram)

s_115_list = as.vector(unlist(s_115_ls))

s_299_list = as.vector(unlist(s_299_ls))


venn.diagram(list(spinalcord_118=s_115_list,spinalcord_D126_304=s_299_list),
             
             resolution = 500,imagetype = "tiff",alpha=c(0.5,0.5),
             fill=c("yellow","green"),
             
             main="Mus_spinalcord difer tissue DEG",
             
             
             filename = "VennDiagram_spinal.tif")
