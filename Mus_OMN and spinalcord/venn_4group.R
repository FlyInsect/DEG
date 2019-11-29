library(VennDiagram)
OMN = read.csv("OMN_106.csv")

OMN_list = as.vector(unlist(OMN["X"]))

spinalcord_118 = read.csv("spinalcord_118.csv")
spinalcord_118_list = as.vector(unlist(spinalcord_118["X"]))

MN = read.csv("../Mus_MN/nrDEG_109.csv")
MN_list =as.vector(unlist(MN[1]))

spinalcord_D126_304 = read.csv("../Mus_spinalcord/D126_304.csv")
spinalcord_D126_304_list =as.vector(unlist(spinalcord_D126_304[1]))

venn.diagram(list(OMN=OMN_list,spinalcord_118=spinalcord_118_list,MN=MN_list,spinalcord_D126_304=spinalcord_D126_304_list),
             
             resolution = 500,imagetype = "tiff",alpha=c(0.5,0.5,0.5,0.5),
             fill=c("red","yellow","blue","green"),
             
             main="Mus difer tissue DEG",
          
             
             filename = "VennDiagram2.tif")

#write.table(spinalcord_D126_304_list,file="spinalcord_D126_304.txt" , sep =" ", row.names =FALSE,col.names =FALSE, quote =FALSE)
