library(limma)

#分组矩阵
design <- model.matrix(~0+factor(group_list))
colnames(design) <- levels(factor(group_list))
rownames(design) <-colnames(new_a1)
design

#比较矩阵
contrast.matrix <-makeContrasts(paste0(unique(group_list),collapse = "-")
                                ,levels = design)
contrast.matrix



#limma分析前对数据log2处理(下载矩阵的标准化处理方法MAS5返回值未经过log2转换)
#log2_a1 <- log2(new_a1)
#step1
fit1_nor <- lmFit(log2_a1_nor,design)
#step2
fit2_nor <- contrasts.fit(fit1_nor,contrast.matrix)
fit2_nor <- eBayes(fit2_nor)
#step3
tempOutput_nor <- topTable(fit2_nor,coef = 1,n = Inf)
nrDEG_nor <- na.omit(tempOutput_nor)
head(nrDEG_nor)

#heatmap 

#以2和0.05分别作为logFC和P.value的阈值筛选校正过的nrDEG_nor得613个gene
new_nrDEG_nor <- nrDEG_nor[which((nrDEG_nor$logFC>=2|nrDEG_nor$logFC<=-2)&nrDEG_nor$P.Value<=0.05),]
dim(new_nrDEG_nor)


library(pheatmap)
choose_gene=head(rownames(new_nrDEG_nor),613) #取多少基因做热图？
choose_matrix=new_a1[choose_gene,]
choose_matrix=t(scale(t(choose_matrix)))
pheatmap(choose_matrix)

#write.csv(nrDEG_nor,"F:/RstudioWorkspace/sod1 and ALS/result/nrDEG_nor.csv")
