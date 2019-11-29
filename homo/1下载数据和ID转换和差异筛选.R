#####################################获得芯片表达矩阵,a1#########################################
library(GEOquery)
gset = getGEO('GSE20589',destdir = '.',AnnotGPL = 'F',getGPL = 'F')#将GSE20589数据下载到R里当前工作目录并赋值给gset
#class(gset) 
#gseet对象里包含着各种各样的信息：表达矩阵、芯片如何设计的、样本如何分组etc.gset是一个大列表，我们需要从中提取出表达矩阵
#数据使用Microarray Suite 5.0 (MAS5)进行分析，使用Affymetrix默认分析设置和全局缩放作为标准化方法.每个阵列的平均目标强度被任意设置为100

b = gset[[1]]
a1 = exprs(b)   #提取eSet列表里的第一个元素：eSet[[1]]；并使用exprs函数把它转化成矩阵

#####################################过滤探针##########################################
#安装包'org.Hs.eg.db'和'hgu133plus2.db'
library(hgu133plus2.db)
ids <- toTable(hgu133plus2SYMBOL)   #得到probe_id和symbol的对应关系要用hgu133plus2SYMBOL数据集，用toTable提取数据集里面的信息
plot(table(sort(table(ids$symbol))))#查看基因对应probe_id数量

###过滤无对应基因探针
#显示分布
table(rownames(a1)%in% ids$probe_id) #12734个探针不对应基因
table(sort(table(ids$symbol)))
tail(sort(table(ids$symbol)))
#过滤
dim(a1)            #过滤前54675行
a1 = a1[rownames(a1)%in% ids$probe_id,]
dim(a1)            #过滤后41922个探针

#改变探针顺序和表达矩阵a1中探针顺序一致
ids <- ids[match(rownames(a1),ids$probe_id),]
head(ids)
a1[1:5,]

#修改a1，保留对应多个探针的symbol中平均表达最高的探针
tmp <- by(a1,
  ids$symbol,
  function(x) rownames(x)[which.max(rowMeans(x))]
  )
#将同一个symbol所对应的多个探针分成不同的组，并对每组探针进行统计：
#1、计算每组中每行探针表达量的平均值（也就是每个探针在10个样本中表达量的均值rowMeans(x)），
#2、再取平均值最大的那个探针作为该symbol所对应的唯一探针，该组中的其它探针过滤掉

probes<- as.character(tmp)
a1 <- a1[rownames(a1) %in% probes ,]
dim(a1)#20174个探针


########################################转换探针ID后的表达矩阵new_a1##################################################
new_ids <- ids[match(rownames(a1),ids$probe_id),]
new_a1 <- a1[new_ids$probe_id,]           #将a1按照 从ids取出probe_id这一列中的每一行组成一个新的new_a1
rownames(new_a1) <-new_ids$symbol         #把ids的symbol这一列中的每一行给new_a1作为new_a1的行名
dim(new_a1)



########################################从b提取样本分组信息################################################
#tmp2 <-pData(b)
#grouplist <- tmp2[,1]
#group_list <- as.character(tmp[,1])
group_list <-c(rep('SOD1',3),rep('control',7))
group_list


#limma分析前对数据log2处理(下载矩阵的标准化处理方法MAS5返回值未经过log2转换)
log2_a1 <- log2(new_a1)
########################################画图看各个样本的表达量-看表达矩阵的分布图##########################
#使用ggplot2画各个样本表达量的boxplot图
# 准备画图所需数据log2_a1_L
library(reshape2)
head(log2_a1)
log2_a1_L = melt(log2_a1)
head(log2_a1_L)
colnames(log2_a1_L)=c('symbol','sample','value')
head(log2_a1_L)

# 获得分组信息
library(stringr)
#group_list = ifelse(str_detect(tmp2$title,"Control")==TRUE,"control","case")
#group_list
log2_a1_L$group = rep(group_list,each=nrow(log2_a1))
head(log2_a1_L)

# ggplot2画图 
library(ggplot2)
p = ggplot(log2_a1_L,
           aes(x=sample,y=value,fill=group))+geom_boxplot()
print(p)

#######################################检查样本分组信息,一般看PCA图，hclust图###################################

#histogram
p=ggplot(log2_a1_L,aes(value,fill=group))+geom_histogram(bins = 200)+facet_wrap(~sample, nrow = 4)
print(p)
#发现GSM480309分布与其他对照组有区别，考虑是否删除,

#hclust
hc<- hclust(dist(t(log2_a1)))
plot(hc)#结果不好

#PCA,没有掌握
pc <- prcomp(t(new_a1),scale=TRUE)
pcx=data.frame(pc$x)
pcr=cbind(samples=rownames(pcx),grouplist, pcx) 
p=ggplot(pcr, aes(PC1, PC2))+geom_point(size=5, aes(color=group_list)) +
  geom_text(aes(label=samples),hjust=-0.1, vjust=-0.3)
print(p)#结果130单独分出？


########################################用limma对基因差异表达做分析，及可视化####################################
#用limma做差异分析##三个矩阵：表达矩阵(log2_a1)、分组矩阵(design)、差异比较矩阵（contrast.matrix）
                  ##三个步骤：lmFit(这类对象通常是每个唯一的探针包含一行),eBayes,topTable

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
log2_a1 <- log2(new_a1)
#step1
fit1 <- lmFit(log2_a1,design)
#step2
fit2 <- contrasts.fit(fit1,contrast.matrix)
fit2 <- eBayes(fit2)
#step3
tempOutput <- topTable(fit2,coef = 1,n = Inf)
nrDEG <- na.omit(tempOutput)
head(nrDEG)

#heatmap 

#以2和0.05分别作为logFC和P.value的阈值筛选nrDEG得554个gene
new_nrDEG <- nrDEG[which((nrDEG$logFC>=2|nrDEG$logFC<=-2)&nrDEG$P.Value<=0.05),]
#以2和0.01分别作为logFC和P.value的阈值筛选nrDEG得198个gene
new_nrDEG2 <- nrDEG[which((nrDEG$logFC>=2|nrDEG$logFC<=-2)&nrDEG$P.Value<=0.01),]
#以2.5和0.01分别作为logFC和P.value的阈值筛选nrDEG得119个gene
new_nrDEG3 <- nrDEG[which((nrDEG$logFC>=2.5|nrDEG$logFC<=-2.5)&nrDEG$P.Value<=0.01),]
#以3和0.01分别作为logFC和P.value的阈值筛选nrDEG得50个gene
new_nrDEG4 <- nrDEG[which((nrDEG$logFC>=3|nrDEG$logFC<=-3)&nrDEG$P.Value<=0.01),]
library(pheatmap)
choose_gene=head(rownames(new_nrDEG),554) #取多少基因做热图？
choose_matrix=new_a1[choose_gene,]
choose_matrix=t(scale(t(choose_matrix)))
pheatmap(choose_matrix)
#pheatmap(choose_matrix,display_numbers = TRUE,number_color = "blue")

#火山图

plot(nrDEG$logFC,nrDEG$P.Value)
plot(nrDEG$logFC,-log10(nrDEG$P.Value))






###############################################################GSEA######################################################################
#KEGG pathway analysis
gene <-head(rownames(nrDEG),1000)
library(clusterProfiler)
#使用clsuterProfiler包需要ENTREZID
gene.df <- bitr(gene,fromType = "SYMBOL",
                toType = c("ENSEMBL","ENTREZID"),
                OrgDb = org.Hs.eg.db)
head(gene.df)#6.96% of input gene IDs are fail to map...

kk<-enrichKEGG(gene = gene.df$ENTREZID,
               organism = 'hsa',
               pvalueCutoff = 0.05)
head(kk@result)[,1:6]

genelist = nrDEG$logFC
names(genelist) = rownames(nrDEG)
genelist = sort(genelist,decreasing = T)

kk2 <-gseKEGG(geneList = genelist,
              organism = 'hsa',
              nPerm = 1000,
              minGSSize = 120,
              pvalueCutoff = 0.05,
              )
