#####################################获得芯片表达矩阵,a1#########################################
library(GEOquery)
gset = getGEO('GSE3343',destdir = '.',AnnotGPL = 'F',getGPL = 'F')#将GSE20589数据下载到R里当前工作目录并赋值给gset
#class(gset) 
#gset对象里包含着各种各样的信息：表达矩阵、芯片如何设计的、样本如何分组etc.gset是一个大列表，我们需要从中提取出表达矩阵
#数据使用Microarray Suite 5.0 (MAS5)进行分析，使用Affymetrix默认分析设置和全局缩放作为标准化方法.每个阵列的平均目标强度被任意设置为100

b = gset[[1]]
a1 = exprs(b)   #提取eSet列表里的第一个元素：eSet[[1]]；并使用exprs函数把它转化成矩阵
a1 = a1[,c(16:18,11:12)]
head(a1)

#####################################过滤探针##########################################
#安装包'org.Mm.eg.db'和'moe430a.db'
library(org.Mm.eg.db)
library(moe430a.db)
ids <- toTable(moe430aSYMBOL)   #得到probe_id和symbol的对应关系要用mouse430a2SYMBOL数据集，用toTable提取数据集里面的信息
plot(table(sort(table(ids$symbol))))#查看基因对应probe_id数量

###过滤无对应基因探针
#显示分布
table(rownames(a1)%in% ids$probe_id) #1460个探针不对应基因
table(sort(table(ids$symbol)))
tail(sort(table(ids$symbol)))
#过滤
dim(a1)            #过滤前22690行
a1 = a1[rownames(a1)%in% ids$probe_id,]
dim(a1)            #过滤后21230个探针

#改变探针顺序和表达矩阵a1中探针顺序一致
ids <- ids[match(rownames(a1),ids$probe_id),]
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
dim(a1)#12993个探针


########################################转换探针ID后的表达矩阵new_a1##################################################
new_ids <- ids[match(rownames(a1),ids$probe_id),]
new_a1 <- a1[new_ids$probe_id,]           #将a1按照 从ids取出probe_id这一列中的每一行组成一个新的new_a1
rownames(new_a1) <-new_ids$symbol         #把ids的symbol这一列中的每一行给new_a1作为new_a1的行名
dim(new_a1)



########################################从b提取样本分组信息################################################
tmp2 <-pData(b)
grouplist <- tmp2[,1]
#group_list <- as.character(tmp[,1])
group_list <-c(rep('sod1_mutant',3),rep('control',2))
group_list


#limma分析前可能需对数据log2处理(下载矩阵的标准化处理方法MAS5返回值未经过log2转换)
new_a1 <- log2(new_a1)
########################################画图看各个样本的表达量-看表达矩阵的分布图##########################
#使用ggplot2画各个样本表达量的boxplot图
# 准备画图所需数据new_a1_L
library(reshape2)
head(new_a1)
new_a1_L = melt(new_a1)
head(new_a1_L)
colnames(new_a1_L)=c('symbol','sample','value')
head(new_a1_L)

# 获得分组信息
#library(stringr)
#group_list = ifelse(str_detect(tmp2$title,"control")==TRUE,"control","sod1_mutant")
#group_list
new_a1_L$group = rep(group_list,each=nrow(new_a1))
head(new_a1_L)

# ggplot2画图 
library(ggplot2)
p = ggplot(new_a1_L,
           aes(x=sample,y=value,fill=group))+geom_boxplot()
print(p)

#######################################检查样本分组信息,一般看PCA图，hclust图###################################

#histogram
p=ggplot(new_a1_L,aes(value,fill=group))+geom_histogram(bins = 200)+facet_wrap(~sample, nrow = 4)
print(p)

#hclust
hc<- hclust(dist(t(new_a1)))
plot(hc)#结果不好

#PCA,没有掌握
pc <- prcomp(t(new_a1),scale=TRUE)
pcx=data.frame(pc$x)
pcr=cbind(samples=rownames(pcx),group_list, pcx) 
p=ggplot(pcr, aes(PC1, PC2))+geom_point(size=5, aes(color=group_list)) +
  geom_text(aes(label=samples),hjust=-0.1, vjust=-0.3)
print(p)



#########################################################校正批次效应并查看分布和分组情况#############################################
#用ggplot2画各个样本表达量的boxplot
#从图中看到两个分组control和treat基本在一条线上，这样的数据才可以进行后续比较，
#如果不在一条线上说明有批次效应batch infect，需要用limma包内置函数normalizeBetweenArrays人工校正一下(Normalization)：

library(limma) 
new_a1_nor = normalizeBetweenArrays(new_a1)
#normalizeBetweenArrays(object, method=NULL, targets=NULL, cyclic.method="fast", ...):参数method=" Scale "将列缩放到具有相同的中位数.
#如果object是矩阵：则假定它包含对数转换的单通道数据
#method：指定要使用的规范化方法。对于单通道对象，默认是“分位数”
#如果对象是EListRaw对象(原始表达列表.一种用于存储标准化之前单通道原始强度的类,强度未取对数.此类对象包含一行代表每个探针,一列代表每个微阵列)
#那么将对表达式值的矩阵对象$E应用标准化，然后对其进行log2转换.
#value:如果object是一个矩阵，那么将生成一个大小相同的矩阵.如果对象是EListRaw对象，那么将生成一个EList对象，其表达式值为log2.

boxplot(new_a1_nor)


# 准备画图所需数据类型new_a1_nor_L
library(reshape2)
head(new_a1_nor)
new_a1_nor_L = melt(new_a1_nor)
head(new_a1_nor_L)
colnames(new_a1_nor_L)=c('symbol','sample','value')
head(new_a1_nor_L)

# 获得分组信息
library(stringr)
#group_list2 = ifelse(str_detect(tmp2$title,"control")==TRUE,"control","treat")
#group_list2
new_a1_nor_L$group = rep(group_list,each=nrow(new_a1_nor))
head(new_a1_nor_L)

# ggplot2画图 
library(ggplot2)
p = ggplot(new_a1_nor_L,
           aes(x=sample,y=value,fill=group))+geom_boxplot()
print(p)


##############################################绘制PCA,hclust看样本分组###############################################
pc <- prcomp(t(new_a1_nor),scale=TRUE)
pcx=data.frame(pc$x)
pcr=cbind(samples=rownames(pcx),group_list, pcx) 
q=ggplot(pcr, aes(PC1, PC2))+geom_point(size=5, aes(color=group_list)) +
  geom_text(aes(label=samples),hjust=-0.1, vjust=-0.3)
print(q)

library(ggfortify)
# 互换行和列，再dim一下
df=as.data.frame(t(new_a1_nor))
# 不要view df，列太多，软件会卡住；
dim(df)
dim(new_a1_nor)

new_a1_nor[1:3,1:2]
df[1:3,1:2]

df$group=group_list 
autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')#校正前后结果都不好



######################################hclust
# 更改表达矩阵列名
head(new_a1_nor)
colnames(new_a1_nor) = paste(group_list,1:5,sep='')
head(new_a1_nor)
# 定义nodePar
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                cex = 0.7, col = "blue")
# 聚类
hc=hclust(dist(t(new_a1_nor)))
par(mar=c(5,5,5,10)) 
# 绘图
plot(as.dendrogram(hc), nodePar = nodePar,  horiz = TRUE) #校正后区分不开case和control





####################################################校正批次效应后差异基因筛选#################################################
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
#new_a1 <- log2(new_a1)
#step1
fit1_nor <- lmFit(new_a1_nor,design)
#step2
fit2_nor <- contrasts.fit(fit1_nor,contrast.matrix)
fit2_nor <- eBayes(fit2_nor)
#step3
tempOutput_nor <- topTable(fit2_nor,coef = 1,n = Inf)
nrDEG_nor <- na.omit(tempOutput_nor)
head(nrDEG_nor)

#heatmap 

#以1.5和0.05分别作为logFC和P.value的阈值筛选校正过的nrDEG_nor得109个gene
new_nrDEG_nor <- nrDEG_nor[which((nrDEG_nor$logFC>=1.5|nrDEG_nor$logFC<=-1.5)&nrDEG_nor$P.Value<=0.05),]
dim(new_nrDEG_nor)

library(pheatmap)
choose_gene=head(rownames(new_nrDEG_nor),118) #取多少基因做热图？
choose_matrix=new_a1[choose_gene,]
choose_matrix=t(scale(t(choose_matrix)))
pheatmap(choose_matrix)

#write.csv(nrDEG_nor,"F:/RstudioWorkspace/sod1 and ALS/result/nrDEG_nor.csv")


