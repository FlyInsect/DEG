#用ggplot2画各个样本表达量的boxplot
#从图中看到两个分组control和treat基本在一条线上，这样的数据才可以进行后续比较，
#如果不在一条线上说明有批次效应batch infect，需要用limma包内置函数normalizeBetweenArrays人工校正一下(Normalization)：

library(limma) 
log2_a1_nor = normalizeBetweenArrays(log2_a1)
#normalizeBetweenArrays(object, method=NULL, targets=NULL, cyclic.method="fast", ...):参数method=" Scale "将列缩放到具有相同的中位数.
#如果object是矩阵：则假定它包含对数转换的单通道数据
#method：指定要使用的规范化方法。对于单通道对象，默认是“分位数”
#如果对象是EListRaw对象(原始表达列表.一种用于存储标准化之前单通道原始强度的类,强度未取对数.此类对象包含一行代表每个探针,一列代表每个微阵列)
                       #那么将对表达式值的矩阵对象$E应用标准化，然后对其进行log2转换.
#value:如果object是一个矩阵，那么将生成一个大小相同的矩阵.如果对象是EListRaw对象，那么将生成一个EList对象，其表达式值为log2.

boxplot(log2_a1_nor)


#使用ggplot2画各个样本表达量的boxplot图
# 准备画图所需数据log2_a1_nor_L
library(reshape2)
head(log2_a1_nor)
log2_a1_nor_L = melt(log2_a1_nor)
head(log2_a1_nor_L)
colnames(log2_a1_nor_L)=c('symbol','sample','value')
head(log2_a1_nor_L)

# 获得分组信息
library(stringr)
group_list2 = ifelse(str_detect(tmp2$title,"Control")==TRUE,"control","treat")
group_list2
log2_a1_nor_L$group = rep(group_list,each=nrow(log2_a1_nor))
head(log2_a1_nor_L)

# ggplot2画图 
library(ggplot2)
p = ggplot(log2_a1_nor_L,
           aes(x=sample,y=value,fill=group))+geom_boxplot()
print(p)


##############################################绘制PCA,hclust看样本分组###############################################
pc <- prcomp(t(log2_a1_nor),scale=TRUE)
pcx=data.frame(pc$x)
pcr=cbind(samples=rownames(pcx),grouplist, pcx) 
q=ggplot(pcr, aes(PC1, PC2))+geom_point(size=5, aes(color=group_list)) +
  geom_text(aes(label=samples),hjust=-0.1, vjust=-0.3)
print(q)

library(ggfortify)
# 互换行和列，再dim一下
df=as.data.frame(t(log2_a1_nor))
# 不要view df，列太多，软件会卡住；
dim(df)
dim(log2_a1_nor)

log2_a1_nor[1:7,1:7]
df[1:7,1:7]

df$group=group_list 
autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')#校正前后结果都不好



######################################hclust
# 更改表达矩阵列名
head(log2_a1_nor)
colnames(log2_a1_nor) = paste(group_list,1:7,sep='')
head(log2_a1_nor)
# 定义nodePar
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                cex = 0.7, col = "blue")
# 聚类
hc=hclust(dist(t(log2_a1_nor)))
par(mar=c(5,5,5,10)) 
# 绘图
plot(as.dendrogram(hc), nodePar = nodePar,  horiz = TRUE) #校正后区分不开case和control
