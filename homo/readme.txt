#代码1是limma筛选DEG前 未去除样本分布差异（批次效应），整个流程
       2为查看样本数据分布和分组情况，及处理后查看
       3为消除样本分布差异normalizeBetweenArrays()处理后的筛选

1、筛选差异表达基因
#芯片数据：control：正常神经元的表达谱（对照）
                SOD1：sod1相关肌萎缩性脊髓侧索硬化症(ALS)患者分离的运动神经元的表达谱

#代码编码为GB2312，用UTF-8中文注释乱码。
第三方包Bioconductor的下载： http://www.bioconductor.org/install/

问题：原始下载的探针基因表达矩阵中，表达值从0.1-几万。

2、用筛选得差异表达基因富集pathways，提取pathway上所有基因的表达值。




