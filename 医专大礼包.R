library(GEOquery)
library(tidyverse)
library(dplyr)
library(limma)
library(glmnet)
library(rjson)
library(XML)
library(methods)
library(pheatmap)
library(stringr)
library(ggplot2)
library(survival)
library(survminer)
library(stringi)	
library(Hmisc)
library(mgcv)
library(nlme)
library(sva)
library(tinyarray)
library(rms)
library(RColorBrewer)
library(gridExtra)
library(pROC)
library(clusterProfiler)
library(enrichplot)
library(GOplot)
library(timeROC)
library(regplot)
library(caret)
library(lattice)
library(org.Hs.eg.db)
library(readxl)
library(sm)
library(zoo)
library(usethis)
library(devtools)
library(e1071)
library(preprocessCore)
library(parallel)
library(CIBERSORT)
library(corrplot)
library(vioplot)
library(reshape2)
library(ggpubr)
library(ggcorrplot)
library(seqinr)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(mixOmics)
library(survcomp)
library(CoxBoost)
library(survivalsvm)
library(BART)
library(snowfall)
library(ComplexHeatmap)
library(autoReg)
library(dcurves)
library(cowplot)
library(GSVA)
library(GSEABase)
library(ggDCA)
library(patchwork)
library(ggsci)
library(ggVolcano)
library(Rtsne)
library(ConsensusClusterPlus)
library(ggrepel)


##################
#1.差异分析

#读取输入文件
data=read.table("1.TCGA原始表达矩阵.txt", header=T, sep="\t", check.names=F,row.names = 1)
dimnames=list(rownames(data), colnames(data))
data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)
data=data[rowMeans(data)>1,]
#正常和肿瘤数目,第14,15字符,01-09是癌症，10-19是正常，20-29是癌旁
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
#将2变为1
group=gsub("2", "1", group)
#获取正常及肿瘤组样本数目
conNum=length(group[group==1])
treatNum=length(group[group==0])
Type=c(rep(1,conNum), rep(2,treatNum))
#根据正常和肿瘤排序
data1 = data[,group == 1]
data2 = data[,group == 0]
data = cbind(data1,data2)

outTab=data.frame()
for(i in row.names(data)){
  rt=data.frame(expression=data[i,], Type=Type)
  wilcoxTest=wilcox.test(expression ~ Type, data=rt)
  pvalue=wilcoxTest$p.value
  conGeneMeans=mean(data[i,1:conNum])
  treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)
  conMed=median(data[i,1:conNum])
  treatMed=median(data[i,(conNum+1):ncol(data)])
  diffMed=treatMed-conMed
  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){
    outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))}
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)), method="fdr")
outTab=cbind(outTab, fdr=fdr)

# 输出差异基因
logFCfilter=1
fdrFilter=0.05
outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(outDiff,file="(2)TCGA.diff.Wilcoxon.txt",sep="\t",row.names=F,quote=F)
write.table(outTab,file="(3)TCGA.all.diff.Wilcoxon.txt",sep="\t",row.names=F,quote=F)

#  热图
#设置展示基因的数目

geneNum=30
outDiff=outDiff[order(as.numeric(as.vector(outDiff$logFC))),]
diffGeneName=as.vector(outDiff[,1])
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(2*geneNum)){
  hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
  hmGene=diffGeneName
}
hmExp=log2(data[hmGene,]+0.01)
Type=c(rep("Normal",conNum),rep("Tumor",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c(rep("#0477A8",5), "white", rep("#9E024E",5)))(100),
         cluster_cols =F,
         cluster_rows =T,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=5,
         fontsize_col=8)


#  火山图

###-----------------------------------------------------------------------------
### 设置参数
###-----------------------------------------------------------------------------
deg_res <- read.table("(3)TCGA.all.diff.Wilcoxon.txt", header=T, sep="\t", check.names=F,row.names = 1)
colnames(deg_res)

###-----------------------------------------------------------------------------
### 设置参数
###-----------------------------------------------------------------------------
x = "logFC"
y = "fdr"
# y = "adj.P.Val"
abs_logFC_cut = 1
p_cut = 0.05

mark_genes <- NULL      # NULL 或者输入要标记的基因组成的向量
sort_by = "fdr"     # 按照什么排序筛选要标记的topgene
split = T               # 上下调是否分开排序
mark_top = 10                # top几
point_size = 1
outpath = "./"


###-----------------------------------------------------------------------------
### 绘制火山图的代码 不用修改！
###-----------------------------------------------------------------------------

value2 <- 1:10
rr <- range(value2)
svals <- (value2-rr[1])/diff(rr)
f <- colorRamp(c("#0477A8","#DFCF0B","#9E024E"))     ## 渐变颜色设置
cols <- rgb(f(svals)/255) 

change_colname = "Change"
deg_res[,change_colname] <- ifelse(deg_res[,y] > p_cut,'Not',                                        
                                   ifelse( deg_res[,x] >= abs_logFC_cut,'Up',                               
                                           ifelse( deg_res[,x] <= -abs_logFC_cut,'Down','Not') )  )    

top = mark_top
color_colname = y
p <- ggplot(data = deg_res, aes(x = .data[[x]], y = -log10(.data[[y]]))) +
  geom_point(aes(color = -log10(.data[[color_colname]]),
                 size = -log10(.data[[color_colname]])
  ) ,alpha = 1
  )+
  
  geom_vline(xintercept = abs_logFC_cut,color = "grey10", linetype="longdash",size = 0.4) +
  geom_vline(xintercept = -abs_logFC_cut,color = "grey10", linetype="longdash",size = 0.4) +
  
  geom_hline(yintercept = -log10(p_cut) ,color = "grey10",linetype="longdash", size = 0.4) +
  
  theme_bw()+
  theme(
    #legend.position = "top",
    #legend.justification = c(0, 1)
  )+
  #下面注释的这行代码会改变color的图例
  # guides(color = guide_legend(override.aes = list(size = 3),title = ""))+
  scale_color_gradientn(values = seq(0,1,0.2),colors = cols)+
  scale_size(range = c(0.5,3))

p
#ggsave(filename = paste0(outpath,"/deg_volcano_plot_clean.pdf"),width = 7,height = 6)


if (length(mark_genes) != 0) {
  mark_genes <-  mark_genes
}else{
  
  if (split == F) {
    # 方法二 不区分上下调 筛选top gene
    if (all(deg_res[,sort_by] >= 0)) {
      # p值的top gene
      mark_genes <- 
        deg_res %>% filter(Change %in% c("Up","Down")) %>% 
        arrange(.data[[sort_by]]) %>% head(top) %>% rownames()
      
    }else{
      # logfc的top gene
      mark_genes <- 
        deg_res %>% filter(Change %in% c("Up","Down")) %>% mutate(tmp = abs(.data[[sort_by]])) %>% arrange(desc(tmp)) %>% head(top) %>% rownames()
    }
    
  }else{
    
    # 方法三 分别筛选出上下调的top的基因
    if (all(deg_res[,sort_by] >= 0)){
      # p值的top gene
      mark_genes_up <- 
        deg_res %>% filter(Change == "Up") %>% arrange(.data[[sort_by]]) %>% head(top) %>% rownames()
      
      mark_genes_down <- 
        deg_res %>% filter(Change == "Down") %>% arrange(.data[[sort_by]]) %>% head(top) %>% rownames()
      
      mark_genes <- c(mark_genes_up,mark_genes_down)
      
    }else{
      # logfc
      mark_genes_up <- 
        deg_res %>% filter(Change == "Up") %>% arrange(desc(.data[[sort_by]])) %>% head(top) %>% rownames()
      
      mark_genes_down <- 
        deg_res %>% filter(Change == "Down") %>% arrange(.data[[sort_by]]) %>% head(top) %>% rownames()
      
      mark_genes <- c(mark_genes_up,mark_genes_down)
      
    }
  }
}

# 标记基因的差异结果
show_df <- deg_res[mark_genes,]
show_df$gene_name <- rownames(show_df)

Up_label_data <- show_df[show_df[,x] > 0,] 
Down_label_data <- show_df[show_df[,x] < 0,]

distance = 0.05
# 设置mark gene 的label到logfc最大的点的距离
distance = distance                      # 点和label的距离

max_logFC <- deg_res[,x] %>% max()
min_logFC <- deg_res[,x] %>% min()

distance
p1 <- p +
  #scale_x_continuous(expand = c(0.3,0.3)) +
  #scale_y_continuous(expand = expansion(add = c(0.5,1))) +
  ggrepel::geom_label_repel(
    data = Up_label_data,
    size = 3,
    aes(label = gene_name,color = -log10(.data[[y]])),
    nudge_y      = 0 ,
    direction    = "y",
    hjust        = 0,            # 0 左对齐
    # segment.color = point_color[3],
    segment.size = 0.1,
    segment.linetype = 1,
    max.overlaps = 10,
    max.iter = 1000000,
    max.time = 10,
    nudge_x =  max_logFC + distance - Up_label_data[,x],
    min.segment.length = 0,
    # fontface = "bold"
  )+
  ggrepel::geom_label_repel(
    
    data = Down_label_data,
    size = 3,
    aes(label = gene_name,color = -log10(.data[[y]])),
    nudge_y      = 0 ,
    direction    = "y",
    hjust        = 1,    # 1 右对齐
    segment.size = 0.1,
    # segment.color = point_color[1],
    segment.linetype = 1,
    max.overlaps = 10,
    max.iter = 1000000,
    max.time = 10,
    nudge_x =  min_logFC - distance - Down_label_data[,x],
    min.segment.length = 0,
    # fontface = "bold"
  ) 

p1
#ggsave(filename = paste0(outpath,"/deg_volcano_plot_add_label.pdf"),width = 8,height = 8)

#2.单因素cox 

#表达矩阵数据
data=read.table("4.三重筛选预后前表达矩阵(log2转换).txt", header=T, sep="\t", check.names=F,row.names = 1)
dimnames=list(rownames(data), colnames(data))
data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
data = data[,group == 0]
data=t(data)
rownames(data)=substr(rownames(data),1,12)

#生存数据
cli=read.table("TCGA临床数据 - (短).txt", header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
rt=cbind(cli, data)

write.table(rt,"5.三重筛选表达矩阵和临床数据合并.txt", sep="\t", quote=F, col.names=T,row.names = T)


p.value = 0.01
outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(time, state) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP < p.value){
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}

rt = cbind(rt[,c(1,2)],rt[,outTab$id])
write.table(outTab, file="(5)单因素cox基因.txt", sep="\t", row.names=F, quote=F)
write.table(rt, file="6.单因素cox基因表达矩阵和临床数据.txt", sep="\t", row.names=T, quote=F)



############定义森林图函数############
bioForest=function(coxFile=null, forestFile=null, forestCol=null){
  #读取输入文件
  rt <- read.table(coxFile, header=T, sep="\t", check.names=F, row.names=1)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  #输出图形
  height=nrow(rt)/12.5+5
  pdf(file=forestFile, width=7, height=height)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #绘制森林图左边的基因信息
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
  
  #对HR进行可视化,绘制右边的森林图
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol[1], forestCol[2])
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.6)
  axis(1)
  dev.off()
}

#调用函数, 对预后相关的基因进行可视化, 得到森林图
bioForest(coxFile="(5)单因素cox基因.txt", forestFile="[3]单因素cox森林图.pdf", forestCol=c("#FF9289","#00DAE0"))



#3.lasso

rt=read.table("7.训练集表达矩阵和临床数据.txt", header=T, sep="\t", check.names=F,row.names = 1)
x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$time, rt$state))
fit=glmnet(x, y, family="cox", nfolds = 10)

#c-index
cvfit <- cv.glmnet(x, y, family = "cox", type.measure = "C", nfolds = 10)
pdf("lasso.c-index.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)), lty="dashed")
dev.off()

#deviance图形
cvfit=cv.glmnet(x, y, family="cox", type.measure = "deviance", nfolds = 10)
pdf("lasso.cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)), lty="dashed")
dev.off()

#Coefficients图形
pdf("lasso.lambda.pdf")
plot(fit, xvar="lambda", label=T)
abline(v=log(cvfit$lambda.min), lty="dashed")
dev.off()

#输出lasso显著基因表达量
coef=coef(fit, s=cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
row.names(coef)[index]
actCoef
index
lassoSigExp=rt[,c("time", "state", lassoGene)]
lassoSigExpOut=cbind(id=row.names(lassoSigExp), lassoSigExp)
write.table(lassoSigExpOut,file="7.lasso回归的基因表达矩阵和临床合并.txt",sep="\t",row.names=F,quote=F)




#4.多因素cox + 风险评分


set.seed(3296)
data <- read.table("6.单因素cox基因表达矩阵和临床数据.txt", header = T, sep = "\t", check.names=F, row.names=1, stringsAsFactors=F)
train_id=sample(1:nrow(data),0.5*nrow(data))#训练集id
train=data[train_id,]  #训练集
test=data[-train_id,]  #测试集
write.table(train, file="7.训练集表达矩阵和临床数据.txt", sep="\t", row.names=T, quote=F)
write.table(test, file="8.测试集表达矩阵和临床数据.txt", sep="\t", row.names=T, quote=F)



Train_data <- read.table("7.训练集表达矩阵和临床数据.txt", header = T, sep = "\t", check.names=F, row.names=1, stringsAsFactors=F)
Test_data <- read.table("8.测试集表达矩阵和临床数据.txt", header=T, sep="\t", check.names=F, row.names=1, stringsAsFactors = F)
rt <- Train_data
rt <- Test_data
#多因素cox
rt=read.table("7.lasso回归的基因表达矩阵和临床合并.txt",header=T,sep="\t",check.names=F,row.names=1)   
multiCox=coxph(Surv(time, state) ~ ., data = rt)
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
outTab <- data.frame(outTab)
outMultiTab<- outTab
#write.table(outMultiTab, file='多因素cox基因.txt',sep = "\t",row.names = T,col.names = NA,quote = F)

#风险及评分计算
actCoef <- round(as.numeric(outMultiTab$coef),3)  #系数提取
lsg <- rownames(outMultiTab)  #每个基因
FinalGeneExp = rt[,lsg] #每个基因的表达量
myFun = function(x){crossprod(as.numeric(x),actCoef)}  #风险评分计算函数
riskScore = apply(FinalGeneExp,1,myFun)   #风险评分
risk = as.vector(ifelse(riskScore > median(riskScore), "high", "low"))    #风险高低
trainRiskOut = cbind(rt[,c("time","state",lsg)], riskScore=as.vector(riskScore), risk)

#输出文件
coef_table1 <- outMultiTab %>% as.data.frame()
colnames(coef_table1)[1:2] <- c("Gene","Coef")
write.table(coef_table1, file='多因素cox基因.txt',sep = "\t",row.names = T,col.names = NA,quote = F)
#write.table(trainRiskOut, file='8.多因素cox基因表达矩阵和临床数据和风险评分.txt',sep = "\t",row.names = T,col.names = NA,quote = F)

Train_data <- trainRiskOut
Test_data <- trainRiskOut
dat <- rbind(Train_data,Test_data)

write.table(dat,file="9.多因素cox基因表达矩阵和临床数据和风险评分合并.txt", sep="\t", quote=F, row.names = T,col.names = T)
write.table(Train_data,file="10.多因素cox基因表达矩阵和临床数据和风险评分合并Train.txt", sep="\t", quote=F, row.names = T,col.names = T)
write.table(Test_data,file="11.多因素cox基因表达矩阵和临床数据和风险评分合并Test.txt", sep="\t", quote=F, row.names = T,col.names = T)











#预后模型评估三连
diff=survdiff(Surv(time, state) ~risk,data = trainRiskOut)
pValue=1-pchisq(diff$chisq,df=1)
fit <- survfit(Surv(time, state)~ risk, data =trainRiskOut )
pround <- round(pValue,9)

#生存曲线
p <- ggsurvplot(fit,
               data = trainRiskOut,
               palette =c("red","#3C5488FF"),
               risk.table = T,
               pval = pround)
p

#ROC曲线
ROC3 <- timeROC(T=trainRiskOut$time,   #结局时间
                delta=trainRiskOut$stat,   #结局指标
                marker=trainRiskOut$riskScore,   #预测变量
                cause=1,   #阳性结局指标数值
                weighting="marginal",   #计算方法，默认为marginal
                times=c(1,2,3,4,5),   #时间点，选取1年，3年和5年的生存率
                iid=TRUE)

plot(ROC3,
     time=1, lwd=2,col="red",main="cnm")
plot(ROC3,
     time=2, lwd=2,col="#3C5488FF", add=TRUE)   
plot(ROC3,
     time=3, lwd=2,col="grey", add=TRUE)
plot(ROC3,
     time=4, lwd=2,col="yellow", add=TRUE)   
plot(ROC3,
     time=5, lwd=2,col="green", add=TRUE)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC3$AUC[1])),
         paste0('AUC at 2 years: ',sprintf("%.03f",ROC3$AUC[2])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC3$AUC[3])),
         paste0('AUC at 4 years: ',sprintf("%.03f",ROC3$AUC[4])),
         paste0('AUC at 5 years: ',sprintf("%.03f",ROC3$AUC[5]))),
       col=c("red",'#3C5488FF','grey',"yellow","green"),lwd=2,bty = 'n')


#风险曲线

#####患者分布可视化##
rt <- trainRiskOut
rt=rt[order(rt$riskScore),]                                            
#绘制风险曲线
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
line=rt[,"riskScore"]
line[line>10]=10
pdf(file="riskScore.pdf",width = 10,height = 3.5)
plot(line, type = "p", pch = 20,
     xlab = "", ylab = "Risk score",
     col = c(rep("blue", lowLength), rep("red", highLength)),
     axes = FALSE)
axis(2)
box()
abline(h=median(rt$riskScore),v=lowLength,lty=2)
legend("topleft", c("High risk", "low Risk"),bty="n",pch=19,col=c("red","blue"),cex=1.2)
dev.off()

#绘制生存图
color=as.vector(rt$state)
color[color==1]="red"
color[color==0]="blue"
pdf(file="survstat.pdf",width = 10,height = 3.5)
plot(rt$time, pch=19,
     xlab="", ylab="Survival time (years)",
     col=color)
legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","blue"),cex=1.2)
abline(v=lowLength,lty=2)
dev.off()

#绘制风险热图
rt1=rt[c(3:(ncol(rt)-2))]
rt1=log2(rt1+1)
rt1=t(rt1)
annotation=data.frame(type=rt[,ncol(rt)])
rownames(annotation)=rownames(rt)
pdf(file="heatmap.pdf",width = 10,height = 4)
pheatmap(rt1, 
         annotation=annotation, 
         cluster_cols = FALSE,
         fontsize_row=11,
         show_colnames = F,
         fontsize_col=3,
         scale = "row",
         color = colorRampPalette(c("blue", "white", "red"))(100) )
dev.off()


#森林图绘制
bioForest=function(coxFile=null, forestFile=null, forestCol=null){
  #读取输入文件
  rt <- read.table(coxFile, header=T, sep="\t", check.names=F, row.names=1)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  #输出图形
  height=nrow(rt)/12.5+5
  pdf(file=forestFile, width=7, height=height)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #绘制森林图左边的基因信息
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
  
  #对HR进行可视化,绘制右边的森林图
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol[1], forestCol[2])
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.6)
  axis(1)
  dev.off()
}

#调用函数, 对预后相关的基因进行可视化, 得到森林图
bioForest(coxFile="多因素cox基因.txt", forestFile="多因素森林图.pdf", forestCol=c("red","blue"))


#5.临床多因素cox回归


#读入
risk=read.table("8.多因素cox基因表达矩阵和临床数据和风险评分.txt", header=T, sep="\t", check.names=F, row.names=1)#读取风险文件
cli=read.table("TCGA临床数据(清洗).txt", header=T, sep="\t", check.names=F, row.names=1)#读取临床文件

#合并
sameSample=intersect(row.names(cli),row.names(risk))
risk=risk[sameSample,]
cli=cli[sameSample,]
rt=cbind(time=risk[,1], state=risk[,2], cli, riskScore=risk[,(ncol(risk)-1)])
colnames(rt)

#设置为数值或分类变量,注意，这里不用列名的原因，
#是因为，T，不仅仅是字符的T分期，也是逻辑的true
rt[,4] = as.numeric(rt[,4])
rt[,3] = factor(rt[,3])
rt[,5] = factor(rt[,5])
rt[,6] = factor(rt[,6])
rt[,7] = factor(rt[,7])
rt[,8] = as.numeric(rt[,8])

#删除空值及NA
rt=rt[apply(rt,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
#paste(colnames(rt), collapse="+")

#cox回归模型构建
coxmod<-coxph(Surv(time, state)~.,
              data=rt)
summary(coxmod)

#可将threshold设置为0.05，既排除单因素大于0.05的
ft3 = autoReg(coxmod,uni=TRUE,threshold=0.05) 
myft(ft3)
#导出
write.table(ft3, file="result.txt", sep="\t", row.names=F, quote=F)

#paste(colnames(rt), collapse="+")
#纳入p小于0.05临床特征
multiCox<-coxph(Surv(time, state)~ gender + age + N_stage + T_stage + Stage + riskScore, #此处修改
                data=rt)
multiCoxSum = summary(multiCox)
multiTab=data.frame()
multiTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
multiTab=cbind(id=row.names(multiTab),multiTab)
write.table(multiTab, file="临床多因素cox结果.txt", sep="\t", row.names=F, quote=F)

#绘制森林图
pdf(file="forest.pdf", width=8, height=6, onefile = FALSE)
ggforest(multiCox, data = rt,
         main = "Hazard ratio",
         cpositions = c(0.02, 0.22, 0.4), 
         fontsize = 0.8, 
         refLabel = "reference", 
         noDigits = 3)
dev.off()



#5.Nomogram

#读取风险输入文件
risk=read.table("8.多因素cox基因表达矩阵和临床数据和风险评分.txt", header=T, sep="\t", check.names=F, row.names=1)
#读取临床数据文件
cli=read.table("TCGA临床数据(清洗).txt", header=T, sep="\t", check.names=F, row.names=1)
#删除包含unknow的样本
cli=cli[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
#将年龄这一列转化为数值型
cli$Age=as.numeric(cli$Age)

#合并数据我
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
colnames(risk1)
#rt=cbind(risk1[,c("time", "state", "risk")], cli)
#rt=cbind(risk1[,c("time", "state", "riskScore")], cli)
rt=cbind(time=risk1[,1], state=risk1[,2], cli, riskScore=risk1[,(ncol(risk1)-1)])
#rt=cbind(risk1, cli)

#cox回归
res.cox=coxph(Surv(time, state) ~ . , data = rt)


#绘制列线图

#单项得分，即图中的Points，表示每个变量在不同取值下所对应的单项分数，
#总得分，即Total Point，表示所有变量取值后对应的单项分数加起来合计的总得分。
#预测概率

nom1=regplot(res.cox,
             clickable=F,
             title="",
             points=TRUE,
             droplines=TRUE,
             observation=rt[50,],
             rank="sd",
             failtime = c(1,3,5),
             prfail = F)

#列线图打分
nomoRisk=predict(res.cox, data=rt, type="risk")
rt=cbind(risk1, Nomogram=nomoRisk)    #重新定义rt了
outTab=rbind(ID=colnames(rt), rt)
write.table(outTab, file="nomoRisk.txt", sep="\t", col.names=F, quote=F)


#校准曲线
pdf(file="calibration.pdf", width=6, height=6)
#1年
f <- cph(Surv(time, state) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=1)
cal <- calibrate(f, cmethod="KM", method="boot", u=1, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=3, col="Firebrick2", sub=F)
#3年
f <- cph(Surv(time, state) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=3)
cal <- calibrate(f, cmethod="KM", method="boot", u=3, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", lwd=3, col="MediumSeaGreen", sub=F, add=T)
#5年
f <- cph(Surv(time, state) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=5)
cal <- calibrate(f, cmethod="KM", method="boot", u=5, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",  lwd=3, col="NavyBlue", sub=F, add=T)
legend('bottomright', c('1-year', '3-year', '5-year'),
       col=c("Firebrick3","MediumSeaGreen","NavyBlue"), lwd=3, bty = 'n')
dev.off()


#累计风险曲线
data=cbind(time=risk1[,1], state=risk1[,2], cli, riskScore=risk1[,(ncol(risk1)-1)],Nomogram=nomoRisk)
data=data[order(data$Nomogram),]  
NomoRisk = as.vector(ifelse(data$Nomogram > median(data$Nomogram), "high", "low"))    #风险高低
data = cbind(data, NomoRisk)

fit <- survfit(Surv(time,state) ~ NomoRisk, 
               data = data) # 拟合方程

ggsurvplot(fit,
           data=data,
           conf.int=TRUE,
           risk.table=TRUE)


ggsurvplot(fit,
           fun = "cumhaz",                  #指定累计风险曲线绘制 
           legend.title="",                 #设置图例标题
           legend=c(0.15,0.90),             #指定图例位置
           conf.int = TRUE,                 #可信区间
           palette = "lancet",              #支持ggsci配色，自定义颜色，brewer palettes中的配色，等
           ggtheme = theme_bw())          #支持ggplot2及其扩展包的主题


ggsurvplot(fit,                             #创建的拟合对象
           data = data,                       #指定变量数据来源
           conf.int=TRUE,                   #显示置信区间
           pval=T,                          #添加P值
           size = 1,                        # 更改线条粗细
           surv.median.line="hv",           #添加中位生存时间线
           risk.table=TRUE,                 #添加风险表
           risk.table.col="strata",         #根据分层更改风险表颜色
           xlab="Follow up time(d)",        #指定x轴标签
           legend=c(0.85,0.85),             #指定图例位置
           legend.title="Kaplan-Meier",     #设置图例标题
           legend.labs=c("pr=1","pr=0"),    #指定图例分组标签
           risk.table.height = 0.25,        # risk table高度
           #break.x.by=10,                  #设置x轴刻度间距
           break.y.by=0.1,                  #设置y轴刻度间距
           #palette=c("#E7B800","#2E9FDF"), #更改线条颜色
           ggtheme=theme_bw())              #添加网格线


#DCA风险曲线绘制

#构建cox回归模型
cox_fit1 <- coxph(Surv(time, state) ~ age, data = data)
cox_fit2 <- coxph(Surv(time, state) ~ gender , data = data)
cox_fit3 <- coxph(Surv(time, state) ~ Stage  , data = data)
cox_fit4 <- coxph(Surv(time, state) ~ riskScore  , data = data)
cox_fit5 <- coxph(Surv(time, state) ~ Nomogram  , data = data)
cox_fit6 <- coxph(Surv(time, state) ~ T_stage  , data = data)
cox_fit7 <- coxph(Surv(time, state) ~ N_stage  , data = data)

#决策曲线分析
df1 <- ggDCA::dca(cox_fit1, cox_fit2, cox_fit3, cox_fit4, cox_fit5, times = 1)  # 1.5年，默认值是中位数

#DCA曲线绘制和美化
ggplot(df1,linetype = F,lwd = 1.2)+
  scale_color_jama(name="Model Type",labels=c("Age","Gender","Stage","riskScore","Nomogram","All","None"))+
  theme_classic(base_size = 14)+
  theme(legend.position = c(0.8,0.75),
        legend.background = element_blank())


#DCA曲线绘制和美化
ggplot(df1,linetype = F,lwd = 1.2) +
  theme_classic()+
  theme_prism(base_size =15)+
  scale_x_continuous(limits = c(0, 1),guide = "prism_minor") +
  scale_y_continuous(limits = c(-0.05, 0.15),guide = "prism_minor")+
  scale_color_manual(values = c("#F4796E", "#C09C06", "#53B400", "#22C094", "#34B4EB", "#A987FF", "#F962D7"),labels = c("Age", "Gender", "Stage", "riskScore", "Nomogram", "All", "None")) 



#临床相关性热图

#读取
data=read.table("1.TCGA原始表达矩阵.txt", header=T, sep="\t", check.names=F,row.names = 1)
#转化为matrix
dimnames=list(rownames(data), colnames(data))
data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)

#去除正常样品
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=t(data[,group==0])
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))

#根据ABCA3表达量的中位值，把样品分为两组
gene = "riskScore"
Type = ifelse(data[,gene]>quantile(data[,gene], seq(0,1,1/2))[2],"High","Low")
Type=factor(Type, levels=c("Low","High"))
data=cbind(as.data.frame(data), Type)
data=data[,c(gene,"Type")]
rownames(data) = gsub("[.]","-",rownames(data))
data=data[order(data[,gene]),] 

#读取临床
cli=read.table("clinical.txt", header=T, sep="\t", check.names=F, row.names=1)
cli[,"Age"]=ifelse(cli[,"Age"]=="unknow", "unknow", ifelse(cli[,"Age"]>65,">65","<=65"))

#合并
samSample=intersect(row.names(data), row.names(cli))
data=data[samSample,"Type",drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(data, cli)

#对临床性状进行循环，观察临床性状在高低表达组之间是否具有差异
sigVec=c(gene)
for(clinical in colnames(rt[,2:ncol(rt)])){
  data=rt[c("Type", clinical)]
  colnames(data)=c("Type", "clinical")
  data=data[(data[,"clinical"]!="unknow"),]
  tableStat=table(data)
  stat=chisq.test(tableStat)
  pvalue=stat$p.value
  Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
  sigVec=c(sigVec, paste0(clinical, Sig))
}
colnames(rt)=sigVec

#颜色
bioCol=c("#0066FF","#FF9900","#FF0000","#ed1299", "#0dbc21", "#246b93", "#cc8e12", "#d561dd", "#c93f00", 
         "#ce2523", "#f7aa5d", "#9ed84e", "#39ba30", "#6ad157", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
         "#1a918f", "#7149af", "#ff66fc", "#2927c4", "#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
         "#4aef7b", "#e86502",  "#99db27", "#e07233", "#8249aa","#cebb10", "#03827f", "#931635", "#ff523f",
         "#edd05e", "#6f25e8", "#0dbc21", "#167275", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
         "#dd27ce", "#07a301", "#ddd53e",  "#391c82", "#2baeb5","#925bea", "#09f9f5",  "#63ff4f")
colorList=list()
colorList[[gene]]=c("Low"="Orange2", "High"="DarkOrchid")
j=0
for(cli in colnames(rt[,2:ncol(rt)])){
  cliLength=length(levels(factor(rt[,cli])))
  cliCol=bioCol[(j+1):(j+cliLength)]
  j=j+cliLength
  names(cliCol)=levels(factor(rt[,cli]))
  cliCol["unknow"]="grey75"
  colorList[[cli]]=cliCol
}

#绘制热图
ha=HeatmapAnnotation(df=rt, col=colorList)
zero_row_mat=matrix(nrow=0, ncol=nrow(rt))
Hm=Heatmap(zero_row_mat, top_annotation=ha)

#输出热图
pdf(file="heatmap.pdf", width=7, height=5)
draw(Hm, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()




#临床相关性图
riskFile="risk.TCGA.txt"                         #风险文件
subtypeFile="Subtype_Immune_Model_Based.txt"     #免疫分型的文件

#定义卡方检验函数
ChischisqTest <- function(tabledata){
  ct = chisq.test(tabledata)
  ct.pvalue = ct$p.value
  ct.pvalue = ifelse(ct.pvalue<0.001,0.001,round(ct.pvalue,3))
}

#定义绘制长方形函数
rect2text <- function(x1,y1,x2,y2,text,rect.col,text.cex,text.col="white",...){
  rect(x1,y1,x2,y2,col=rect.col,border=NA,...)
  text((x1+x2)/2,(y1+y2)/2,text,cex=text.cex,col=text.col,...)
}

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
colnames(risk)[ncol(risk)]="Risk"

#读取免疫分型的文件
subtype=read.table(subtypeFile, header=T, sep="\t", check.names=F, row.names=1)
rownames(subtype)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(subtype))

#样品取交集
sameSample=intersect(row.names(subtype), row.names(risk))
subtype=subtype[sameSample,]
subtype=gsub(".+Immune |\\)","",subtype)
risk=risk[sameSample,]
data=cbind(subtype, as.data.frame(risk))
typeTab=table(data$subtype)
typeName=names(typeTab[typeTab>6])
merge.data=data[which(data[,"subtype"] %in% typeName),]
cliName=colnames(merge.data)[1]
colnames(merge.data)[1]="Clinical"

#进行卡方检验，观察免疫分型在高低风险组之间是否具有差异
groups = sort(unique(merge.data$Clinical))
groups.table = table(merge.data$Clinical)
groups.num = length(groups)
samples.num = nrow(merge.data)
groups.lownum = sum(merge.data$Risk=="Low")
groups.highnum = sum(merge.data$Risk=="High")
pergrouptable = table(merge.data$Clinical, merge.data$Risk)
ct.pvalue = ChischisqTest(pergrouptable)

#定义图形的坐标和图形的颜色
xlim1=0; xlim2=samples.num; ylim1=1; ylim2=10
space = dist.up = 0.1
space.x = 1/100*xlim2
groups.col = col = colorRampPalette(brewer.pal(9, "Paired"))(groups.num)
high_low.col = c(rgb(235/255,116/255,106/255),rgb(30/255,181/255,184/255))
bottom.unitwidth = xlim2 / (groups.num+1+1.3)
bottom.left.width = bottom.unitwidth*1.3

#输出图形
pdf(file="immSubtype.pdf", width=10, height=7)
plot(1,type="n",xlim=c(xlim1,xlim2),ylim=c(ylim1,ylim2),axes=F,xlab="",ylab="")

#绘制图形的表头(免疫分型的结果)
header.text.cex = 1.5
header.text = gettextf("%s TCGA patients",samples.num)
header.text.width = strwidth(header.text,cex=header.text.cex)
arrows.len = (xlim2 - header.text.width - space.x*2)/2
text(xlim2/2,9.5,header.text,cex=header.text.cex,font=2)
arrows(0,9.5,arrows.len,9.5,code=1,lwd=3,length=0.2)
arrows(xlim2-arrows.len,9.5,xlim2,9.5,code=2,lwd=3,length=0.2)
rect2text(x1=0,y1=4.5+space,x2=bottom.left.width,y2=7-space,text.cex=1.2,text = 
            "Risk",rect.col=rgb(64/255,109/255,181/255),text.col="white",font=2)
#draw bottom header
bottom.header.text = paste0(cliName, " group (n=" , samples.num, ")")
rect2text(x1=bottom.left.width+space.x,y1=6+space/2,x2=xlim2,y2=7-space,text.cex=1.2,
          text= bottom.header.text,rect.col=rgb(64/255,109/255,181/255),text.col="white",font=2)

#绘制图形左边的样品分组信息(低风险组和高风险组)
#draw bottom left 2
bottom.left2.text = gettextf("Low risk\n(n=%s)",groups.lownum)
rect2text(x1=0,y1=3+space/2,x2=bottom.left.width,y2=4.5-space/2,text.cex=1.2,text= 
            bottom.left2.text,rect.col=high_low.col[2],text.col="white",font=2)
#draw bottom left 3
bottom.left3.text = gettextf("High risk\n(n=%s)",groups.highnum)
rect2text(x1=0,y1=1.5+space/2,x2=bottom.left.width,y2=3-space/2,text.cex=1.2,text= 
            bottom.left3.text,rect.col=high_low.col[1],text.col="white",font=2)

#绘制免疫分型的长方形
start = 0
for(i in 1:length(groups.table)){
  end = groups.table[i]+start
  namesi = names(groups.table)[i]
  # up rect
  rect2text(x1=start,y1=8+space,x2=end,y2=9-space,text.cex=1.1,text= namesi, 
            rect.col=groups.col[i],text.col="white")
  merge.datai = merge.data[merge.data$Clinical==namesi,,drop=F]
  high.num = sum(merge.datai$Risk=="High")
  low.num = sum(merge.datai$Risk=="Low")
  # middle low rect
  rect(start,7+dist.up,start+low.num,8-dist.up,col=high_low.col[2],border=NA)
  # middle high rect
  rect(start+low.num,7+dist.up,end,8-dist.up,col=high_low.col[1],border=NA)
  #展示所有样品中免疫分型的数目
  bottom.starti = bottom.left.width+(i-1)* bottom.unitwidth
  bottom.endi = bottom.starti+bottom.unitwidth
  subheader.text = gettextf("%s\n(n=%s, %s)",namesi,nrow(merge.datai),paste0(round(nrow(merge.datai)/samples.num*100),"%"))
  rect2text(x1=bottom.starti+space.x,y1=4.5+space,x2=bottom.endi,y2=6-space,text.cex=1.1,text= subheader.text,rect.col=groups.col[i],text.col="white")
  #展示低风险组中免疫分型的数目
  low.texti = gettextf("%s(%s)",low.num,paste0(round(low.num/groups.lownum*100),"%"))
  rect2text(x1=bottom.starti+space.x,y1=3+space/2,x2=bottom.endi,y2=4.5-space/2,text.cex=1.1,text= low.texti,rect.col="grey90",text.col="black")
  #展示高风险组中免疫分型的数目
  high.texti = gettextf("%s(%s)",high.num,paste0(round(high.num/groups.highnum*100),"%"))
  rect2text(x1=bottom.starti+space.x,y1=1.5+space/2,x2=bottom.endi,y2=3-space/2,text.cex=1.1,text= high.texti,rect.col="grey70",text.col="black")
  start = end
}

#绘制pvalue的长方形
rect2text(x1=bottom.endi+space.x,y1=4.5+space,x2=xlim2,y2=6-space,text.cex=1.1,text= "P-value",rect.col="grey70",text.col="black",font=3)
rect2text(x1=bottom.endi+space.x,y1=1.5+space/2,x2=xlim2,y2=4.5-space/2,text.cex=1.1,text= ct.pvalue,rect.col="grey70",text.col="black")

#给图形加上边框
rect(0,8+space,xlim2,9-space,lwd=2,border="grey30")
rect(0,7+space,xlim2,8-space,lwd=2,border="grey30")
rect(0,1.5+space/2,xlim2,7-space,lwd=2,border="grey30")

#绘制图例
legend("bottom",legend=c("Risk-low","Risk-high"),col=rev(high_low.col),bty="n",ncol=2,pch=15,pt.cex=1.3,cex=1.3)
dev.off()







#6.免疫评估
#6.1 Cibersort

#导入R包内置数据
data(LM22)
#读入
data=read.table("1.TCGA原始表达矩阵.txt", header=T, sep="\t", check.names=F,row.names = 1)
#转化为matrix
dimnames = list(rownames(data), colnames(data))
data = matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)
#进行cibersort分析
results <- cibersort(sig_matrix = LM22, mixture_file = data)
#保存
#results=as.matrix(results[,1:(ncol(results)-3)])
results=as.matrix(results[,1:ncol(results)])
write.table(results, file="cibersort结果.txt", sep="\t", quote=F, col.names=F, row.names = T)

### 准备数据
#exp_data = read.table("1.TCGA原始表达矩阵.txt",sep="\t",header=T,row.names=1,check.names=F)
#cibersort.results = read.csv("./cibersort.results.csv",row.names = 1)

exp_data = data
cibersort.results = results
group_data = read.table("group.txt",sep="\t",header=T,row.names=1,check.names=F)
cols = c('#0072b5',"#bc2c29","#20854e","#ef7c1c")

# 添加样本列和分组列
cibersort.results$sample <- rownames(cibersort.results)
cibersort.results$group <- group_data$group
# 转化适合ggplot画图的格式
cibersort.results_melt = melt(cibersort.results)
# 修改列名
cibersort.results_melt <- cibersort.results_melt %>% dplyr::rename( Celltype = variable,Composition = value)

### 箱线图
p_boxplot  <- ggplot(cibersort.results_melt, aes(x = Celltype, y = Composition))+
                     geom_boxplot(aes(fill = group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+
                     stat_compare_means(aes(group =  group),
                     size = 5,
                     label = "p.signif",
                     # label = "p.format",
                     method = "t.test",
                     hide.ns = T,
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), symbols = c("***", "**", "*", "")))+
                     scale_fill_manual(values = cols)+
                     theme_bw()+
                     theme(axis.text.x = element_text(hjust = 1,angle = 45),
                           axis.text = element_text(colour = "black"))+
                     labs(y="Cell composition",x= NULL)

p_boxplot
ggsave("免疫细胞箱线图.pdf",width = 8,height = 5)

### 小提琴图
mean_df <- cibersort.results_melt %>% group_by(Celltype,group) %>% summarise(mean = mean(Composition))

p_violin  <-  ggplot(data = cibersort.results_melt,aes(x = Celltype,y = Composition))+
              geom_violin(width=0.6,size = 0.2,
              trim = T,
              scale = "width",
              aes(fill = group),
              position=position_dodge(0.9))+
              geom_point(data = mean_df,
              aes(x = Celltype,y = mean,group = group),
              color = "white",
              position=position_dodge(0.9),
              size = 1.3)+
              scale_y_continuous(expand = expansion(mult = c(0.03, 0.07)))+
              theme_bw()+
              theme( axis.text.x = element_text(angle = 45, hjust = 1),
                     axis.text = element_text(colour = "black"))+
              scale_fill_manual(values = cols)+
              stat_pwc(
              aes(group = group), 
              tip.length = 0,
              y.position = 0.55,    # 尽量让p值跟提琴图有距离
              method = "t_test", label = "p.format",
              bracket.nudge.y = -0.08)

p_violin
ggsave("免疫细胞提琴图.pdf",width = 14,height = 5)

### 免疫细胞丰度堆叠图

# 设置样本的因子水平，按分组排序
group_data <- group_data %>% arrange(group)
cibersort.results_melt$sample <- factor(cibersort.results_melt$sample,levels =rownames(group_data))

# 颜色
mypalette <- colorRampPalette(brewer.pal(8,"Paired"))

p_barplot = ggplot(cibersort.results_melt,aes(sample, Composition,fill = Celltype)) +
                   geom_bar(stat = "identity") +
                   labs(fill = "Cell Type",x = "",y = "Estiamted Proportion") +
                   theme_bw() +
                   theme(axis.text.x = element_blank(),
                         axis.ticks.x = element_blank()) +
                   scale_y_continuous(expand = c(0.01,0)) +
                   scale_fill_manual(values = mypalette(22))+
                   theme(legend.position = "right")+
                   guides(fill = guide_legend(ncol = 1) )+
                   theme(plot.margin = margin(b = 0,unit =  "pt"))

p_barplot

### group_data 用于绘制下方的注释
group_data$group <- fct_inorder(group_data$group)
group_data$sample <- rownames(group_data)
group_data$value <- 1
group_data$sample <- factor(group_data$sample,levels =rownames(group_data))

### label_data用来添加分组文字
label_data <- group_data %>% group_by(group) %>% summarise(n = n())
# i = 2
for (i in 1:nrow(label_data)) {
  pre_n <- 0
  if (i == 1) {
    label_data[i,"pos"] <- label_data[i,"n"] %>% pull(n)/2
  }else{
    pre_n <- sum(label_data[1:(i-1),"n",drop = T]) + pre_n
    
    label_data[i,"pos"] <- label_data[i,"n"] %>% pull(n)/2+pre_n
  }
}

label_data$group
label_data$value <- c(1,1)


p <- ggplot(group_data,aes(x=sample,y=value))+
     geom_tile(aes(fill=group))+
     geom_text(data = label_data,inherit.aes = F,
               aes(x = pos,y = value,label = group),size = 6,nudge_x = 0.5)+ 
     scale_x_discrete(expand = c(0,0))+
     theme(panel.background = element_blank(),
           axis.ticks = element_blank(),
           axis.title = element_blank(),
           axis.text.x = element_blank(),
           axis.text.y = element_blank(),
           legend.position = "none",
           legend.title = element_blank())+
     scale_fill_manual(values = c("#0072b5","#ef7c1c"))+
     theme(plot.margin = margin(t = 0,unit =  "pt"))
p

library(aplot)  
p_barplot %>% insert_bottom(p,height = 0.08)
ggsave(filename = "./免疫细胞丰度堆叠条形图.pdf",width = 9,height = 6)

### 免疫细胞之间的相关性热图

# 加载计算相关性的函数
source("cal_cor.R")

cibersort.results2 <- cibersort.results[,1:22]
cibersort.results2 <- cibersort.results2[,!colSums(cibersort.results2) == 0]
cor_res <- cal_cor(df_x = cibersort.results2,df_y = cibersort.results2,cor_type = "pearson")

cor_mat <- cor_res$r
cell_corr_plot <- ggcorrplot(corr = cor_mat,show.legend = T,colors = c("#2166AC","white","#B2182B"),hc.order = T,lab = T)+coord_cartesian(clip = "off")
cell_corr_plot

ggsave(filename = "./免疫细胞之间的相关性热图.pdf",width = 12,height = 10)



###风险评分和细胞相关性
# 选择基因(风险屏风)和细胞
genes = c("riskScore")
cell_type = colnames(Macrophages.M0) 
#cell_type = colnames(cibersort.results)      # 这里用全部的 也可以只选择感兴趣的来绘制

### 单个基因与免疫细胞的相关性散点图带侧边密度柱状图
dir.create("单个基因和免疫细胞的相关性散点图")

for (gene in genes) {
  p_list <- list()
  for(i in 1:length(cell_type)){
    gene_exp = exp_data[gene,]
    gene_exp <- gene_exp %>% t() %>% as.data.frame()
    #gene_exp$riskScore <- as.numeric(gene_exp$riskScore)
    gene_cell_df <- cbind(gene_exp,cibersort.results[,cell_type[i],drop = F])
    colname <- colnames(gene_cell_df)
    pearson_R <- cor(gene_cell_df[,1],gene_cell_df[,2],method = "pearson")
    pearson_P <- cor.test(gene_cell_df[,1],gene_cell_df[,2],method = "pearson")
    pearson_P <- pearson_P$p.value
    # 散点图
    p <- ggplot(data=gene_cell_df,aes(x=.data[[colname[1]]],
                                      y=.data[[colname[2]]]))+
    geom_point(size = 2)+
    geom_smooth(method = "lm",linewidth = 1.3,
                formula = "y~x",)+
    theme_bw()
    # 添加密度  
    p <- ggExtra::ggMarginal(p, type = "density",
                             # type = "density",  #指定添加类型
                             xparams=list(fill = "#1d4a9f",color = "black",linewidth = 0.5),  #指定颜色
                             yparams = list(fill="#f4a01a",color = "black",linewidth = 0.5),  )#指定颜色
    
    # 在坐标轴绝对添加p值 cor值
    p <- cowplot::ggdraw(p) + 
      cowplot::draw_text(paste0("R= ",pearson_R,"\nP= ",pearson_P), x = 0.15, y = 0.75, hjust = 0, vjust = 0, size = 10)
    p
    p_list[[i]] <- p
  }
  pdf(paste0("./单个基因和免疫细胞的相关性散点图/",gene,".pdf"),width = 15,height = 32)
  print(ggpubr::ggarrange(plotlist = p_list,ncol = 3,nrow = ceiling(length(cell_type)/3)))
  dev.off()
}

#6.2 ssGSEA

#1.细胞比例分析
#读取矩阵文件
data=read.table("1.TCGA原始表达矩阵.txt", header=T, sep="\t", check.names=F,row.names = 1)
#转化为matrix
dimnames=list(rownames(data), colnames(data))
mat=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)
#读取数据集文件
geneSet=getGmt("immune1.gmt", geneIdType=SymbolIdentifier())
#ssGSEA分析
ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
#定义ssGSEA score矫正函数
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
#对ssGSEA score进行矫正
rt=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(rt), rt)
write.table(ssgseaOut, file="1.ssGSEA免疫细胞评分.txt", sep="\t", quote=F, col.names=F)

#去除正常样品
group=sapply(strsplit(colnames(rt),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
rt=t(rt[,group==0])
rownames(rt)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(rt))

#读入
risk=read.table("风险评分.txt",header=T,sep="\t",row.names=1,check.names=F)

#合并
sameSample=intersect(row.names(rt),row.names(risk))
rt=rt[sameSample,,drop=F]
risk=risk[sameSample,"risk",drop=F]
rt1=cbind(rt, risk)

#画图  高低分组
rt=melt(rt1,id.vars=c("risk"))
colnames(rt)=c("Risk","Type","Score")
rt$Risk=factor(rt$Risk, levels=c("low","high"))
p=ggboxplot(rt, x="Type", y="Score", color = "Risk",
            ylab="Score",add = "none",xlab="",palette = c("Orange2","DarkOrchid") )
p=p+rotate_x_text(50)
p=p+stat_compare_means(aes(group=Risk),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                                        symbols = c("***", "**", "*", "")),label = "p.signif")
pdf(file="1.免疫细胞比例.pdf", width=8, height=5)
print(p)
dev.off()


#2.功能分析
#读取矩阵文件
data=read.table("TCGA_LUSC_TPM.txt", header=T, sep="\t", check.names=F,row.names = 1)
#转化为matrix
dimnames=list(rownames(data), colnames(data))
mat=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)

#读取数据集文件
geneSet=getGmt("immune2.gmt", geneIdType=SymbolIdentifier())

#ssgsea分析
ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
#定义ssGSEA score矫正函数
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
#对ssGSEA score进行矫正
rt=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(rt), rt)
write.table(ssgseaOut, file="2.ssGSEA免疫功能评分.txt", sep="\t", quote=F, col.names=F)

#去除正常样品
group=sapply(strsplit(colnames(rt),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
rt=t(rt[,group==0])
rownames(rt)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(rt))

#读入
risk=read.table("一致性样本分类.txt",header=T,sep="\t",row.names=1,check.names=F)

#合并
sameSample=intersect(row.names(rt),row.names(risk))
rt=rt[sameSample,,drop=F]
risk=risk[sameSample,"risk",drop=F]
rt1=cbind(rt, risk)

#画图  高低分组
rt=melt(rt1,id.vars=c("risk"))
colnames(rt)=c("Risk","Type","Score")
rt$Risk=factor(rt$Risk, levels=c("low","high"))
p=ggboxplot(rt, x="Type", y="Score", color = "Risk",
            ylab="Score",add = "none",xlab="",palette = c("Orange2","DarkOrchid") )
p=p+rotate_x_text(50)
p=p+stat_compare_means(aes(group=Risk),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                                        symbols = c("***", "**", "*", "")),label = "p.signif")
pdf(file="2.免疫功能评分.pdf", width=8, height=5)
print(p)
dev.off()

#画图  AB分组
rt=melt(rt1,id.vars=c("cluster"))
colnames(rt)=c("Cluster","Type","Score")
rt$Cluster=factor(rt$Cluster, levels=c("A","B"))
p=ggboxplot(rt, x="Type", y="Score", color = "Cluster",
            ylab="Score",add = "none",xlab="",palette = c("red","blue") )
p=p+rotate_x_text(50)
p=p+stat_compare_means(aes(group=Cluster),
                       symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("***", "**", "*", "")),
                       label = "p.signif")

pdf(file="ssGSVE富集功能分析.pdf", width=8, height=5)
print(p)
dev.off()



#6.3 Estimate

#基因过滤，将表达谱文件中的基因和软件的common gene进行比对过滤
filterCommonGenes(input.f="1.TCGA原始表达矩阵.txt", output.f="commonGenes.gct", id="GeneSymbol")
#运行estimate包
estimateScore(input.ds = "commonGenes.gct",output.ds="estimateScore.gct")
#输出每个样品的打分
scores=read.table("estimateScore.gct", skip=2, header=T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)=gsub("\\.", "\\-", rownames(scores))
out=rbind(ID=colnames(scores), scores)
write.table(out, file="ESTIMATE评分.txt", sep="\t", quote=F, col.names=F)
#读取风险文件
Risk=read.table("risk.txt", header=T, sep="\t", check.names=F, row.names=1)
Risk$risk=factor(Risk$risk, levels=c("low","high"))
#仅保留前3列,TumorPurity的值与其它差距太大，且与Estimate socre相反,可单独画图
score=scores[,1:3]

#正常和肿瘤数目,第14,15字符,01-09是癌症，10-19是正常，20-29是癌旁
group=sapply(strsplit(rownames(score),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
#保留肿瘤样本
score = score[group==0,]
#样本名仅保留前12字符
rownames(score)=substr(rownames(score),1,12)
#将.改为-
rownames(score) = gsub('[.]', '-', rownames(score))
#获取共同样本
score=score[row.names(Risk),,drop=F]
#数据合并
rt=cbind(Risk[,"risk",drop=F], score)

#转化
data=melt(rt, id.vars=c("risk"))
colnames(data)=c("Risk", "scoreType", "Score")

#小提琴图
p=ggviolin(data, x="scoreType", y="Score", fill = "Risk",
           xlab="",
           ylab="TME score",
           legend.title="Risk",
           add = "boxplot", add.params = list(color="white"),
           palette = c("SteelBlue3","Firebrick3","green"), width=1)
p=p+rotate_x_text(45)
p1=p+stat_compare_means(aes(group=Risk),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")

#输出图形
pdf(file="ESTIMATE高低风险组小提琴图.pdf", width=8, height=5)
print(p1)
dev.off()

###风险评分和细胞相关性
# 选择基因(风险屏风)和细胞
genes = c("riskScore")        #加在表达矩阵里面,看作一个基因
cell_type = "ImmuneScore"     #本来这块是cibersort矩阵里面的细胞名称，
                              #现在换成estimate的两个列名，读取estimate结果文件。

### 单个基因与免疫细胞的相关性散点图带侧边密度柱状图
dir.create("单个基因和免疫细胞的相关性散点图")

for (gene in genes) {
  p_list <- list()
  for(i in 1:length(cell_type)){
    gene_exp = exp_data[gene,]
    gene_exp <- gene_exp %>% t() %>% as.data.frame()
    #gene_exp$riskScore <- as.numeric(gene_exp$riskScore)
    gene_cell_df <- cbind(gene_exp,cibersort.results[,cell_type[i],drop = F])
    colname <- colnames(gene_cell_df)
    pearson_R <- cor(gene_cell_df[,1],gene_cell_df[,2],method = "pearson")
    pearson_P <- cor.test(gene_cell_df[,1],gene_cell_df[,2],method = "pearson")
    pearson_P <- pearson_P$p.value
    # 散点图
    p <- ggplot(data=gene_cell_df,aes(x=.data[[colname[1]]],
                                      y=.data[[colname[2]]]))+
      geom_point(size = 2)+
      geom_smooth(method = "lm",linewidth = 1.3,
                  formula = "y~x",)+
      theme_bw()
    # 添加密度  
    p <- ggExtra::ggMarginal(p, type = "density",
                             # type = "density",  #指定添加类型
                             xparams=list(fill = "#1d4a9f",color = "black",linewidth = 0.5),  #指定颜色
                             yparams = list(fill="#f4a01a",color = "black",linewidth = 0.5),  )#指定颜色
    
    # 在坐标轴绝对添加p值 cor值
    p <- cowplot::ggdraw(p) + 
      cowplot::draw_text(paste0("R= ",pearson_R,"\nP= ",pearson_P), x = 0.15, y = 0.75, hjust = 0, vjust = 0, size = 10)
    p
    p_list[[i]] <- p
  }
  pdf(paste0("./单个基因和免疫细胞的相关性散点图/",gene,".pdf"),width = 15,height = 32)
  print(ggpubr::ggarrange(plotlist = p_list,ncol = 3,nrow = ceiling(length(cell_type)/3)))
  dev.off()
}


#6.4 免疫检查点相关性分析
#说明: 
#1.在原始表达矩阵中手动提取相关检查点基因，把cibersort表格替换。
#2.group不变。
#3.genes为筛选出来的基因,原始矩阵需要有riskScore

exp_data = read.table("1.TCGA原始表达矩阵.txt", header=T, sep="\t", check.names=F,row.names = 1)
cibersort.results =  read.csv("./cibersort.results.csv",row.names = 1)
group_data = read.table("group.txt",sep="\t",header=T,row.names=1,check.names=F)
cols = c('#0072b5',"#bc2c29","#20854e","#ef7c1c")



exp_data =  read.csv("./data/expdata.csv",row.names = 1,check.names = F)
group_data = read.csv("./data/group.csv",row.names = 1)
cibersort.results = read.csv("./cibersort.results.csv",row.names = 1)

genes = c("A","B","C","D","E","F","risk")




###-----------------------------------------------------------------------------
### 计算免疫细胞与基因的相关性和R方
###-----------------------------------------------------------------------------

gene_exp = exp_data[genes,]
data_x = t(gene_exp) %>% as.data.frame() %>% mutate_all(as.numeric)
data_y = cibersort.results

# 计算相关性
r_mat <- cor(x = data_x,y = data_y,method = "pearson") %>% as.data.frame()

# 计算pvalue
p_mat <- map_df(data_x ,function(x){
  map(data_y,function(y){
    cor.test(as.numeric(x),method = "pearson",
             as.numeric(y))$p.value
  }
  )
}
) %>% as.data.frame() %>% `rownames<-`(colnames(data_x))

# 计算拟合R方
rsq_mat <- map_df(data_x ,function(x){
  map(data_y,function(y){
    data = data.frame(x = as.numeric(x),y = as.numeric(y))
    model = lm(x~y, data = data)
    r.squared = summary(model)$r.squared
    r.squared
  }
  )
}
) %>% as.data.frame() %>% `rownames<-`(colnames(data_x))


# 合并信息
d1 = r_mat %>% rownames_to_column() %>% dplyr::rename(gene_symbol = rowname ) %>% melt()
d2 = p_mat %>% rownames_to_column() %>% dplyr::rename(gene_symbol = rowname ) %>% melt()
d3 = rsq_mat %>% rownames_to_column() %>% dplyr::rename(gene_symbol = rowname ) %>% melt()
df <- inner_join(d1,d2,by = c("gene_symbol","variable")) %>% 
  inner_join(d3,by = c("gene_symbol","variable"))
colnames(df) <-c("gene_symbol","Celltype","correlation","P.value","R.squared")
df
write.csv(df,file = "免疫细胞与基因的相关性表格.csv",row.names = F)


###-----------------------------------------------------------------------------
### 免疫细胞与基因的相关性热图 热图格子展示相关系数
###-----------------------------------------------------------------------------
p_corr_plot <- 
  ggcorrplot(corr = r_mat,show.legend = T,colors = c("#2166AC","white","#B2182B"),
             lab = T)+
  coord_cartesian(clip = "off")
p_corr_plot
ggsave(filename = "./免疫细胞与基因的相关性热图_展示相关系数.pdf",width = 8,height = 6)




###-----------------------------------------------------------------------------
### 免疫细胞与基因的相关性热图 热图格子展示显著性
###-----------------------------------------------------------------------------

#添加显著性**
df  <-  
  df %>% mutate(star = case_when(P.value > 0.05 ~ "",
                                 P.value < 0.05 & P.value > 0.01 ~ "*",
                                 P.value < 0.01 & P.value > 0.001 ~ "**",
                                 P.value < 0.001 ~ "***")
  )

p_corr_plot_showsig <- 
  ggplot(df, aes(gene_symbol, Celltype, fill=correlation)) +
  geom_tile(height=0.95, width=0.95) +
  scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B",limits = c(-1,1)) +
  geom_text(aes(label = star),size = 5)+
  theme_minimal() +
  coord_equal() +
  labs(x="",y="",fill="Corr") +
  theme(axis.text.x=element_text(size=13, angle=45, vjust=1, hjust=1, 
                                 margin=margin(-3,0,0,0)),
        axis.text.y=element_text(size=13, margin=margin(0,-3,0,0)),
        # legend.key.height = grid::unit(1,"lines"),
        # legend.key.width  = grid::unit(1,"lines"),
        panel.grid.major=element_blank())+
  coord_cartesian(clip = "off")

p_corr_plot_showsig
ggsave(filename = "./免疫细胞与基因的相关性热图_展示显著性.pdf",width = 6,height = 6)




#6.5 TIDE

#准备TIDE输入文件
#读取表达文件，并对输入文件整理
rt=read.table("1.TCGA原始表达矩阵.txt", header=T, sep="\t", check.names=F)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]

#标准化，按照行计算平均值，再拿这一行所有表达量-该平均值。
#既基因的均值为0
Expr <- t(apply(exp, 1, function(x)x-(mean(x)))) 
write.table(data.frame(ID=rownames(Expr),Expr,check.names = F),'tcga_normalize.txt', sep="\t", quote=F, row.names = TRUE)


#读取TIDE数据
tide=read.table("TIDE.txt", header=T, sep="\t", check.names=F, row.names=1)
trde<-read.csv("export.csv",header=T)
#仅仅保留肿瘤样本
group=sapply(strsplit(row.names(tide),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
tide=tide[group==0,,drop=F]
#修改样本名
tide[,2] = gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(tide))
#去除重复
tide1=tide[!duplicated(tide$V2),]
#整理
rownames(tide1) = tide1[,2]
tide = as.data.frame(tide1[,-2])
rownames(tide) = rownames(tide1)
colnames(tide) = "TIDE"
tide=avereps(tide)

#读取风险数据文件
risk=read.table("risk.txt", header=T, sep="\t", check.names=F, row.names=1)

#合并数据
sameSample=intersect(row.names(tide), row.names(risk))
tide=tide[sameSample, , drop=F]
risk=risk[sameSample, "risk", drop=F]
data=cbind(tide, risk)

#设置比较组
data$risk=ifelse(data$risk=="high", "High-risk", "Low-risk")
group=levels(factor(data$risk))

#1.TIDE小提琴图：
my_comparisons <- list( c("Low_Risk", "High_Risk")) #添加比较分组
p1 <- ggviolin(result, x = 'Risk', y = 'TIDE', fill = 'Risk',
               palette = c("#4979b6","#d9352a"),
               add = 'boxplot', add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", bracket.size=0.5, tip.length = 0.02, method = 't.test')
p1
#2.Dysfunction小提琴图：
p2 <- ggviolin(result, x = 'Risk', y = 'Dysfunction', fill = 'Risk',
               palette = c("#4979b6","#d9352a"),
               add = 'boxplot', add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", bracket.size=0.5, tip.length = 0.02, method = 't.test')
p2
#Exclusion小提琴图：
p3 <- ggviolin(result, x = 'Risk', y = 'Exclusion', fill = 'Risk',
               palette = c("#4979b6","#d9352a"),
               add = 'boxplot', add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", bracket.size=0.5, tip.length = 0.02, method = 't.test')
p3
#MSI小提琴图：
colnames(result)[6]
colnames(result)[6] <- c('MSI') #简化一下列名
p4 <- ggviolin(result, x = 'Risk', y = 'MSI', fill = 'Risk',
               palette = c("#4979b6","#d9352a"),
               add = 'boxplot', add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", bracket.size=0.5, tip.length = 0.02, method = 't.test')
p4


#画图
pdf(file="TIDE.pdf", width=5, height=4.5)
ggviolin(data, x="risk", y="TIDE", fill = "risk", 
         xlab="", ylab="TIDE",
         palette=c("Firebrick2","DodgerBlue1"),
         legend.title="Risk",
         add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(comparisons = my_comparisons)
dev.off()




#6.6 TCIA

#免疫打分
tcia=read.table("TCIA-ClinicalData.tsv", header=T, sep="\t", check.names=F, row.names=1)
tcia=tcia[,c("ips_ctla4_neg_pd1_neg","ips_ctla4_neg_pd1_pos","ips_ctla4_pos_pd1_neg","ips_ctla4_pos_pd1_pos")]

#读取文件
data=read.table("TCGA_LUSC_TPM.txt", header=T, sep="\t", check.names=F, row.names=1)

#正常和肿瘤数目,第14,15字符,01-09是癌症，10-19是正常，20-29是癌旁
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
#仅保留肿瘤样本的
data = data[,group == 0]
#前12个字符
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*)", "\\1\\-\\2\\-\\3", colnames(data))
#转置
data=t(data)
#合并
sameSample=intersect(row.names(tcia), row.names(data))
tcia=tcia[sameSample, , drop=F]
colnames(data)
#以CTLA4为例
data=data[sameSample, "CTLA4", drop=F]
data = cbind(tcia,data)

#设置比较组
data$group = ifelse(data[,"CTLA4"]>quantile(data[,"CTLA4"], seq(0,1,1/2))[2],"high","low")
data = data[,-5]
rt=melt(data,id.vars=c("group"))
colnames(rt)=c("group","Genesets","Expression")
group=levels(factor(rt$group))
rt$group=factor(rt$group, levels=c("low","high"))

#提琴图
violin=ggviolin(rt, x="Genesets", y="Expression", fill="group",
                xlab="",
                ylab="Score",#修改Y轴标题名称
                legend.title="group",
                #width=0.8,
                palette = c("DodgerBlue1","Firebrick2") )+#修改颜色
  rotate_x_text(50)+
  stat_compare_means(aes(group=group),
                     method="wilcox.test",
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")+
  theme(axis.text= element_text(face = "bold.italic",colour = "#441718",size = 16),
        axis.title = element_text(face = "bold.italic",colour = "#441718",size = 16),#Font face ("plain", "italic", "bold", "bold.italic")
        axis.line = element_blank(),
        plot.title = element_text(face = "bold.italic",colour = "#441718",size = 16),
        legend.text = element_text(face ="bold.italic"),
        panel.border = element_rect(fill=NA,color="#35A79D",size=1.5,linetype="solid"),
        panel.background = element_rect(fill = "#F1F6FC"),
        panel.grid.major = element_line(color = "#CFD3D6", size =.5,linetype ="dotdash" ),
        legend.title = element_text(face ="bold.italic",size = 13)
  )

#输出
pdf(file="ggviolin.pdf", width=8, height=7)
print(violin)
dev.off()



#7. GO KEGG富集分析
#https://www.bilibili.com/video/BV1V24y14732/?vd_source=143606cb3caa8f5a10d942320b27b3b5


#8.一致性聚类

#8.1 聚类分组
#读入
data=read.table("1.TCGA原始表达矩阵", header=T, sep="\t", check.names=F,row.names = 1)

#转化为matrix
dimnames=list(rownames(data), colnames(data))
data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)

#正常和肿瘤数目,第14,15字符,01-09是癌症，10-19是正常，20-29是癌旁
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
#根据正常和肿瘤排序
data = data[,group == 0]

#读取多因素cox基因的表达矩阵
multicoxgene=read.table("单因素cox基因.txt", header=T, sep="\t", check.names=F,row.names = 1)
data = data[rownames(multicoxgene),]
sameSample=intersect(row.names(multicoxgene), row.names(data))
data=data[sameSample,]


#对样品进行分型
#maxK, 最大的K值，形成一系列梯度
#pItem, 选择80%的样本进行重复抽样
#pfeature, 选择80%的基因进行重复抽样
#reps, 重复抽样的数目，先设置为100，结果不错再设置为1000
#clusterAlg, 聚类的算法,hc,pam,km
#distanc, 距离矩阵的算法,pearson,spearman,euclidean
#title, 输出结果的文件夹名字，包含了输出的图片
#seed, 随机种子，用于固定结果

results=ConsensusClusterPlus(data,
                             maxK=9,
                             reps=100,
                             pItem=0.8,
                             pFeature=1,
                             title=workDir,
                             clusterAlg="pam",
                             distance="euclidean",
                             seed=123,
                             plot="png")

#输出分型结果
clusterNum=2      #分成几个亚型
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("cluster")
letter=c("A","B","C","D","E","F","G","H","I","J")
uniqClu=levels(factor(cluster$cluster))
cluster$cluster=letter[match(cluster$cluster, uniqClu)]
clusterOut=rbind(ID=colnames(cluster), cluster)
write.table(clusterOut, file="一致性样本分类.txt", sep="\t", quote=F, col.names=F)
RunPCA(cluster)


#8.2 不同组别ssGSEA分析


#1.细胞比例分析
#读取矩阵文件
data=read.table("1.TCGA原始表达矩阵.txt", header=T, sep="\t", check.names=F,row.names = 1)
#转化为matrix
dimnames=list(rownames(data), colnames(data))
mat=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)
#读取数据集文件
geneSet=getGmt("immune1.gmt", geneIdType=SymbolIdentifier())
#ssGSEA分析
ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
#定义ssGSEA score矫正函数
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
#对ssGSEA score进行矫正
rt=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(rt), rt)
write.table(ssgseaOut, file="1.ssGSEA免疫细胞评分.txt", sep="\t", quote=F, col.names=F)

#去除正常样品
group=sapply(strsplit(colnames(rt),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
rt=t(rt[,group==0])
rownames(rt)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(rt))

#读入
risk=read.table("风险评分.txt",header=T,sep="\t",row.names=1,check.names=F)

#合并
sameSample=intersect(row.names(rt),row.names(risk))
rt=rt[sameSample,,drop=F]
risk=risk[sameSample,"risk",drop=F]
rt1=cbind(rt, risk)


#画图  AB分组
rt=melt(rt1,id.vars=c("cluster"))
colnames(rt)=c("Cluster","Type","Score")
rt$Cluster=factor(rt$Cluster, levels=c("A","B"))
p=ggboxplot(rt, x="Type", y="Score", color = "Cluster",
            ylab="Score",add = "none",xlab="",palette = c("red","blue") )
p=p+rotate_x_text(50)
p=p+stat_compare_means(aes(group=Cluster),
                       symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("***", "**", "*", "")),
                       label = "p.signif")

pdf(file="AB组间免疫细胞比例分析.pdf", width=8, height=5)
print(p)
dev.off()

#2.功能分析
#读取矩阵文件
data=read.table("1.TCGA原始表达矩阵.txt", header=T, sep="\t", check.names=F,row.names = 1)
#转化为matrix
dimnames=list(rownames(data), colnames(data))
mat=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)

#读取数据集文件
geneSet=getGmt("immune2.gmt", geneIdType=SymbolIdentifier())

#ssgsea分析
ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
#定义ssGSEA score矫正函数
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
#对ssGSEA score进行矫正
rt=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(rt), rt)
write.table(ssgseaOut, file="2.ssGSEA免疫功能评分.txt", sep="\t", quote=F, col.names=F)

#去除正常样品
group=sapply(strsplit(colnames(rt),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
rt=t(rt[,group==0])
rownames(rt)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(rt))

#读入
risk=read.table("一致性样本分类.txt",header=T,sep="\t",row.names=1,check.names=F)

#合并
sameSample=intersect(row.names(rt),row.names(risk))
rt=rt[sameSample,,drop=F]
risk=risk[sameSample,"risk",drop=F]
rt1=cbind(rt, risk)

#画图  AB分组
rt=melt(rt1,id.vars=c("cluster"))
colnames(rt)=c("Cluster","Type","Score")
rt$Cluster=factor(rt$Cluster, levels=c("A","B"))
p=ggboxplot(rt, x="Type", y="Score", color = "Cluster",
            ylab="Score",add = "none",xlab="",palette = c("red","blue") )
p=p+rotate_x_text(50)
p=p+stat_compare_means(aes(group=Cluster),
                       symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("***", "**", "*", "")),
                       label = "p.signif")

pdf(file="AB组间免疫细胞功能分析.pdf", width=8, height=5)
print(p)
dev.off()


#3.两组间临床相关性分析
#https://www.bilibili.com/video/BV1GF41127xb/?vd_source=143606cb3caa8f5a10d942320b27b3b5


#4.两组间分险得分比较箱线图
library(ggpubr)
rt$ID <- rownames(rt)
rt <- rt[,c(3,1,2)]

colnames(rt)=c("id","Cluster","RiskScore")
x=colnames(rt)[2]
y=colnames(rt)[3]

group=levels(factor(rt$Cluster))
rt$Cluster=factor(rt$Cluster, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

boxplot=ggboxplot(rt, x="Cluster", y="RiskScore", color="Cluster",
                  xlab=x,
                  ylab=y,
                  legend.title=x,
                  palette = c("blue","red"),
                  add = "jitter")+ 
stat_compare_means(comparisons = my_comparisons)





