getEigenGeneValues<-function(datRef,colorh1,datAll)
{
  eigenGenesCoef<-list()
  i<-0
  for (c in unique(colorh1))
  {
    i<-i+1
    eigenGenesCoef[[i]]<-prcomp(scale(datRef[,which(colorh1 == c)]))$rotation[,1]
  }
  names(eigenGenesCoef)<-unique(colorh1)
  values<-NULL
  for( c in unique(colorh1))
  {
    v<-rbind(datAll)[,which(colorh1 == c)] %*%  eigenGenesCoef[[c]]
    values<-cbind(values,sign(mean(v))*v)
  }
  colnames(values)<-unique(colorh1)
  values
}

plot_C_AA_Heatmap<-function(colorh1_C_AA,AdjMatC,AdjMatAA, expr_caucasian_cpmlog, expr_african_american_cpmlog,ordering=NULL,file="AA_C_NFD_DifferentialPlot.png")
{
  if (is.null(ordering))
  {
    h<-hclust(as.dist(1-abs(cor(getEigenGeneValues(expr_caucasian_cpmlog[,which(colorh1_C_AA!="grey")],colorh1_C_AA[which(colorh1_C_AA!="grey")],
                                                   rbind(expr_caucasian_cpmlog,expr_african_american_cpmlog)[,which(colorh1_C_AA!="grey")])))))
    for (c in h$label[h$order])
    {
      ordering<-c(ordering,which(colorh1_C_AA ==c))
    }
  }
  mat_tmp<-(AdjMatC[ordering,ordering])
  mat_tmp[which(row(mat_tmp)>col(mat_tmp))]<-(AdjMatAA[ordering,ordering][which(row(mat_tmp)>col(mat_tmp))])
  diag(mat_tmp)<-0
  mat_tmp<-sign(mat_tmp)*abs(mat_tmp)^(1/2)
  png(file=file,height=1000,width=1000)
  image(mat_tmp,col=rev(brewer.pal(11,"RdYlBu")),axes=F,asp=1,breaks=seq(-1,1,length.out=12))
  dev.off()
  unique(colorh1_C_AA[ordering])
}

plotExprChange<-function(datC1,datC2, colorh1_C_AA,ordering=NULL)
{
  if (is.null(ordering))
  {
    h<-hclust(as.dist(1-abs(cor(getEigenGeneValues(datC1[,which(colorh1_C_AA!="grey")],colorh1_C_AA[which(colorh1_C_AA!="grey")],rbind(datC1,datC2)[,which(colorh1_C_AA!="grey")])))))
    for (c in h$label[h$order])
    {
      ordering<-c(ordering,which(colorh1_C_AA ==c))
    }
  }
  mycolors<-colorh1_C_AA[ordering]
  plot(x=0:length(which(mycolors!="grey")),y=rep(1,length(which(mycolors!="grey"))+1),col="white",axes=F,xlab="",ylab="",ylim=c(0,1))
  rr=c(244,239,225,215,209,193,181,166,151,130,110)
  gg=c(228,204,174,160,146,117,94,58,44,45,45)
  bb=c(176,140,109,105,102,91,84,74,70,68,66)
  MyColours<-NULL
  for ( i in 1:11)
  {
    MyColours=c(MyColours,rgb(rr[i],gg[i],bb[i],maxColorValue=255)  )
  }
  exprDiff<-NULL
  l<-0
  for (c in setdiff(unique(mycolors),"grey"))
  {
    meanC1<-mean(t(datC1)[colnames(datC1)[which(colorh1_C_AA == c)],])
    meanC2<-mean(t(datC2)[colnames(datC2)[which(colorh1_C_AA == c)],])
    exprDiff<-rbind(exprDiff,c(meanC1,meanC2))
    r<-l+length(which(mycolors==c))
    rect(l,0.85,r,1,col=c,border=F)
    rect(l,0,r,.4,col=MyColours[floor(meanC2*2)-10],border="white",lwd=2)
    rect(l,0.4,r,.8,col=MyColours[floor(meanC1*2)-10],border="white",lwd=2)
    l<-r
  }
  exprDiff
}


extractModules<-function(colorh1,datExpr,anno,write=F,file_prefix="",dir=NULL)
{
  module<-list()
  if (!is.null(dir))
  {
    dir.create(dir)
    file_prefix=paste(dir,"/",file_prefix,sep="")
  }
  i<-1
  for (c in unique(colorh1))
  {
    module[[i]]<-(anno[colnames(datExpr)[which(colorh1==c)],1])
    if (write) {write.table(rownames(anno)[which(colorh1==c)],file=paste(file_prefix,"_",c,".txt",sep=""),quote=F,row.names=F,col.names=F)}
    i<-i+1
  }
  names(module)<-unique(colorh1)
  module
}


library(WGCNA)
library(edgeR)
library(flashClust)
library(RColorBrewer)

expr_caucasian <- read.csv("~/MSc/project/DiffCoEx/C_NFD_matrix.csv",row.names = 1)
expr_african_american <- read.csv("~/MSc/project/DiffCoEx/AA_NFD_matrix.csv", row.names = 1)

common_genes <- intersect(rownames(expr_caucasian), rownames(expr_african_american))
expr_caucasian <- expr_caucasian[common_genes, ]
expr_african_american <- expr_african_american[common_genes, ]

expr_caucasian_cpmlog <- cpm(expr_caucasian,log=TRUE)
expr_african_american_cpmlog <- cpm(expr_african_american,log=TRUE)

#transposal
expr_caucasian_cpmlog <- t(expr_caucasian_cpmlog)
expr_african_american_cpmlog <- t(expr_african_american_cpmlog)

##following DiffCoEx
#Constructing the adjacency matrix
beta1 <- 6

AdjMatC <- sign(cor(expr_caucasian_cpmlog, method = "spearman"))*(cor(expr_caucasian_cpmlog, method = "spearman"))^2
AdjMatAA <- sign(cor(expr_african_american_cpmlog, method = "spearman"))*(cor(expr_african_american_cpmlog, method = "spearman"))^2
diag(AdjMatC) <- 0
diag(AdjMatAA) <- 0
collectGarbage()
dissTOM_C_AA <- TOMdist((abs(AdjMatC - AdjMatAA)/2)^beta1/2)

collectGarbage()

geneTree_C_AA <- flashClust(as.dist(dissTOM_C_AA), method = "average")

png("NFD_c_aa_geneTree.png")
plot(geneTree_C_AA, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

##Extract modules
dynamicModsHybrid_C_AA <- cutreeDynamic(dendro = geneTree_C_AA,
                                        distM = dissTOM_C_AA, method = "hybrid", cutHeight = 0.9999,
                                        deepSplit = T, pamRespectsDendro = F, minClusterSize = 30)
print("Module extraction complete")

dynamicColorsHybrid_C_AA <- labels2colors(dynamicModsHybrid_C_AA)

mergedColour_C_AA <- mergeCloseModules(rbind(expr_caucasian_cpmlog,expr_african_american_cpmlog),dynamicColorsHybrid_C_AA,cutHeight = 0.2)$color
colorh1_C_AA <- mergedColour_C_AA

colorh1_C_AA[which(colorh1_C_AA=="midnightblue")]<-"red"
colorh1_C_AA[which(colorh1_C_AA=="lightgreen")]<-"yellow"
colorh1_C_AA[which(colorh1_C_AA=="cyan")]<-"orange"
colorh1_C_AA[which(colorh1_C_AA=="lightcyan")]<-"green"

png("DendroTree_DCM_aa_c.png")
plotDendroAndColors(geneTree_C_AA,colorh1_C_AA, "Hybrid Tree Cut", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colours cells")
dev.off()

cormat1= AdjMatC
cormat2= AdjMatAA
perm_num=100
datC1=expr_caucasian_cpmlog
datC2=expr_african_american_cpmlog

dispersionModule2Module<-function(cormat1,cormat2)
{   
        index=rownames(cormat1)
       difCor<-(cormat1[index,index]-cormat2[index,index])^2
       n<-length(index)
      (1/((n^2 -n)/2)*(sum(difCor)/2))^(.5)
}

perm_dispersionModule2Module <- function(cormat1,cormat2,perm_num,datC1,datC2)
{   
    index <- rownames(cormat1)
    difCor<-(cormat1[index,index]-cormat2[index,index])^2
    n<-length(index)
    obs=sqrt(1/((n^2 -n)/2)*(sum(difCor)/2))
    perm_dif <- numeric(perm_num) 
    for (i in 1:perm_num) {
        perm_indices <- sample(colnames(datC1),n)
        perm_cormat1 <- cor(datC1[,perm_indices],method = "spearman")
        perm_cormat2 <- cor(datC2[,perm_indices],method = "spearman")
        perm_dif[i] <- dispersionModule2Module(perm_cormat1,perm_cormat2)
    }
    p_value <- mean(perm_dif > obs)
    return(p_value)
}

modules <- setdiff(unique(colorh1_C_AA),"grey")
num_modules <- length(modules)

dispersionMatrix <- matrix(nrow=num_modules,ncol=num_modules)
nullDistrib <- list()

i <- 0

for (c1 in modules) {
    i <- i+1
    j<-0
    nullDistrib[[c1]] <- list()
    for (c2 in modules) {
        j <- j + 1
        dispersionMatrix[i,j] <- dispersionModule2Module(AdjMatC,AdjMatAA)
        nullDistrib[[c1]][[c2]] <- perm_dispersionModule2Module(AdjMatC,AdjMatAA,perm_num,expr_caucasian_cpmlog,expr_african_american_cpmlog)
    }
}

permutationSummary <- matrix(nrow=num_modules,ncol=num_modules)
colnames(permutationSummary) <- modules
rownames(permutationSummary) <- modules

for (i in 1:num_modules) {
    for (j in 1:num_modules) {
        c1 <- modules[i]
        c2 <- modules[j]
        permutationSummary[i,j] <- length(which(nullDistrib[[c1]][[c2]] >= dispersionMatrix[i,j]))
    }
}

##interpret as p-values
print("proportional p-vals")
permSum <- permutationSummary/perm_num
permSum

plotMatrix<-function(mat)
{
  mat[which(row(mat)>col(mat))]<-1001  image(mat,col=c(gray.colors(4),"white"),breaks=c(0,0.1,50,100,1000,1001),xaxt='n',yaxt='n',xlim=c(-0.2,1.2),ylim=c(-0.2,1.2),bty='n',asp=1)
  text(0:(nrow(mat)-1)/(nrow(mat)-1),1.1,rownames(mat),cex=1,col=rownames(mat))
  text(-0.15,0:(ncol(mat)-1)/(ncol(mat)-1),colnames(mat),cex=1,col=colnames(mat))
  text(apply(matrix(0:(nrow(mat)-1)/(nrow(mat)-1)),1,rep,ncol(mat)),rep(0:(ncol(mat)-1)/(ncol(mat)-1),nrow(mat)),as.numeric(t(mat)),col="white",cex=1.5)
}

png("NFD_aa_c_PermSumm.png")
plotMatrix(permutationSummary)
dev.off()

print(permutationSummary)
