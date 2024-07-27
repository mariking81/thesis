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

plot_NFD_DCM_Heatmap<-function(colorh1_NFD_DCM,AdjMatNFD,AdjMatDCM, expr_nfd_cpmlog, expr_dcm_cpmlog,ordering=NULL,file="DifferentialPlot.png")
{
  if (is.null(ordering))
  {
    h<-hclust(as.dist(1-abs(cor(getEigenGeneValues(expr_nfd_cpmlog[,which(colorh1_NFD_DCM!="grey")],colorh1_NFD_DCM[which(colorh1_NFD_DCM!="grey")],
                                                   rbind(expr_nfd_cpmlog,expr_dcm_cpmlog)[,which(colorh1_NFD_DCM!="grey")])))))
    for (c in h$label[h$order])
    {
      ordering<-c(ordering,which(colorh1_NFD_DCM ==c))
    }
  }
  mat_tmp<-(AdjMatNFD[ordering,ordering])
  mat_tmp[which(row(mat_tmp)>col(mat_tmp))]<-(AdjMatDCM[ordering,ordering][which(row(mat_tmp)>col(mat_tmp))])
  diag(mat_tmp)<-0
  mat_tmp<-sign(mat_tmp)*abs(mat_tmp)^(1/2)
  png(file=file,height=1000,width=1000)
  image(mat_tmp,col=rev(brewer.pal(11,"RdYlBu")),axes=F,asp=1,breaks=seq(-1,1,length.out=12))
  dev.off()
  unique(colorh1_NFD_DCM[ordering])
}

plotExprChange<-function(datC1,datC2, colorh1_NFD_DCM,ordering=NULL)
{
  if (is.null(ordering))
  {
    h<-hclust(as.dist(1-abs(cor(getEigenGeneValues(datC1[,which(colorh1_NFD_DCM!="grey")],colorh1_NFD_DCM[which(colorh1_NFD_DCM!="grey")],rbind(datC1,datC2)[,which(colorh1_NFD_DCM!="grey")])))))
    for (c in h$label[h$order])
    {
      ordering<-c(ordering,which(colorh1_NFD_DCM ==c))
    }
  }
  mycolors<-colorh1_NFD_DCM[ordering]
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
    meanC1<-mean(t(datC1)[colnames(datC1)[which(colorh1_NFD_DCM == c)],])
    meanC2<-mean(t(datC2)[colnames(datC2)[which(colorh1_NFD_DCM == c)],])
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

options(pkgType = "source")

#Load dependencies
library(WGCNA)
library(edgeR)
library(flashClust)
library(RColorBrewer)

expr_nfd <- read.csv("./MSc/project/DiffCoEx/AA_NFD_matrix.csv",row.names = 1)
expr_dcm <- read.csv("./MSc/project/DiffCoEx/AA_DCM_matrix.csv", row.names = 1)

common_genes <- intersect(rownames(expr_nfd), rownames(expr_dcm))
expr_nfd <- expr_nfd[common_genes, ]
expr_dcm <- expr_dcm[common_genes, ]

expr_nfd_cpmlog <- cpm(expr_nfd,log=TRUE)
expr_dcm_cpmlog <- cpm(expr_dcm,log=TRUE)

#transposal
expr_nfd_cpmlog <- t(expr_nfd_cpmlog)
expr_dcm_cpmlog <- t(expr_dcm_cpmlog)

beta1 <- 6

AdjMatNFD <- sign(cor(expr_nfd_cpmlog, method = "spearman"))*(cor(expr_nfd_cpmlog, method = "spearman"))^2
AdjMatDCM <- sign(cor(expr_dcm_cpmlog, method = "spearman"))*(cor(expr_dcm_cpmlog, method = "spearman"))^2
diag(AdjMatNFD) <- 0
diag(AdjMatDCM) <- 0
collectGarbage()

dissTOM_NFD_DCM <- TOMdist((abs(AdjMatNFD - AdjMatDCM)/2)^beta1/2)
print("Dissimilarity calculation complete")
collectGarbage()
geneTree_NFD_DCM <- flashClust(as.dist(dissTOM_NFD_DCM), method = "average")
print("Hierarchical clustering complete")

plot(geneTree_NFD_DCM, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)

dynamicModsHybrid_NFD_DCM <- cutreeDynamic(dendro = geneTree_NFD_DCM,
                                        distM = dissTOM_NFD_DCM, method = "hybrid", cutHeight = 0.9999,
                                       deepSplit = T, pamRespectsDendro = F, minClusterSize = 30)
print("Module extraction complete")

dynamicColorsHybrid_NFD_DCM <- labels2colors(dynamicModsHybrid_NFD_DCM)
table(dynamicColorsHybrid_NFD_DCM)

#merge clusters too close together
mergedColour_NFD_DCM <- mergeCloseModules(rbind(expr_nfd_cpmlog,expr_dcm_cpmlog),dynamicColorsHybrid_NFD_DCM,cutHeight = 0.2)$color
colorh1_NFD_DCM <- mergedColour_NFD_DCM

#reassign better colors
colorh1_NFD_DCM[which(colorh1_NFD_DCM=="midnightblue")]<-"red"
colorh1_NFD_DCM[which(colorh1_NFD_DCM=="lightgreen")]<-"yellow"
colorh1_NFD_DCM[which(colorh1_NFD_DCM=="cyan")]<-"orange"
colorh1_NFD_DCM[which(colorh1_NFD_DCM=="lightcyan")]<-"green"

plotDendroAndColors(geneTree_NFD_DCM,colorh1_NFD_DCM, "Hybrid Tree Cut", dendroLabels = FALSE, hang = 0.03,
                   addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colours cells")

anno <- read.csv("annotation_matrix.csv")
#write each module to an individual file containing affymetrix probeset IDs
modules_NFD_DCM_Merged <- extractModules(colorh1_NFD_DCM, expr_nfd_cpmlog,anno, dir = "AA_NFD_DCM_modules",
                                      file_prefix = paste("Output","Specific_module",sep=''),write=T)
write.table(colorh1_NFD_DCM, file = "AA_module_assignment.txt", row.names = F, col.names = F, quote = F)
plot_NFD_DCM_Heatmap(colorh1_NFD_DCM,AdjMatNFD,AdjMatDCM,expr_nfd_cpmlog,expr_dcm_cpmlog)
plotExprChange(expr_nfd_cpmlog,expr_dcm_cpmlog,colorh1_NFD_DCM)

#permutation procedure for significance testing
dispersionModule2Module <- function(c1,c2,expr_nfd_cpmlog,expr_dcm_cpmlog,colorh1_NFD_DCM)
    {
    if (c1==c2)
        {
        difCor <- (cor(expr_nfd_cpmlog[,which(colorh1_NFD_DCM == c1)],method="spearman") - cor(expr_dcm_cpmlog[,which(colorh1_NFD_DCM == c1)],method="spearman"))^2
        n <- length(which(colorh1_NFD_DCM == c1))
        (1/((n^2 -n)/2)*(sum(difCor)/2))^(0.5)
    }
    else if (c1 !=2)
        {
        difCor<-(cor(expr_nfd_cpmlog[,which(colorh1_NFD_DCM == c1)],expr_nfd_cpmlog[,which(colorh1_NFD_DCM==c2)],method="spearman")-
              cor(expr_dcm_cpmlog[,which(colorh1_NFD_DCM == c1)],expr_dcm_cpmlog[,which(colorh1_NFD_DCM==c2)],method="spearman"))^2
     n1<-length(which(colorh1_NFD_DCM  ==c1))
     n2<-length(which(colorh1_NFD_DCM  ==c2))
     (1/((n1*n2))*(sum(difCor)))^(.5)
    }
}

permutations<-NULL
for (i in 1:1000)
{
   permutations<-rbind(permutations,sample(1:(nrow(expr_nfd_cpmlog)+nrow(expr_dcm_cpmlog)),nrow(expr_nfd_cpmlog)))
}

d<-rbind(scale(expr_nfd_cpmlog),scale(expr_dcm_cpmlog))
permutationProcedureModule2Module<-function(permutation,d,c1,c2,colorh1_NFD_DCM)
{
  d1<-d[permutation,]
  d2<-d[-permutation,]
  dispersionModule2Module(c1,c2,d1,d2,colorh1_NFD_DCM)
}

dispersionMatrix<-matrix(nrow=length(unique(colorh1_NFD_DCM))-1,ncol=length(unique(colorh1_NFD_DCM))-1)
nullDistrib<-list()
i<-j<-0
for (c1 in setdiff(unique(colorh1_NFD_DCM),"grey"))
{
  i<-i+1
  j<-0
  nullDistrib[[c1]]<-list()
  for (c2 in setdiff(unique(colorh1_NFD_DCM),"grey"))
  {
    j<-j+1
    dispersionMatrix[i,j]<-dispersionModule2Module(c1,c2,expr_nfd_cpmlog,expr_dcm_cpmlog,colorh1_NFD_DCM)
    nullDistrib[[c1]][[c2]]<-apply(permutations,1,permutationProcedureModule2Module,d,c2,c1,colorh1_NFD_DCM)
  }
}

permutationSummary<-matrix(nrow=8,ncol=8)
colnames(permutationSummary)<-setdiff(unique(colorh1_NFD_DCM),"grey")
rownames(permutationSummary)<-setdiff(unique(colorh1_NFD_DCM),"grey")
for (i in 1:8) { for (j in 1:8) {permutationSummary[i,j]<-length(which(nullDistrib[[i]][[j]] >= dispersionMatrix[i,j]))}}

plotMatrix(permutationSummary)
