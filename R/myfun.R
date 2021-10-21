## Functions of WGCNA shiny App
suppressMessages(library(shinyjs))
suppressMessages(library(dashboardthemes))
suppressMessages(library(shinydashboard))
suppressMessages(library(DT))
suppressMessages(library(shiny))
suppressMessages(library(DESeq2))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(WGCNA))
suppressMessages(library(stringr))
suppressMessages(library(ape))
suppressMessages(library(reshape2))
suppressMessages(library(edgeR))
suppressMessages(library(shinythemes))
suppressMessages(library(ggplotify))
suppressMessages(library(ggprism))
suppressMessages(library(patchwork))
suppressMessages(library(tidyverse))
suppressMessages(library(shinyjqui))
#' Data cleaning of expression matrix.
#' Follow the rules of https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html.
#' @param rawdata Input expression data, such as FPKM, CPM, TPM or read count.
#' @param RcCutoff Cutoff of rawdata, if readcount or FPKM less than this "RcCutoff" in x samples, it should be treat as the noise.
#' @param samplePerc If readcount or FPKM less than RcCutoff "insamplePerc"  samples, it should be treat as the noise.
#' @param datatype Data type of your input rawdata. readcount or FPKM. CPM,TPM,RPM and RPKM = "FPKM"
#' @param method Conversion method of input raw data. readcount: "varianceStabilizingTransformation" or "lgcpm", FPKM: "rawFPKM" or "lgFPKM"
#' @return dx a matrix like datExpr
#' @references https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html
#' @export
#' @import DESeq2
#' @import edgeR
#' @import dplyr
#' @example
#' @author Shawn Wang <url\{http://www.shawnlearnbioinfo.top}>

getdatExpr = function(rawdata,RcCutoff,samplePerc,datatype,method){
  ## Avoid numeric id
  rawdata <- data.frame(row.names = as.character(rawdata[,1]),
                        rawdata[,-1])
  countvarTran = function(rawcount,RcCutoff,samplePerc) {
    samnum <- ncol(rawcount)
    casenum = ceiling(samnum/2)
    controlnum = samnum - casenum
    condition <- factor(c(rep("case",casenum),rep("control",controlnum)),levels = c("case","control"))
    colData <- data.frame(row.names = colnames(rawcount), condition)
    ## remove background noise
    x <- rawcount[apply(rawcount,1,function(x) sum(x > RcCutoff) > (samplePerc*ncol(rawcount))),]

    ## readcount standardization by DESeq2
    dds <- DESeqDataSetFromMatrix(x, colData, design = ~ condition)
    dds <- DESeq(dds)
    vsd <- assay(varianceStabilizingTransformation(dds))
    dx = data.frame(vsd)
    return(dx)
  }
  ## func2 count 2 cpm
  countCPM = function(rawcount, RcCutoff,samplePerc) {
    x <- rawcount[apply(rawcount,1,function(x) sum(x > RcCutoff) > (samplePerc*ncol(rawcount))),]
    dx  <- log10(edgeR::cpm(x)+1)
    return(dx)
  }
  ## raw fpkm filter
  fpkmfilter = function(rawcount, RcCutoff,samplePerc) {
    x <- rawcount[apply(rawcount,1,function(x) sum(x > RcCutoff) > (samplePerc*ncol(rawcount))),]
    dx  <- x
    return(dx)
  }

  ## log fpkm filter
  lgfpkmfilter = function(rawcount, RcCutoff,samplePerc) {
    x <- rawcount[apply(rawcount,1,function(x) sum(x > RcCutoff) > (samplePerc*ncol(rawcount))),]
    dx  <- log10(x+1)
    return(dx)
  }
  ## dx final
  if (
    datatype == "count" & method == "varianceStabilizingTransformation"
  ) {
    dx = countvarTran(rawcount = rawdata,RcCutoff = RcCutoff, samplePerc = samplePerc)
  } else if (
    datatype == "count" & method == "lgcpm"
  ) {
    dx = countCPM(rawcount = rawdata,RcCutoff = RcCutoff, samplePerc = samplePerc)
  } else if (
    datatype == "FPKM" & method == "rawFPKM"
  ) {
    dx = fpkmfilter(rawcount = rawdata,RcCutoff = RcCutoff, samplePerc = samplePerc)
  } else if (
    datatype == "FPKM" & method == "lgFPKM"
  ) {
    dx = lgfpkmfilter(rawcount = rawdata,RcCutoff = RcCutoff, samplePerc = samplePerc)
  }
  return(dx)
}

#' Data cleaning of expression matrix step2.
#' Follow the rules of https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html.
#' @param datExpr Expression data matrix generate from 1st step.
#' @param GeneNumCut How many transcripts you want to retain for WGCNA.
#' @param cutmethod The datExpr will ranked by MAD or Var. which method you want.
#' @return dx a matrix like datExpr
#' @references https://www.jianshu.com/p/e9cc3f43441d
#' @export
#' @import WGCNA
#' @import dplyr
#' @example
#' @author Shawn Wang <url\{http://www.shawnlearnbioinfo.top}>

getdatExpr2 = function(datExpr,GeneNumCut,cutmethod){
  datExpr <- datExpr
  type = "unsigned"
  corType = "pearson"
  corFnc = cor
  maxPOutliers = ifelse(corType=="pearson",1,0.05)
  robustY = ifelse(corType=="pearson",T,F)
  m.mad <- apply(datExpr,1,mad)
  m.var <- apply(datExpr,1,var)
  if(GeneNumCut == 0){
    datExprVar = datExpr
  } else if(cutmethod == "MAD"){
    datExprVar <- datExpr[which(m.mad >
                                  max(quantile(m.mad, probs=seq(0, 1, GeneNumCut))[2],0.01)),]
  } else if(cutmethod == "Var"){
    datExprVar <- datExpr[which(m.var >
                                  max(quantile(m.var, probs=seq(0, 1, GeneNumCut))[2],0.01)),]
  } else {
    print("Error: wrong filter method")
  }
  dim(datExprVar)
  datExpr <- as.data.frame(t(datExprVar))
  gsg = goodSamplesGenes(datExpr, verbose = 3)
  if (!gsg$allOK){
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
      printFlush(paste("Removing genes:",
                       paste(names(datExpr)[!gsg$goodGenes], collapse = ",")));
    if (sum(!gsg$goodSamples)>0)
      printFlush(paste("Removing samples:",
                       paste(rownames(datExpr)[!gsg$goodSamples], collapse = ",")));
    # Remove the offending genes and samples from the data:
    datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
  }
  return(datExpr)
}

#' Plot the sample tree.
#' The function of the number of sample clusters is to remove outliers.
#' @param datExpr Expression data matrix generate from 2nd step.
#' @param layout Layout of your sample tree cluster.
#' @return A list contains Tree file, and tree plot.
#' @references https://www.jianshu.com/p/e9cc3f43441d
#' @export
#' @import ape
#' @import tidyverse
#' @example
#' @author Shawn Wang <url\{http://www.shawnlearnbioinfo.top}>

getsampleTree = function(datExpr,layout = "circular"){
  ## sample cluster based on expression values
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  sampleTree = hclust(dist(datExpr), method = "average")
  treenew= as.phylo(sampleTree)
  # anno = data.frame(lab = treenew$tip.label,
  #                   group = factor(gsub(".$","",treenew$tip.label),levels = unique(gsub(".$","",treenew$tip.label))))
  # p =ggtree(tr = treenew,layout = layout)
  # p2 = p %<+% anno +geom_tiplab(aes(color=group))+theme(legend.position = "")
  out = list(tree = treenew,
             sampleTree = sampleTree,
             nGenes = nGenes,
             nSamples = nSamples)
  return(out)
}


#' SFT detection.
#' Calculate the appropriate power for scale-free network construction.
#' @param datExpr Expression data matrix generate from 2nd step.
#' @param rscut cutoff of R2.
#' @return A list contains power test detials, power and two plots which demonstrated the sft and mean connectivity.
#' @references https://www.jianshu.com/p/e9cc3f43441d
#' @export
#' @import WGCNA
#' @import tidyverse
#' @import ggprism
#' @import patchwork
#' @example
#' @author Shawn Wang <url\{http://www.shawnlearnbioinfo.top}>
getpower = function(datExpr,rscut){
  nSamples = nrow(datExpr)
  type = 'unsigned'
  powers = c(c(1:10), seq(from = 12, to=30, by=2))
  sft = pickSoftThreshold(datExpr, powerVector=powers, RsquaredCut = rscut,
                          networkType=type, verbose=5)
  fitIn = sft$fitIndices
  min.sft = min(-sign(sft$fitIndices[,3])*sft$fitIndices[,2])
  fitIn %>%
    mutate(r2 = -sign(slope)*SFT.R.sq) %>%
    ggplot(.,aes(x = Power,y = r2,color = "red"))+
    geom_text(aes(label = Power))+
    theme_bw()+
    ylim(min.sft,1)+
    geom_hline(aes(yintercept = rscut),color = "green",linetype = "dashed")+
    geom_hline(aes(yintercept = 0.8),color = "blue",linetype = "dashed")+
    geom_hline(aes(yintercept = 0.9),color = "red",linetype = "solid")+
    theme_prism(border = T)+
    xlab("Soft Threshold (Power)")+
    ylab(expression("Scale Free Topology Model Fit,signed R"^2))+
    ggtitle("Scale independence")+
    theme(legend.position = "")->p1
  fitIn %>%
    ggplot(.,aes(x = Power,y = mean.k.,color = "red"))+
    geom_text(aes(label = Power))+
    theme_bw()+
    theme_prism(border = T)+
    xlab("Soft Threshold (Power)")+
    ylab("Mean Connectivity")+
    ggtitle("Mean Connectivity")+
    theme(legend.position = "")->p2
  p.all = p1+p2
  power = sft$powerEstimate
  ## exper power
  if (is.na(power)){
    power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                   ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                          ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                                 ifelse(type == "unsigned", 6, 12))
                   )
    )
  }
  out = list(sft = fitIn,
             plot = p.all,
             power = power)

  return(out)
}

#' ggplot2 version of scaleFreePlot.
#' Transformed from the WGCNA::scaleFreePlot function.
#' @param connectivity See WGCNA::scaleFreePlot.
#' @param nBreaks See WGCNA::scaleFreePlot.
#' @param truncated See WGCNA::scaleFreePlot
#' @return A list contains Tree file, and tree plot.
#' @references https://rdrr.io/cran/WGCNA/man/scaleFreePlot.html
#' @export
#' @import WGCNA
#' @import tidyverse
#' @import ggprism
#' @import stringr
#' @example
#' @author Shawn Wang <url\{http://www.shawnlearnbioinfo.top}>

ggscaleFreePlot = function (connectivity, nBreaks = 10, truncated = FALSE, removeFirst = FALSE, main = "", ...) {
  k = connectivity
  discretized.k = cut(k, nBreaks)
  dk = tapply(k, discretized.k, mean)
  p.dk = as.vector(tapply(k, discretized.k, length)/length(k))
  breaks1 = seq(from = min(k), to = max(k), length = nBreaks +
                  1)
  hist1 = suppressWarnings(hist(k, breaks = breaks1, equidist = FALSE,
                                plot = FALSE, right = TRUE, ...))
  dk2 = hist1$mids
  dk = ifelse(is.na(dk), dk2, dk)
  dk = ifelse(dk == 0, dk2, dk)
  p.dk = ifelse(is.na(p.dk), 0, p.dk)
  log.dk = as.vector(log10(dk))
  if (removeFirst) {
    p.dk = p.dk[-1]
    log.dk = log.dk[-1]
  }
  log.p.dk = as.numeric(log10(p.dk + 1e-09))
  lm1 = lm(log.p.dk ~ log.dk)
  if (truncated == TRUE) {
    lm2 = lm(log.p.dk ~ log.dk + I(10^log.dk))
    OUTPUT = data.frame(scaleFreeRsquared = round(summary(lm1)$adj.r.squared,
                                                  2), slope = round(lm1$coefficients[[2]], 2), TruncatedRsquared = round(summary(lm2)$adj.r.squared,
                                                                                                                         2))
    printFlush("the red line corresponds to the truncated exponential fit")
    title = paste(main, " scale free R^2=", as.character(round(summary(lm1)$adj.r.squared,
                                                               2)), ", slope=", round(lm1$coefficients[[2]], 2),
                  ", trunc.R^2=", as.character(round(summary(lm2)$adj.r.squared,
                                                     2)))
  }
  else {
    title = paste(main, " scale R^2=", as.character(round(summary(lm1)$adj.r.squared,
                                                          2)), ", slope=", round(lm1$coefficients[[2]], 2))
    OUTPUT = data.frame(scaleFreeRsquared = round(summary(lm1)$adj.r.squared,
                                                  2), slope = round(lm1$coefficients[[2]], 2))
  }
  pdata = data.frame(x = log.dk,
                     y = log.p.dk)
  lm = if(truncated == F){
    lm1
  } else {
    lm2
  }

  p = ggplot(pdata ,aes(x = x,y = y))+
    geom_point(shape  = 1,size = 5)+
    geom_abline(intercept = as.numeric(lm$coefficients[1]),slope =as.numeric(lm$coefficients[2]) )+
    xlab("log10(k)")+
    ylab("log10(p(k))")+
    ggtitle(title)+theme_prism(border = T)

}


#' Power test.
#' Test scale-free of corresponding power
#' @param power.test Power you want test.
#' @param datExpr datExpr generate in second step.
#' @param nGenes Select gene number
#' @return Plot of scale-free test result
#' @references https://www.jianshu.com/p/25905a905086 jiashu Author: sxzt
#' @export
#' @import WGCNA
#' @import tidyverse
#' @import ggprism
#' @import patchwork
#' @example
#' @author Shawn Wang <url\{http://www.shawnlearnbioinfo.top}>

powertest = function(power.test,datExpr,nGenes){
  ## reference: https://www.jianshu.com/p/25905a905086 jiashu Author: sxzt
  softPower <- power.test
  adjacency = adjacency(datExpr = datExpr, power = softPower)
  # convert adj 2 TOM
  TOM = TOMsimilarity(adjacency)
  # disTOM
  dissTOM = 1-TOM
  hierTOM = hclust(as.dist(dissTOM),method="average")
  ADJ1_cor <- abs(WGCNA::cor( datExpr,use = "p" ))^softPower
  if(nGenes< 5000){
    k <- as.vector(apply(ADJ1_cor,2,sum,na.rm=T))
  } else {
    k <- softConnectivity(datE=datExpr,power=softPower)
  }
  k.df = data.frame(value = k)
  p1 = ggplot(k.df,mapping = aes(x = value)) +
    geom_histogram(color = "black",fill = "white")+
    theme_prism(border = T)+
    xlab("k")+
    ylab("Frequence")+
    ggtitle(paste("soft connectivity (power = ",power.test,")"))+
    theme(legend.position = "")
  p2 = ggscaleFreePlot(k,main = "Check Scale free topology\n")
  p = p1+p2
  return(p)
}

#' Scale-free network generation.
#' By WGCNA::blockwiseModules
#' @param datExpr datExpr produced by 2nd step.
#' @param power power you select after scale-free test.
#' @param minModuleSize minimal size of module.
#' @param mergeCutHeight Dynamic cutting height, clusters less than this value will be merged
#' @param maxBlocksize number of select genes
#' @return A list contains multiple parameter for next step.
#' @references https://www.jianshu.com/p/f0409a045dab
#' @export
#' @import WGCNA
#' @import tidyverse
#' @import stringr
#' @example
#' @author Shawn Wang <url\{http://www.shawnlearnbioinfo.top}>


getnetwork = function(datExpr,power,minModuleSize,mergeCutHeight,maxBlocksize){
  cor <- WGCNA::cor
  net = blockwiseModules(datExpr, power = power, maxBlockSize = maxBlocksize,
                         TOMType = "unsigned", minModuleSize = minModuleSize,
                         reassignThreshold = 0, mergeCutHeight = mergeCutHeight,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs=F, corType = 'pearson',
                         maxPOutliers=1,
                         verbose = 3)
  moduleLabels = net$colors
  moduleColors = labels2colors(moduleLabels)
  MEs = net$MEs
  MEs_col = MEs
  colnames(MEs_col) = paste0("ME", labels2colors(
    as.numeric(str_replace_all(colnames(MEs),"ME",""))))
  MEs_col = orderMEs(MEs_col)
  x = table(moduleColors)
  Gene2module <- data.frame(GID = colnames(datExpr),
                            Module = moduleColors)
  out = list(net = net,
             moduleLabels =moduleLabels,
             moduleColors = moduleColors,
             MEs_col = MEs_col,
             MEs = MEs,
             Gene2module = Gene2module
  )
  return(out)
}

#' Module-trait prepare.
#' Generate the variables needed by the module-trait relationships.
#' @param phenotype Phenotype matrix, with sample in row and traits in column.
#' @param MEs_col module eigengene.
#' @param moduleColors module to color.
#' @param datExpr datExpr produced by 2nd step.
#' @return A list contains multiple parameter for next step.
#' @references https://www.jianshu.com/p/f0409a045dab
#' @export
#' @import WGCNA
#' @import tidyverse
#' @example
#' @author Shawn Wang <url\{http://www.shawnlearnbioinfo.top}>
getMt = function(phenotype,MEs_col,nSamples,moduleColors,datExpr){
  traitData <- phenotype
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  modTraitCor = cor(MEs, traitData , use = "p")
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
  ## bugfix make mod cor and mod p as df
  df_modTraitCor = as.data.frame(modTraitCor)
  df_modTraitP = as.data.frame(modTraitP)
  mod_cor_R2 = df_modTraitCor %>%
    mutate(across(where(is.numeric), signif, 2)) %>%
    mutate_if(is.numeric,as.character)
  mod_P_R2 = df_modTraitP %>%
    mutate(across(where(is.numeric), signif, 2)) %>%
    mutate_if(is.numeric,as.character)
  ## for heatmap cluster
  tmp3 <- mod_cor_R2
  for (i in 1:nrow(mod_cor_R2)) {
    for (j in 1:ncol(mod_cor_R2)) {
      tmp3[i,j] = paste0(mod_cor_R2[i,j],"\n(",mod_P_R2[i,j],")")
    }
  }

  out = list(modTraitCor = modTraitCor,
             modTraitP = modTraitP,
             textMatrix = tmp3)
  return(out)

}

#' Get kME.
#' kME values export.
#' @param datExpr datExpr produced by 2nd step.
#' @param MEs_col module eigengene.
#' @param moduleColors module to color.
#' @return A list contains multiple parameter for next step.
#' @references https://www.jianshu.com/p/f0409a045dab
#' @export
#' @import WGCNA
#' @import tidyverse
#' @example
#' @author Shawn Wang <url\{http://www.shawnlearnbioinfo.top}>

getKME = function(datExpr,moduleColors,MEs_col){
  connet=abs(cor(datExpr,use="p"))^6
  Alldegrees1=intramodularConnectivity(connet, moduleColors)
  datKME=signedKME(datExpr, MEs_col, outputColumnName="MM.")
  return(datKME)
}

#' Get MM.
#' module membership.
#' @param datExpr datExpr produced by 2nd step.
#' @param MEs_col module eigengene.
#' @param nSamples module to color.
#' @param corType correlation types. such as "pearson"
#' @return A list contains gene module membership and Pvalue.
#' @references https://www.jianshu.com/p/f0409a045dab
#' @export
#' @import WGCNA
#' @import tidyverse
#' @example
#' @author Shawn Wang <url\{http://www.shawnlearnbioinfo.top}>


getMM = function(datExpr,MEs_col,nSamples,corType){
  if (corType=="pearson") {
    geneModuleMembership = as.data.frame(cor(datExpr, MEs_col, use = "p"))
    MMPvalue = as.data.frame(corPvalueStudent(
      as.matrix(geneModuleMembership), nSamples))
  } else {
    geneModuleMembershipA = bicorAndPvalue(datExpr, MEs_col, robustY=robustY)
    geneModuleMembership = geneModuleMembershipA$bicor
    MMPvalue   = geneModuleMembershipA$p
  }
  out = list(MM= geneModuleMembership,
             MMP = MMPvalue)
}

#' GS vs MM.
#' module membership.
#' @param datExpr datExpr produced by 2nd step.
#' @param module which module you choose to observation.
#' @param pheno which traits you choose to observation.
#' @param traitData trait matrix.
#' @param moduleColors param generated by step getnetwork.
#' @param geneModuleMembership param generated by getMM
#' @param MEs param generated by getnetwork
#' @param nSamples how many samples.
#' @return A list contains gene module membership and Pvalue.
#' @references https://www.jianshu.com/p/f0409a045dab
#' @export
#' @import WGCNA
#' @import tidyverse
#' @import ggprism
#' @import stringr
#' @example
#' @author Shawn Wang <url\{http://www.shawnlearnbioinfo.top}>

getverboseplot = function(datExpr,module,pheno,traitData,moduleColors,geneModuleMembership,MEs,nSamples){
  modNames = substring(names(MEs), 3)
  tmp1 = traitData %>%
    select(pheno);
  #  print(dim(tmp1))
  #  print(dim(datExpr))
  geneTraitSignificance = as.data.frame(cor(datExpr, tmp1, use = "p"));
  #  print(dim(geneTraitSignificance))
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  names(geneTraitSignificance) = paste("GS.", names(tmp1), sep="");
  #  print(dim(geneTraitSignificance))
  names(GSPvalue) = paste("p.GS.", names(tmp1), sep="")
  module = module
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  plt.df = data.frame(x = abs(geneModuleMembership[moduleGenes, column]),
                      y = abs(geneTraitSignificance[moduleGenes, 1]))
  #formula <- y ~ x
  corExpr = parse(text = paste("cor", "(plt.df$x, plt.df$y ", prepComma("use = 'p'"),
                               ")"))
  cor = signif(eval(corExpr), 2)
  if (is.finite(cor)) {
    if (abs(cor) < 0.00001){ cor = 0 }
  }
  corp = signif(corPvalueStudent(cor, sum(is.finite(plt.df$x) & is.finite(plt.df$y))), 2)
  if (is.finite(corp) && corp < 10^(-200)) {
    corp = "<1e-200"
  } else {corp = paste("=", corp, sep = "")}

  ggplot(plt.df,aes(x = x, y = y))+
    geom_point(shape = 21,color = "black",fill = module,size = 3)+
    geom_smooth(method = 'lm', formula = y ~ x)+
    xlab(paste("Module Membership in", module, "module"))+
    ylab(paste("Gene significance for ",pheno))+
    ggtitle(paste("Module membership vs. gene significance\n","cor=",cor," ,p",corp,sep = ""))+
    theme_bw()+
    theme_prism(border = T)
}

#' Module heatmap.
#' Gene in module heatmap.
#' @param datExpr datExpr produced by 2nd step.
#' @param which.module which module you choose to observation.
#' @param moduleColors param generated by step getnetwork.
#' @param MEs param generated by getnetwork
#' @return A heatmap of selected module.
#' @export
#' @import WGCNA
#' @import tidyverse
#' @import ggprism
#' @import stringr
#' @import patchwork
#' @example
#' @author Shawn Wang <url\{http://www.shawnlearnbioinfo.top}>
moduleheatmap = function(datExpr,MEs,which.module,moduleColors){
  ME=MEs[,paste("ME",which.module,sep = "")]
  ## heatmap data
  h <- t(scale(datExpr[,moduleColors == which.module]))
  x = reshape2::melt(data = h,id.var = rownames(h))
  ## barplot data
  b <- data.frame(Sample = factor(rownames(datExpr),levels = rownames(datExpr)) ,
                  Eigengene = ME)
  p1 = ggplot(x, aes(Var2, Var1)) +
    geom_tile(aes(fill = value),colour = alpha("white",0)) +
    scale_fill_gradient2(low = "green",mid = "black",high = "red")+
    theme_classic()+
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 90,face = "bold",vjust = 0.5,size = 15),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.x = element_blank())
  ## barplot
  p2 = ggplot(b,mapping  = aes(x = Sample,y = Eigengene)) + geom_bar(stat='identity',fill = which.module)+
    theme_classic()+
    theme(axis.title = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  p3 = p1/p2
  return(p3)
}


#' select hubgene
#' select hubgenes with GS and MM
#' @param datExpr datExpr produced by 2nd step.
#' @param mdl which module you choose to observation.
#' @param power power you selected.
#' @param trt trait you selected.
#' @param KME param generated by.
#' @param GS.cut cutoff Gene significance.
#' @param kME.cut cutoff kME.
#' @param datTrait trait matrix.
#' @return A list of hubgenes.
#' @export
#' @import WGCNA
#' @import tidyverse
#' @import stringr
#' @example
#' @author Shawn Wang <url\{http://www.shawnlearnbioinfo.top}>
hubgenes = function(datExpr,mdl,power,trt,KME,GS.cut,kME.cut,datTrait){
  KME.clean = KME
  colnames(KME.clean) = gsub("MM.","",colnames(KME))
  hub1 = chooseTopHubInEachModule(datExpr = datExpr,colorh = mdl,power = power,type = "unsigned")
  ## GS calculate
  which.module = mdl
  which.trait = datTrait %>%
    select(trt);
  GS1 = as.numeric(cor(which.trait,datExpr, use="p"))
  GS.df = data.frame(GeneID = colnames(datExpr),
                     GS = GS1)
  GS2=abs(GS1)
  ## MM
  kME.mod = KME.clean %>%
    select(mdl)
  kME.mod.df = data.frame(GeneID = rownames(kME.mod),
                          kME = kME.mod[,1])
  ## judge
  hub2 = abs(kME.mod) > kME.cut & GS2 > GS.cut
  hub2.num = table(hub2)
  print(paste0("MM and GS method selected out " ,as.numeric(hub2.num[2])," hub genes."))
  colnames(hub2) = "judge"
  hub3 <- hub2 %>%as.data.frame() %>%  filter(judge == 'TRUE') %>%
    rownames_to_column(var = "GeneID") %>%
    left_join(.,GS.df,by = "GeneID") %>%
    left_join(.,kME.mod.df,by = "GeneID")
  hub1 = data.frame(GeneID = as.character(hub1)) %>%
    left_join(.,GS.df,by = "GeneID") %>%
    left_join(.,kME.mod.df,by = "GeneID")


  out = list(hub1 = hub1,
             hub3 = hub3)
}


#' GS vs MM.all
#' module membership.
#' @param datExpr datExpr produced by 2nd step.
#' @param which.trait which trait you choose to observation.
#' @param traitData trait matrix.
#' @param moduleColors param generated by step getnetwork.
#' @param geneModuleMembership param generated by getMM
#' @param MEs param generated by getnetwork
#' @param nSamples how many samples.
#' @return A list contains gene module membership and Pvalue.
#' @references https://www.jianshu.com/p/f0409a045dab
#' @export
#' @import WGCNA
#' @import tidyverse
#' @import ggprism
#' @import stringr
#' @example
#' @author Shawn Wang <url\{http://www.shawnlearnbioinfo.top}>
MMvsGSall = function(
  which.trait = exp.ds$st,
  traitData = exp.ds$phen,
  datExpr = exp.ds$table2,
  moduleColors = exp.ds$moduleColors,
  geneModuleMembership = exp.ds$MM,
  MEs = exp.ds$MEs_col,
  nSamples = exp.ds$nSample
){
  connet=abs(cor(datExpr,use="p"))^6
  Alldegrees1=intramodularConnectivity(connet, moduleColors)
  s.trait = traitData %>%
    select(which.trait);
  names(s.trait) = which.trait
  GS1 = as.numeric(cor(s.trait,datExpr, use="p"))
  GeneSignificance=abs(GS1)
  colorlevels=unique(moduleColors)
  # par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2)))
  # par(mar = c(4,5,3,1))
  # for (i in c(1:length(colorlevels)))
  # {
  #   whichmodule=colorlevels[[i]];
  #   restrict1 = (moduleColors==whichmodule);
  #   verboseScatterplot(Alldegrees1$kWithin[restrict1],
  #                      GeneSignificance[restrict1], col=moduleColors[restrict1],
  #                      main=whichmodule,
  #                      xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
  # }
  p = list()
  for (i in c(1:length(colorlevels))) {
    p[[i]] = getverboseplot.all(
      datExpr = datExpr,
      module = colorlevels[i],
      pheno = which.trait,
      traitData = traitData,
      moduleColors = moduleColors,
      geneModuleMembership = geneModuleMembership,
      MEs = MEs,
      nSamples = nSamples)
  }
  p.out = ggpubr::ggarrange(plotlist = p,ncol = as.integer(0.5+length(colorlevels)/2),nrow = 2)
  return(p.out)
}

#' getverboseplot.all
#' module membership.
#' @param datExpr datExpr produced by 2nd step.
#' @param module which module you choose to observation.
#' @param pheno which traits you choose to observation.
#' @param traitData trait matrix.
#' @param moduleColors param generated by step getnetwork.
#' @param geneModuleMembership param generated by getMM
#' @param MEs param generated by getnetwork
#' @param nSamples how many samples.
#' @return A list contains gene module membership and Pvalue.
#' @references https://www.jianshu.com/p/f0409a045dab
#' @export
#' @import WGCNA
#' @import tidyverse
#' @import ggprism
#' @import stringr
#' @example
#' @author Shawn Wang <url\{http://www.shawnlearnbioinfo.top}>
getverboseplot.all = function(datExpr,module,pheno,traitData,moduleColors,geneModuleMembership,MEs,nSamples){
  modNames = substring(names(MEs), 3)
  tmp1 = traitData %>%
    select(pheno);
  #  print(dim(tmp1))
  #  print(dim(datExpr))
  geneTraitSignificance = as.data.frame(cor(datExpr, tmp1, use = "p"));
  #  print(dim(geneTraitSignificance))
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  names(geneTraitSignificance) = paste("GS.", names(tmp1), sep="");
  #  print(dim(geneTraitSignificance))
  names(GSPvalue) = paste("p.GS.", names(tmp1), sep="")
  module = module
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  plt.df = data.frame(x = abs(geneModuleMembership[moduleGenes, column]),
                      y = abs(geneTraitSignificance[moduleGenes, 1]))
  #formula <- y ~ x
  corExpr = parse(text = paste("cor", "(plt.df$x, plt.df$y ", prepComma("use = 'p'"),
                               ")"))
  cor = signif(eval(corExpr), 2)
  if (is.finite(cor)) {
    if (abs(cor) < 0.00001){ cor = 0 }
  }
  corp = signif(corPvalueStudent(cor, sum(is.finite(plt.df$x) & is.finite(plt.df$y))), 2)
  if (is.finite(corp) && corp < 10^(-200)) {
    corp = "<1e-200"
  } else {corp = paste("=", corp, sep = "")}

  ggplot(plt.df,aes(x = x, y = y))+
    geom_point(shape = 21,color = "black",fill = module,size = 3)+
    geom_smooth(method = 'lm', formula = y ~ x)+
    xlab(paste("MM in", module, "module"))+
    ylab(paste("GS for ",pheno))+
    ggtitle(paste("cor=",cor," ,p",corp,sep = ""))+
    theme_bw()+
    theme_prism(border = T)


}



#' getverboseplot.all
#' module membership.
#' @param datExpr datExpr produced by 2nd step.
#' @param power power value.
#' @param module which module you choose to observation.
#' @param moduleColors param generated by step getnetwork.
#' @param threshold threshold of weight
#' @return A list contains gene module membership and Pvalue.
#' @references https://www.jianshu.com/p/f0409a045dab
#' @export
#' @import WGCNA
#' @import tidyverse
#' @import stringr
#' @example
#' @author Shawn Wang <url\{http://www.shawnlearnbioinfo.top}>
cytoscapeout = function(
  datExpr = exp.ds$table2,
  power = exp.ds$power,
  module = exp.ds$mdl,
  moduleColors = exp.ds$moduleColors,
  threshold = exp.ds$threshold
){
  ##
  TOM = TOMsimilarityFromExpr(datExpr, power = power);
  # Select module
  module = module;
  # Select module probes
  probes = colnames(datExpr) ## 我们例子里面的probe就是基因名
  inModule = (moduleColors==module);
  modProbes = probes[inModule];
  ## 也是提取指定模块的基因名
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)
  head(modTOM)
  cyt = exportNetworkToCytoscape(
    modTOM,
    #edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
    #nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
    weighted = TRUE,
    threshold = threshold,
    nodeNames = modProbes,
    nodeAttr = moduleColors[inModule]
  )
  return(cyt)
}
