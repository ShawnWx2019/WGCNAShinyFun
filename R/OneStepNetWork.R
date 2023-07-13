# oneStepNetwork  -----------------------------------------------------------
#' @title oneStepNetwork
#' @description one step get scale-free network via WGCNA.
#' @author Shawn Wang <http://www.shawnlearnbioinfo.top>.
#' \email{shawnwang2016@@126.com}
#' @param exp.mat 'character' or 'data.frame'; expression matix.
#' @param lazy_modle "logical'; run oneStepWGCNA in lazy_model or not. Default is FALSE, Enabling lazy mode will result in skipping interactive steps, primarily the removal of outliers and the selection of the optimal power value. This can lead to significant differences in the final results. However, other parameters can be modified by changing the default values. If not set, the process will run according to the default parameters.
#' @param save_fig 'logical'; how many cores you want to use. Default is TRUE
#' @param F.path 'character'; output file path.
#' @param RcCutoff 'numeric'; Cutoff of rawdata, if readcount or FPKM less than this "RcCutoff" in x samples, it should be treat as the noise.
#' @param samplePerc 'numeric'; If readcount or FPKM less than RcCutoff "insamplePerc"  samples, it should be treat as the noise.
#' @param datatype 'character'; Data type of your input rawdata. "count", "expected count". "normalized count", "peak area" or "protein abundance"
#' @param method 'character'; Conversion method of input raw data. "vst" for count and expected count, "raw" and "logarithm" for other datatype.
#' @param GeneNumCut 'numeric'; How many transcripts you want to retain for WGCNA.
#' @param cutmethod 'character'; The datExpr will ranked by MAD or Var. which method you want.
#' @param rscut 'numeric'; cutoff of R2.
#' @param minModuleSize 'numeric'; minimal size of module.
#' @param mergeCutHeight 'numeric'; Dynamic cutting height, clusters less than this value will be merged
#' @param maxBlocksize 'numeric'; number of select genes, default is 9999999 means all input genes were select for next steps.
#' @return Scale-free network result.
#' @importFrom magrittr %>%
#' @importFrom plotly ggplotly
#' @importFrom dplyr select
#' @importFrom ape write.tree
#' @importFrom ggdendro ggdendrogram
#' @importFrom WGCNA plotDendroAndColors
#' @importFrom crayon green bold italic red yellow
#' @export
oneStepNetwork <- function(
    exp.mat,
    lazy_modle = FALSE,save_fig = TRUE,F.path = "./",
    datatype = 'peak area',RcCutoff = 1,samplePerc = 0.3,method = 'raw',
    GeneNumCut = 0, cutmethod = "MAD",
    rscut = 0.8,
    minModuleSize = 30,mergeCutHeight = 0.25,maxBlocksize = 999999
) {
  ##> setting
  if(isTRUE(save_fig)) {
    dir.create(F.path,showWarnings = F,recursive = T)
    message("Your result will export at: ",F.path)
  }
  msg_yes = green$bold$italic;
  msg_no = red$bold$italic;
  msg_warning = yellow$bold$italic;
  if(isTRUE(lazy_modle)){
    Sys.sleep(7)
    message(msg_warning("The WGCNA is running under lazy_modle. \nNote that under lazy_modle, the entire WGCNA process will run with default parameters. \nThis means there is a chance of overlooking issues such as sample heterogeneity and the appropriateness of the power value."))
  }
  message(msg_yes("==========================================\n"))
  message(msg_yes("Step1. Gene or Protein Expression / Metabolites Accumulation Matrix Filtring \n"))
  message(msg_yes("==========================================\n"))
  Step1 <- getdatExpr(rawdata = exp.mat,RcCutoff = RcCutoff,samplePerc = samplePerc,datatype = datatype,method = method)
  Step2 <- getdatExpr2(datExpr = Step1,GeneNumCut = GeneNumCut,cutmethod = cutmethod)
  Step3 <- getsampleTree(datExpr = Step2)
  if(isTRUE(lazy_modle)) {
    Step4 = Step2
    outlier <- NA
  } else {
    ggdendrogram(Step3$sampleTree) %>% ggplotly() %>% print()
    cat("Enter Sample_ID as a vector (e.g., c('x','b','a','d')):\n")
    outlier <- scan(what = character(), n = 0, quiet = TRUE)
    if(length(outlier) > 0) {
      Step1 <- Step1 %>% select(-all_of(outlier))
      Step2 <- getdatExpr2(datExpr = Step1,GeneNumCut = GeneNumCut,cutmethod = cutmethod)
      msg_yes(paste(outlier, "removed from raw dataset!\n"))
      Step3 <- getsampleTree(datExpr = Step2)
      ggdendrogram(Step3$sampleTree) %>%
        ggplotly()  %>%
        print()
      Step4 = Step2
    } else {
      Step4 <- Step2
      message(msg_warning("All samples reserved for next steps! Be careful of outliers!"))
    }
  }
  if(isTRUE(save_fig)) {
    pdf(file = paste0(F.path,"SampleCluster.pdf"),width = 30,height = 6)
    ggdendrogram(Step3$sampleTree) %>% print()
    dev.off()
    write.tree(phy = Step3$tree,file = paste0(F.path,"SampleCluster.nwk"))
  }
  message(msg_yes("==========================================\n"))
  message(msg_yes("Step2. Power Value Selection and Scale-Free Network Verification  \n"))
  message(msg_yes("==========================================\n"))
  Step5 <- getpower(datExpr = Step4,rscut = rscut)
  if (isTRUE(lazy_modle)) {
    power_rec <- Step5$power
  } else {
    Step5$plot %>% print
    message(msg_yes(paste0("Default SFT cutoff: ",rscut)))
    message(msg_yes(paste0("The recommended power: ",Step5$power)))
    power_rec <- readline("Choose power: ") %>% as.numeric()
    if(power_rec =="") {
      power_rec = Step5$power
    }
  }

  if(isTRUE(save_fig)) {
    pdf(file = paste0(F.path,"PowerSelection.pdf"),width = 11,height = 6)
    Step5$plot %>% print()
    dev.off()
  }

  Step6 = powertest(power.test = power_rec,datExpr = Step4,nGenes = ncol(Step4))
  Step6 %>% print()
  if(isTRUE(save_fig)) {
    pdf(file = paste0(F.path,"ScaleFreeCheck.pdf"),width = 11,height = 6)
    Step6 %>% print()
    dev.off()
  }
  message(msg_yes("==========================================\n"))
  message(msg_yes("Step3. Construction network  \n"))
  message(msg_yes("==========================================\n"))
  if(maxBlocksize == 999999){
    maxBlocksize = ncol(Step4)
  }
  Step7 = getnetwork(datExpr = Step4,power = power_rec,minModuleSize = minModuleSize,mergeCutHeight = mergeCutHeight,maxBlocksize = maxBlocksize )
  plotDendroAndColors(Step7$net$dendrograms[[1]], Step7$moduleColors[Step7$net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  if(isTRUE(save_fig)) {
    pdf(file = paste0(F.path,"ClusterDendrogram.pdf"),width = 10,height = 6)
    plotDendroAndColors(Step7$net$dendrograms[[1]], Step7$moduleColors[Step7$net$blockGenes[[1]]],
                        "Module colors",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05)
    dev.off()
  }
  out = list(
    expmat = Step4,
    power = power_rec,
    net = Step7,
    outlier = outlier
  )
  return(out)
}
