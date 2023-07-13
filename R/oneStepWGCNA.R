# oneStepWGCNA  -----------------------------------------------------------
#' @title oneStepWGCNA
#' @description one step WGCNA pipeline
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
#' @param r_cutoff 'numeric'; cutoff of module-trait correlation value r, default is 0.8
#' @param weight_cutoff 'numeric'; cutoff of cytoscape weight value. default is 0.02
#' @param pheofile 'character'; phenotyp file.
#' @return All result of WGCNA
#' @importFrom magrittr %>%
#' @importFrom plotly ggplotly
#' @importFrom dplyr select filter mutate group_by summarise pull left_join inner_join rename arrange
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom writexl write_xlsx
#' @importFrom purrr map2 map map2_dfr
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom ComplexHeatmap rowAnnotation Heatmap row_order column_order draw
#' @importFrom circlize colorRamp2
#' @importFrom stringr str_remove
#' @importFrom readr write_csv
#' @importFrom crayon green bold italic red yellow
#' @export

oneStepWGCNA <- function(
    exp.mat,
    lazy_modle = TRUE, save_fig = T, F.path = "./", iterative = T, grey_cut = 10,
    datatype = 'peak area', RcCutoff = 1, samplePerc = 0.3, method = 'raw',
    GeneNumCut = 0, cutmethod = "MAD",
    rscut = 0.8,
    minModuleSize = 30, mergeCutHeight = 0.25, maxBlocksize = 999999,
    r_cutoff = 0.8,weight_cutoff = 0.02,
    pheofile
) {
  ##> 01. Parameters
  msg_yes <- green$bold$italic;
  msg_no <- red$bold$italic;
  msg_warning <- yellow$bold$italic;
  ##> 02. Set default parameters
  expmat_new <- NULL
  grey_num <- NA
  round_x <- 0
  tmp_network <- NULL
  ##> 03. Output file path
  if(isTRUE(save_fig)) {
    if(str_detect(F.path,pattern = "/$")){
      out_path = paste0(F.path,"WGCNA_out/")
    } else {
      out_path = paste0(F.path,"/WGCNA_out/")
    }
    dir.create(paste0(out_path),showWarnings = F,recursive = T)
    message("Your result will export at: ",out_path)
  }
  ##> 04. run wgcna
  ##> 4.1 iterative
  path_01 = paste0(out_path,"01.ScaleFreeNetwork/")
  dir.create(path = path_01,showWarnings = F,recursive = T)
  if(isTRUE(iterative)) {
    while (is.na(grey_num) || grey_num > grey_cut ) {
      Sys.sleep(time = 2)
      round_x <- round_x + 1
      message(msg_yes(paste0("iterative WGCNA round ", round_x)))
      if (is.null(expmat_new)) {
        expmat_iteract <- exp.mat
      } else {
        expmat_iteract <- expmat_new
      }
      ##> constraction network
      tmp_network <- suppressMessages(
        oneStepNetwork(
          exp.mat = expmat_iteract,
          lazy_modle = lazy_modle,
          save_fig = save_fig,
          F.path = paste0(path_01, "round_", round_x,"/"),
          method = method,
          datatype = datatype,
          RcCutoff = RcCutoff,
          samplePerc = samplePerc,
          GeneNumCut = GeneNumCut,
          cutmethod = cutmethod,
          rscut = rscut,
          minModuleSize = minModuleSize,
          mergeCutHeight = mergeCutHeight,
          maxBlocksize = maxBlocksize
        )
      )
      ##> check unclassified network
      if("grey" %in% (tmp_network$net$Gene2module %>% pull(Module))) {
        grey_num <- tmp_network$net$Gene2module %>%
          group_by(Module) %>%
          summarise(n = n()) %>%
          filter(Module == "grey") %>%
          pull(n)
      } else {
        grey_num  = 0
      }
      ##> filter genes
      GID_rest <- data.frame(GID = tmp_network$net$Gene2module %>%
                               filter(Module != "grey") %>%
                               pull(GID))
      ##> expression matrix
      colnames(exp.mat)[1] <- "GID"
      if(round_x == 1) {
        exp.mat2 <- exp.mat
      } else {
        exp.mat2 <- expmat_new
      }
      ##> generate expmat matrix
      expmat_new <- GID_rest %>% left_join(exp.mat2, by = "GID")
      if(!is.na(tmp_network$outlier[1])) {
        expmat_new <- expmat_new %>% select(-all_of(tmp_network$outlier))
      }
      message(msg_warning(paste0("Grey module num:", grey_num)))
    }
  } else {
    ##> 4.2 tmp_network
    tmp_network <- suppressMessages(
      oneStepNetwork(
        exp.mat = exp.mat,
        lazy_modle = lazy_modle,
        save_fig = save_fig,
        F.path = path_01,
        method = method,
        datatype = datatype,
        RcCutoff = RcCutoff,
        samplePerc = samplePerc,
        GeneNumCut = GeneNumCut,
        cutmethod = cutmethod,
        rscut = rscut,
        minModuleSize = minModuleSize,
        mergeCutHeight = mergeCutHeight,
        maxBlocksize = maxBlocksize
      )
    )
  }




  ##> 05. Run module trait
  #return(tmp_network)
  message(msg_yes("==========================================\n"))
  message(msg_yes("Module trait  \n"))
  message(msg_yes("==========================================\n"))
  path_02 = paste0(out_path,"02.ModuleTrait/")
  dir.create(path = path_02,showWarnings = F,recursive = T)
  ##> 5.1 parameters
  MEs_col = tmp_network$net$MEs_col
  moduleColors = tmp_network$net$moduleColors
  Gene2Module = tmp_network$net$Gene2module %>% arrange(Module)
  datExpr = tmp_network$expmat
  nSamples = nrow(datExpr)
  nGenes = ncol(datExpr)
  corType = "pearson"
  power = tmp_network$power
  ## export tbls
  writexl::write_xlsx(x = Gene2Module,path = paste0(path_02,"Gene2Module.xlsx"))
  ##> 5.2 generate phenotype file
  if (ncol(pheofile) == 2) {
    phenotype <-
      pheofile %>%
      setNames(c("ID","Class")) %>%
      mutate(num = 1) %>%
      pivot_wider(names_from = Class,values_from = num,values_fill = 0) %>%
      column_to_rownames("ID")
  } else {
    phenotype <-
      pheofile %>%
      rename("ID" = colnames(.)[1]) %>%
      column_to_rownames("ID")
  }
  ##> 5.3 fix phenotype file
  phenotype = phenotype[match(rownames(tmp_network$expmat),rownames(phenotype)),]
  ##> 6.0 Generate trait-module heatmap
  res_kme <- getKME(datExpr = datExpr,moduleColors = moduleColors,
                    MEs_col = MEs_col)
  res_MM <- getMM(datExpr = datExpr,MEs_col = MEs_col,nSamples = nSamples,corType = corType)
  res_GS <- getMt(phenotype = phenotype,nSamples = nSamples,moduleColors = moduleColors,datExpr = datExpr)
  modTraitCor = res_GS$modTraitCor
  mod_color = gsub(pattern = "^..",replacement = "",rownames(modTraitCor))
  mod_color_anno = setNames(mod_color,rownames(modTraitCor))
  writexl::write_xlsx(x = res_kme %>% rownames_to_column("ID"),path = paste0(path_02,"kME.xlsx"))
  writexl::write_xlsx(x = list(
    r = res_GS$modTraitCor %>% as.data.frame() %>% rownames_to_column("ID"),
    p = res_GS$modTraitP %>% as.data.frame() %>% rownames_to_column("ID")
  ),path = paste0(path_02,"Module-traitCorr.xlsx"))
  Left_anno = rowAnnotation(
    Module = rownames(modTraitCor),
    col = list(
      Module = mod_color_anno
    ),
    show_legend = F,
    show_annotation_name = F
  )
  ht =
    Heatmap(
      matrix = modTraitCor,
      cluster_rows = F, cluster_columns = F,
      left_annotation = Left_anno,
      cell_fun = function(j,i,x,y,width,height,fill) {
        grid.text(sprintf(res_GS$textMatrix[i,j]),x,y,gp = gpar(fontsize = 12))
      },
      row_names_side = "left",
      column_names_rot = 90,
      heatmap_legend_param = list(
        at = c(-1,-0.5,0,0.5, 1),
        labels = c("-1","-0.5", "0","0.5", "1"),
        title = "",
        legend_height = unit(9, "cm"),
        title_position = "lefttop-rot"
      ),
      rect_gp = gpar(col = "black", lwd = 1.2),
      column_title = "Module-trait relationships",
      column_title_gp = gpar(fontsize = 15, fontface = "bold"),
      col = colorRamp2(c(-1, 0, 1), c("blue","white","red"))
    )
  draw(ht) %>% print()
  if(isTRUE(save_fig)) {
    pdf(file = paste0(path_02,"Module_trait.pdf"),width = 10,height = 10)
    draw(ht) %>% print()%>% print()
    dev.off()
  }
  ##> 7.0 output MM vs GS
  ##> 7.1 filter sig module-trait
  path_03 <- paste0(out_path,"03.HubGene/")
  dir.create(path = path_03,showWarnings = F,recursive = T)
  res_GS.r_long <-
    res_GS$modTraitCor %>%
    as.data.frame() %>%
    rownames_to_column("MODULE") %>%
    pivot_longer(!MODULE,names_to = "Cluster",values_to = "r")
  res_GS.p_long <-
    res_GS$modTraitP %>%
    as.data.frame() %>%
    rownames_to_column("MODULE") %>%
    pivot_longer(!MODULE,names_to = "Cluster",values_to = "p")
  res_GS_sig <- left_join(res_GS.r_long,res_GS.p_long) %>%
    filter(p < 0.05 & r > r_cutoff)
  r.max <- res_GS.r_long %>% pull(r) %>% max
  if(ncol(res_GS_sig) == 0) {
    message(msg_no("No significant module-trait relationships detected under r cutoff: ",r_cutoff,"\n",
                   "The max r value is: "),r.max)
  } else {
    purrr::map2(.x = res_GS_sig %>% pull(MODULE),.y = res_GS_sig %>% pull(Cluster),.f = function(.x,.y){
      mod <- .x %>% str_remove("ME")
      pdf(file = paste0(path_03,mod,"_modHeatEig.pdf"),width = 10,height = 5)
      moduleheatmap(datExpr = datExpr,MEs = MEs_col,which.module = mod,moduleColors = moduleColors) %>% print
      dev.off()
      pdf(file = paste0(path_03,mod,"-",.y,"_verboseplot.pdf"),width = 10,height = 5)
      getverboseplot(datExpr = datExpr,module = mod,pheno = .y,MEs = MEs_col,
                     traitData = phenotype,moduleColors = moduleColors,
                     geneModuleMembership = res_MM$MM,nSamples = nSamples) %>% print()
      dev.off()
    })
    KeyModule_Trait <- purrr::map2_dfr(.x = res_GS_sig %>% pull(MODULE),.y = res_GS_sig %>% pull(Cluster),.f = function(.x,.y) {
      mod <- .x %>% str_remove("ME")
      gsvsmm <-
        hubgenes(datExpr = datExpr,
                 mdl = mod,
                 power = power,
                 trt = .y,
                 KME = res_kme,
                 GS.cut = 0,
                 kME.cut = 0,
                 datTrait = phenotype,
                 g2m = Gene2Module
        )
      out <-
        gsvsmm$hub3 %>%
        mutate(hubGene.GS_MM = case_when(
          abs_MM > 0.8 & abs_GS > 0.75  ~ "*",
          abs_MM > 0.85 & abs_GS > 0.8  ~ "**",
          abs_MM > 0.9 & abs_GS > 0.8  ~ "**",
          abs_MM > 0.9 & abs_GS > 0.9  ~ "***",
          TRUE ~ ""
        )) %>% arrange(hubGene.GS_MM)
    })
    if(nrow(KeyModule_Trait) != 0) {
      writexl::write_xlsx(KeyModule_Trait,path = paste0(path_03,"KeyModule-trait-GSMM.xlsx"))
    } else {
      message(msg_no("There is no siginifcant correlated trait and module, please pick hubgenes by kME only!"))
    }

  }
  ##> 8.0 ModuleHeatmap
  path_04 <- paste0(out_path,"04.ExpressionModel/")
  dir.create(path = path_04,showWarnings = F,recursive = T)
  ExpdatT <- datExpr %>% scale() %>% t() %>% as.data.frame() %>%
    rownames_to_column("GID") %>%
    inner_join(Gene2Module) %>% arrange(Module)
  tmp_ht.mat <- ExpdatT %>%
    select(-Module) %>%
    column_to_rownames("GID")
  tmp_ht <- tmp_ht.mat %>%
    Heatmap(show_row_names = F,show_column_names = F)
  c_order = column_order(tmp_ht)
  ht_mat_all = tmp_ht.mat %>% select(c_order)
  module = Gene2Module %>% filter(Module != "grey") %>% pull(Module) %>%  unique()

  ordered_mat <- purrr::map_dfr(.x = module,.f = function(.x){
    mode_mat <-
      Gene2Module %>%
      filter(Module == .x) %>%
      left_join(ht_mat_all %>% rownames_to_column("GID")) %>%
      select(-Module) %>%
      column_to_rownames("GID")
    tmp2_ht <-
      Heatmap(mode_mat,show_row_names = F,show_column_names = F,cluster_columns = T,column_title = .x)
    pdf(file = paste0(path_04,.x,"_ExpProfile.pdf"),width = 10,height = 6)
    draw(tmp2_ht)
    dev.off()
    r_order = row_order(tmp2_ht)
    out_mat <- mode_mat %>% t() %>% as.data.frame() %>% select(r_order) %>% t() %>% as.data.frame()
    return(out_mat)
  })
  r_order_new <-
    data.frame(
      GID = rownames(ordered_mat)
    ) %>% left_join(Gene2Module)
  left_anno <- rowAnnotation(
    Module = r_order_new %>%
      column_to_rownames("GID") %>% as.matrix(),
    col = list(
      Module = module %>% setNames(module)
    )
  )
  pdf(file = paste0(path_04,"ExpressionHeatmap.pdf"),height = 12,width = 12)
  Heatmap(ordered_mat,show_row_names = F,show_column_names = F,cluster_columns = F,cluster_rows = F,
          left_annotation = left_anno)
  dev.off()
  ##> 9.0 export cytoscape file
  path_05 <- paste0(out_path,"05.Network/")
  dir.create(path = path_05,showWarnings = F,recursive = T)
  purrr::map(.x = module,.f = function(.x){
    cyto_out <- cytoscapeout(datExpr = datExpr,power = power,module = .x,moduleColors = moduleColors,threshold = weight_cutoff)
    node_out <- cyto_out$nodeData %>% select(1,3) %>% setNames(c("ID","Module"))
    edge_out <- cyto_out$edgeData %>% select(1:4)
    write_csv(node_out,file = paste0(path_05,.x,"_node.csv"))
    write_csv(edge_out,file = paste0(path_05,.x,"_edge.csv"))
  })
}
