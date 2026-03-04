#' MANGO multi-case circular plot (tree-level)
#'
#' Draw a circular (polar) summary plot of MANGO trees across multiple cases.
#'
#' @param filepath Directory path containing input files.
#' @param DEG_list_name DEG list table filename under \code{filepath}.
#' @param status One of \code{"common"}, \code{"specific"}, \code{"DA"}.
#' @param trends One of \code{"UP"}, \code{"DOWN"}, \code{"SIG"}.
#' @param LABEL Character vector of case labels.
#' @param project Project label shown on the plot.
#' @param condition Condition index (used in \code{"specific"}).
#' @param tree_size Text size for tree labels.
#' @param label_size Text size for individual labels.
#' @param project_size Text size for project label.
#' @param width Plot width (used by repr.plot.width option in notebooks).
#' @param height Plot height (used by repr.plot.height option in notebooks).
#'
#' @return A ggplot object.
#' @export
#'
#' @importFrom dplyr full_join inner_join anti_join group_by ungroup mutate summarize rowwise select lag n_distinct
#' @importFrom ggplot2 ggplot geom_bar geom_segment geom_point geom_text scale_size_continuous annotate ylim theme_minimal theme coord_polar scale_fill_gradientn aes element_blank
#' @importFrom stringr str_wrap
#' @importFrom magrittr %>%
#' @importFrom utils read.table read.delim
#' @importFrom stats na.omit
MANGO_TREE_cirPLOT_forMULTI = function(filepath,
                                       DEG_list_name,
                                       status,
                                       trends,
                                       LABEL,
                                       project = "plot\nname",
                                       condition = 1,
                                       tree_size,
                                       label_size,
                                       project_size,
                                       width,
                                       height){

  if(status == "common"){

    blacklist = c("regulation","positive","negative","a","the","of","gene","type","I","II","to","by","an")
    blacklist = as.data.frame(blacklist)
    colnames(blacklist) = "Description"

    DEG_list_path = paste0(filepath,"/",DEG_list_name)
    DEG_list_pathDF = as.data.frame(DEG_list_path)

    TXTbased = nrow(as.data.frame(DEG_list_pathDF[grepl("txt", DEG_list_pathDF[,1]),]))
    CSVbased = nrow(as.data.frame(DEG_list_pathDF[grepl("csv", DEG_list_pathDF[,1]),]))
    BEDbased = nrow(as.data.frame(DEG_list_pathDF[grepl("bed", DEG_list_pathDF[,1]),]))

    if(TXTbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(DEG_list_pathDF[,1], split = "_list.txt"))[1,1],"_GO_list.txt")}
    if(CSVbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(DEG_list_pathDF[,1], split = "_list.csv"))[1,1],"_GO_list.csv")}
    if(BEDbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(DEG_list_pathDF[,1], split = "_list.bed"))[1,1],"_GO_list.bed")}

    GO_list_path = OUTPUTpath

    if(TXTbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(DEG_list_pathDF[,1], split = "input_DEG_list.txt"))[1,1],"MANGO_PREPROCESSING_list.txt")}
    if(CSVbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(DEG_list_pathDF[,1], split = "input_DEG_list.csv"))[1,1],"MANGO_PREPROCESSING_list.csv")}
    if(BEDbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(DEG_list_pathDF[,1], split = "input_DEG_list.bed"))[1,1],"MANGO_PREPROCESSING_list.bed")}

    MANGO_PREPROCESSING_list_path = OUTPUTpath
    MANGO_PREPROCESSING_list_pathDF = as.data.frame(MANGO_PREPROCESSING_list_path)

    if(TXTbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.txt"))[1,1],"MANGO_SEPERATE.txt")}
    if(CSVbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.csv"))[1,1],"MANGO_SEPERATE.csv")}
    if(BEDbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.bed"))[1,1],"MANGO_SEPERATE.bed")}
    MANGO_SEPERATE_path = OUTPUTpath

    if(TXTbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.txt"))[1,1],"MANGO_SEPERATE_forMULTI_range.txt")}
    if(CSVbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.csv"))[1,1],"MANGO_SEPERATE_forMULTI_range.csv")}
    if(BEDbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.bed"))[1,1],"MANGO_SEPERATE_forMULTI_range.bed")}
    MANGO_SEPERATE_forMULTI_range_path = OUTPUTpath

    fileTABLE_path = MANGO_PREPROCESSING_list_path
    GOTABLE_path = GO_list_path

    fileTABLE = utils::read.table(fileTABLE_path, sep="\t", header=TRUE, fill=TRUE)

    MANGO_SEPERATE_INPUT = MANGO_SEPERATE_path
    SEPERATE_INPUT = utils::read.table(MANGO_SEPERATE_INPUT, sep="\t", header=TRUE, fill=TRUE)
    SEPERATE_INPUT = SEPERATE_INPUT[order(-SEPERATE_INPUT[,2]),]
    SEPERATE_INPUT = stats::na.omit(SEPERATE_INPUT)

    fileTABLE = utils::read.table(fileTABLE_path, sep="\t", header=TRUE, fill=TRUE)
    GOTABLE  = utils::read.table(GOTABLE_path,  sep="\t", header=TRUE, fill=TRUE)

    LABEL_DF = as.data.frame(LABEL)
    condition = nrow(LABEL_DF) +1

    SEPERATE_INPUT = SEPERATE_INPUT[which(SEPERATE_INPUT[,ncol(SEPERATE_INPUT)] == condition),]

    SEPERATE_INPUT_RF_mean = SEPERATE_INPUT[,c(-1,-ncol(SEPERATE_INPUT))]
    SEPERATE_INPUT_RF_mean$mean = 0

    for(i in 1:nrow(LABEL_DF)){
      SEPERATE_INPUT_RF_mean[,ncol(SEPERATE_INPUT_RF_mean)] =
        SEPERATE_INPUT_RF_mean[,i] + SEPERATE_INPUT_RF_mean[,ncol(SEPERATE_INPUT_RF_mean)]
    }

    SEPERATE_INPUT_RF_mean = SEPERATE_INPUT_RF_mean[,3]
    SEPERATE_INPUT_RF_mean = as.data.frame(SEPERATE_INPUT_RF_mean)
    SEPERATE_INPUT = cbind(SEPERATE_INPUT,SEPERATE_INPUT_RF_mean)

    SEPERATE_INPUT = SEPERATE_INPUT[order(-SEPERATE_INPUT[,2]),]
    SEPERATE_INPUT = SEPERATE_INPUT[which(SEPERATE_INPUT[,ncol(SEPERATE_INPUT)] > 0),-ncol(SEPERATE_INPUT)]

    i = 1

    GO_RAW = utils::read.delim(GOTABLE[i,1], sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
    GO_DF  = GO_RAW[,c(2,5)]
    colnames(GO_DF)[2] = paste0("RF_",i)
    GO_DF[is.na(GO_DF)] = 0

    for(i in 2:nrow(GOTABLE)){
      GO_RAW = utils::read.delim(GOTABLE[i,1], sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
      GO_DF_frag  = GO_RAW[,c(2,5)]
      colnames(GO_DF_frag)[2] = paste0("RF_",i)
      GO_DF_frag[is.na(GO_DF_frag)] = 0

      GO_DF = full_join(GO_DF,GO_DF_frag,by = "Description")
      GO_DF[is.na(GO_DF)] = 0
    }

    i = 1

    SEPERATE_INPUT_DF = as.data.frame(strsplit(SEPERATE_INPUT[i,1], split = "//"))
    colnames(SEPERATE_INPUT_DF) = "Description"

    for(s in 1:nrow(SEPERATE_INPUT_DF)){
      SEPERATE_INPUT_DF[s,1] = gsub("_", " ", SEPERATE_INPUT_DF[s,1])
    }

    DATA_RAW_DF = inner_join(SEPERATE_INPUT_DF,GO_DF,by = "Description")

    termset = DATA_RAW_DF[,1]
    termset = as.data.frame(termset)

    termset_frag = as.data.frame(strsplit(termset[1,1], split = " "))
    colnames(termset_frag) = "Description"

    for(m in 2:nrow(termset)){
      termset_fragfrag = as.data.frame(strsplit(termset[m,1], split = " "))
      colnames(termset_fragfrag) = "Description"
      termset_frag = rbind(termset_frag,termset_fragfrag)
    }

    termset_frag = anti_join(termset_frag,blacklist,by="Description")
    termset_fraguniq = termset_frag
    termset_fraguniq$hit = 0

    termset_fraguniq = unique(termset_fraguniq)

    for(s in 1:nrow(termset_fraguniq)){
      termset_fraguniq[s,2] = nrow(as.data.frame(termset_frag[which(termset_frag[,1] == termset_fraguniq[s,1]),]))
    }

    CN = paste(termset_fraguniq[1,1],termset_fraguniq[2,1],sep=" ")

    SEPERATE_INPUT[i,1] = CN

    for(i in 2:nrow(SEPERATE_INPUT)){
      SEPERATE_INPUT_DF = as.data.frame(strsplit(SEPERATE_INPUT[i,1], split = "//"))
      colnames(SEPERATE_INPUT_DF) = "Description"

      for(s in 1:nrow(SEPERATE_INPUT_DF)){
        SEPERATE_INPUT_DF[s,1] = gsub("_", " ", SEPERATE_INPUT_DF[s,1])
      }

      DATA_RAW_DF = inner_join(SEPERATE_INPUT_DF,GO_DF,by = "Description")

      termset = DATA_RAW_DF[,1]
      termset = as.data.frame(termset)

      termset_frag = as.data.frame(strsplit(termset[1,1], split = " "))
      colnames(termset_frag) = "Description"

      for(m in 2:nrow(termset)){
        termset_fragfrag = as.data.frame(strsplit(termset[m,1], split = " "))
        colnames(termset_fragfrag) = "Description"
        termset_frag = rbind(termset_frag,termset_fragfrag)
      }

      termset_frag = anti_join(termset_frag,blacklist,by="Description")
      termset_fraguniq = termset_frag
      termset_fraguniq$hit = 0

      termset_fraguniq = unique(termset_fraguniq)

      for(s in 1:nrow(termset_fraguniq)){
        termset_fraguniq[s,2] = nrow(as.data.frame(termset_frag[which(termset_frag[,1] == termset_fraguniq[s,1]),]))
      }

      CN = paste(termset_fraguniq[1,1],termset_fraguniq[2,1],sep=" ")

      SEPERATE_INPUT[i,1] = CN
    }

    SEPERATE_INPUT_fit = SEPERATE_INPUT[,c(-1,-ncol(SEPERATE_INPUT))]

    DF_tmp = cbind("","",0)
    colnames(DF_tmp) = c("individual","group","value")
    DF_tmp = DF_tmp[-1,]
    DF_tmp = as.data.frame(DF_tmp)
    DF_tmp[,3] = as.numeric(DF_tmp[,3])

    for(i in 1:nrow(SEPERATE_INPUT)){
      for(s in 1:ncol(SEPERATE_INPUT_fit)){
        DF_tmp_frag = cbind(LABEL[s],SEPERATE_INPUT[i,1],SEPERATE_INPUT_fit[i,s])
        DF_tmp_frag = as.data.frame(DF_tmp_frag)
        DF_tmp_frag[,3] = as.numeric(DF_tmp_frag[,3])
        colnames(DF_tmp_frag) = c("individual","group","value")

        DF_tmp = rbind(DF_tmp,DF_tmp_frag)
      }

      DF_tmp_frag = cbind("",SEPERATE_INPUT[i,1],0)
      colnames(DF_tmp_frag) = c("individual","group","value")
      DF_tmp_frag = as.data.frame(DF_tmp_frag)
      DF_tmp_frag[,3] = as.numeric(DF_tmp_frag[,3])

      DF_tmp = rbind(DF_tmp,DF_tmp_frag)
      DF_tmp = rbind(DF_tmp,DF_tmp_frag)
    }

    DFcir = DF_tmp
    DFcir$id = c(1:nrow(DFcir))

    DFcir <- DFcir %>%
      group_by(group) %>%
      mutate(
        sep = (individual == "" & value == 0),
        sep_n = cumsum(sep),
        boundary = sep & (sep_n %% 2 == 0),
        block_id = cumsum(dplyr::lag(boundary, default=0)) + 1
      ) %>%
      ungroup() %>%
      group_by(group) %>%
      mutate(
        n_blocks = n_distinct(block_id),
        group = ifelse(n_blocks == 1, group, paste0(group, " tree", block_id))
      ) %>%
      ungroup() %>%
      select(-sep, -sep_n, -boundary, -block_id, -n_blocks)

    label_data <- DFcir
    number_of_bar <- nrow(label_data)
    angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar
    label_data$hjust <- ifelse( angle < -90, 1, 0)
    label_data$angle <- ifelse(angle < -90, angle+180, angle)

    empty_bar = 2
    base_data <- DFcir %>%
      group_by(group) %>%
      summarize(start=min(id), end=max(id) - empty_bar) %>%
      rowwise() %>%
      mutate(title=mean(c(start, end)))

    grid_data <- base_data
    grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
    grid_data$start <- grid_data$start - 1
    grid_data <- grid_data[-1,]

    MAX = ceiling(ceiling(max(DFcir[,3]))/10)*10
    SET = MAX/2
    MAX = SET + MAX

    terms = DFcir[,2]
    terms = as.data.frame(terms)
    terms = unique(terms)
    terms$num = 0

    SEPERATE_INPUT = utils::read.table(MANGO_SEPERATE_INPUT, sep="\t", header=TRUE, fill=TRUE)
    SEPERATE_INPUT = SEPERATE_INPUT[order(-SEPERATE_INPUT[,2]),]
    SEPERATE_INPUT = stats::na.omit(SEPERATE_INPUT)

    fileTABLE = utils::read.table(fileTABLE_path, sep="\t", header=TRUE, fill=TRUE)
    GOTABLE  = utils::read.table(GOTABLE_path,  sep="\t", header=TRUE, fill=TRUE)

    SEPERATE_INPUT = SEPERATE_INPUT[which(SEPERATE_INPUT[,ncol(SEPERATE_INPUT)] == condition),]

    LABEL_DF = as.data.frame(LABEL)

    SEPERATE_INPUT_RF_mean = SEPERATE_INPUT[,c(-1,-ncol(SEPERATE_INPUT))]
    SEPERATE_INPUT_RF_mean$mean = 0

    for(i in 1:nrow(LABEL_DF)){
      SEPERATE_INPUT_RF_mean[,ncol(SEPERATE_INPUT_RF_mean)] =
        SEPERATE_INPUT_RF_mean[,i] + SEPERATE_INPUT_RF_mean[,ncol(SEPERATE_INPUT_RF_mean)]
    }

    SEPERATE_INPUT_RF_mean = SEPERATE_INPUT_RF_mean[,3]
    SEPERATE_INPUT_RF_mean = as.data.frame(SEPERATE_INPUT_RF_mean)
    SEPERATE_INPUT = cbind(SEPERATE_INPUT,SEPERATE_INPUT_RF_mean)

    SEPERATE_INPUT = SEPERATE_INPUT[order(-SEPERATE_INPUT[,2]),]

    SEPERATE_INPUT = SEPERATE_INPUT[which(SEPERATE_INPUT[,ncol(SEPERATE_INPUT)] > 0),-ncol(SEPERATE_INPUT)]

    fin = nrow(LABEL_DF) -1

    BB = as.data.frame(cbind("X",0))

    for(q in 1:fin){
      BB_frag = as.data.frame(cbind("X",0))
      BB = rbind(BB,BB_frag)
    }

    for(q in 1:2){
      BB_frag = as.data.frame(cbind("X",0))
      BB = rbind(BB,BB_frag)
    }

    BB[,2] = as.numeric(BB[,2])
    colnames(BB) = c("group","num")

    m = 1

    SEPERATE_INPUT_DF = as.data.frame(strsplit(SEPERATE_INPUT[m,1], split = "//"))
    colnames(SEPERATE_INPUT_DF) = "Description"

    for(s in 1:nrow(SEPERATE_INPUT_DF)){
      SEPERATE_INPUT_DF[s,1] = gsub("_", " ", SEPERATE_INPUT_DF[s,1])
    }

    DFcir_add = BB
    DFcir_add[,1] = terms[m,1]

    for(i in 1:nrow(LABEL_DF)){
      GO_RAW = utils::read.delim(GOTABLE[i,1], sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
      GO_DF  = GO_RAW[,c(2,5)]
      colnames(GO_DF)[2] = paste0("RF_",i)
      GO_DF[is.na(GO_DF)] = 0

      DATA_RAW_DF = inner_join(SEPERATE_INPUT_DF,GO_DF,by = "Description")
      DFcir_add[i,2] = nrow(DATA_RAW_DF)
    }

    for(m in 2:nrow(terms)){
      SEPERATE_INPUT_DF = as.data.frame(strsplit(SEPERATE_INPUT[m,1], split = "//"))
      colnames(SEPERATE_INPUT_DF) = "Description"

      for(s in 1:nrow(SEPERATE_INPUT_DF)){
        SEPERATE_INPUT_DF[s,1] = gsub("_", " ", SEPERATE_INPUT_DF[s,1])
      }

      DFcir_add_frag = BB
      DFcir_add_frag[,1] = terms[m,1]

      for(i in 1:nrow(LABEL_DF)){
        GO_RAW = utils::read.delim(GOTABLE[i,1], sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
        GO_DF  = GO_RAW[,c(2,5)]
        colnames(GO_DF)[2] = paste0("RF_",i)
        GO_DF[is.na(GO_DF)] = 0

        DATA_RAW_DF = inner_join(SEPERATE_INPUT_DF,GO_DF,by = "Description")
        DFcir_add_frag[i,2] = nrow(DATA_RAW_DF)
      }

      DFcir_add = rbind(DFcir_add,DFcir_add_frag)
    }

    DFcir_add = DFcir_add[,2]
    DFcir_add = as.data.frame(DFcir_add)
    colnames(DFcir_add) = "num"

    DFcir = cbind(DFcir,DFcir_add)

    DEG_LIST = utils::read.delim(DEG_list_pathDF[1,1], sep = "\t", header = FALSE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)

    i = 1

    DEG = utils::read.delim(DEG_LIST[i,1], sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
    DEG = DEG[,c(1,3)]
    DEG[,2] = abs(DEG[,2])
    DEG[,2] = 2^(DEG[,2])
    colnames(DEG)[2] = paste0("RF_",i)

    for(i in 2:nrow(LABEL_DF)){
      DEG_frag = utils::read.delim(DEG_LIST[i,1], sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
      DEG_frag = DEG_frag[,c(1,3)]
      DEG_frag[,2] = abs(DEG_frag[,2])
      DEG_frag[,2] = 2^(DEG_frag[,2])
      colnames(DEG_frag)[2] = paste0("RF_",i)

      DEG = full_join(DEG,DEG_frag,by = "gene")
      DEG[is.na(DEG)] = 0
    }

    DEG$median = DEG[,2]
    for(i in 2:nrow(LABEL_DF)){
      DEG$median = DEG[,(i+1)] + DEG$median
    }

    DEG$median = DEG$median/nrow(LABEL_DF)

    DEG = DEG[,c(1,ncol(DEG))]

    i = 1

    GO_RAW = utils::read.delim(GOTABLE[i,1], sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
    GO_DF  = GO_RAW[,c(2,11)]
    colnames(GO_DF)[2] = paste0("RF_",i)

    for(i in 2:nrow(LABEL_DF)){
      GO_RAW = utils::read.delim(GOTABLE[i,1], sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
      GO_DF_frag  = GO_RAW[,c(2,11)]
      colnames(GO_DF_frag)[2] = paste0("RF_",i)

      GO_DF = full_join(GO_DF,GO_DF_frag,by = "Description")
    }

    GO_DF$FULL = GO_DF[,2]

    for(i in 2:nrow(LABEL_DF)){
      GO_DF$FULL = paste(GO_DF[,(i+1)],GO_DF[,ncol(GO_DF)],sep="/")
    }

    GO_DF = GO_DF[,c(1,ncol(GO_DF))]

    terms$genenum = 0
    terms$geneS = 0

    for(i in 1:nrow(SEPERATE_INPUT)){

      SEPERATE_INPUT_DF = as.data.frame(strsplit(SEPERATE_INPUT[i,1], split = "//"))
      colnames(SEPERATE_INPUT_DF) = "Description"

      for(s in 1:nrow(SEPERATE_INPUT_DF)){
        SEPERATE_INPUT_DF[s,1] = gsub("_", " ", SEPERATE_INPUT_DF[s,1])
      }

      DATA_RAW_DF = inner_join(SEPERATE_INPUT_DF,GO_DF,by = "Description")

      geneset = DATA_RAW_DF[,2]
      geneset = as.data.frame(geneset)
      colnames(geneset) = "gene"

      geneset_s = "NOT"
      geneset_s = as.data.frame(geneset_s)
      colnames(geneset_s) = "gene"
      geneset_s = geneset_s[-1,]

      for(w in 1:nrow(geneset)){
        geneset_s_frag = as.data.frame(strsplit(geneset[w,1], split = "/"))
        colnames(geneset_s_frag) = "gene"
        geneset_s = rbind(geneset_s,geneset_s_frag)
        geneset_s = unique(geneset_s)
      }

      geneset_s = inner_join(geneset_s,DEG,by="gene")

      terms[i,3] = nrow(geneset_s)/nrow(DEG)*100
      terms[i,4] = median(geneset_s[,2])
    }

    ratio_max = ceiling(ceiling(max(terms[,3]))/10)*10

    R_max = -SET*0.7
    R_min = -SET*1.3

    for(i in 1:nrow(terms)){
      terms[i,3] = ((terms[i,3]/10)*(R_max-R_min))+R_min
    }

    base_data = inner_join(base_data,terms,by="group")

    terms$order = c(1:nrow(terms))

    base_data_add = terms[,c(1,ncol(terms))]

    base_data = inner_join(base_data,base_data_add,by="group")
    base_data$project = ""
    base_data[1,ncol(base_data)]=project

    if(trends == "UP"){
      basecol = "#DD686CFF"
      fillcol = "pink2"
      fontcol = "darkred"
    }

    if(trends == "DOWN"){
      basecol = "#33E4FFFF"
      fillcol = "skyblue2"
      fontcol = "darkblue"
    }

    if(trends == "SIG"){
      basecol = "tan"
      fillcol = "tan2"
      fontcol = "brown"
    }

    if(trends == "UP"){colorset = c("#F5B78EFF", "#F19C7CFF", "#EA8171FF", "#DD686CFF", "#CA5268FF", "#B13F64FF")}
    if(trends == "DOWN"){colorset = c("#CCFDFFFF","#99F8FFFF","#66F0FFFF","#33E4FFFF","#00AACCFF","#007A99FF")}

    plot_multi = ggplot(DFcir) +
      geom_bar(aes(x=as.factor(id), y=value, fill=num), stat="identity", alpha=1) +
      geom_segment(data=grid_data, aes(x = end-1, y = SET*2, xend = 0, yend = SET*2), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
      geom_segment(data=grid_data, aes(x = end-1, y = SET, xend = 0, yend = SET), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
      geom_segment(data=grid_data, aes(x = end-1, y = 0, xend = 0, yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +

      geom_segment(data=base_data, aes(x = title, y = 0, xend = title, yend = SET*2), linetype = "dashed", color = "gray40") +
      geom_segment(data=base_data, aes(x = title, y = -SET*1.3, xend = title, yend = -SET*0.7), linetype = "dashed", color = "gray40") +

      geom_segment(data=grid_data, aes(x = end-1, y = -SET*1.3, xend = 0, yend = -SET*1.3), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
      geom_segment(data=grid_data, aes(x = end-1, y = -SET*1, xend = 0, yend = -SET*1), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
      geom_segment(data=grid_data, aes(x = end-1, y = -SET*0.7, xend = 0, yend = -SET*0.7), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +

      geom_bar(aes(x=as.factor(id), y=value, fill=num), stat="identity", alpha=1) +

      geom_segment(data=base_data, aes(x = start, y = 0, xend = end, yend = 0), colour = "black", alpha=1, size=1 , inherit.aes = FALSE ) +
      geom_segment(data=base_data, aes(x = start, y = -SET*0.7, xend = end, yend = -SET*0.7), colour = "black", alpha=1, size=1 , inherit.aes = FALSE ) +
      geom_segment(data=base_data, aes(x = start, y = -SET*1.3, xend = end, yend = -SET*1.3), colour = "black", alpha=1, size=1 , inherit.aes = FALSE ) +

      geom_point(data=base_data, aes(x = title, y = genenum, size = geneS), color = fontcol,show.legend = TRUE) +

      scale_size_continuous("FC of genes",range = c(1, 5)) +

      annotate("text", x = rep(max(DFcir$id)+1,4), y = c(0,SET,SET*2,SET*(2.2)), label = c("0",SET ,SET*2,"HWES") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +

      annotate("text", x = rep(max(DFcir$id)-1.5,4), y = c(-SET*1.3,-SET*1,-SET*0.7,-SET*0.55), label = c("0%",paste0(ratio_max/2,"%"),paste0(ratio_max,"%"),"DEG ratio") , color="grey", size=2 , angle=0, fontface="bold", hjust=0) +

      ylim(-MAX*(0.9),MAX*(1.05)) +
      theme_minimal() +
      theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank()
      ) +
      coord_polar() +
      scale_fill_gradientn("Number of terms", colours = colorset) +

      geom_text(data=label_data, aes(x=id, y=-SET*(0.6), label=individual, hjust=hjust),
                color=fillcol, fontface="bold",alpha=1, size=label_size, angle= label_data$angle, inherit.aes = FALSE ) +

      geom_text(data=base_data, aes(x = title, y = MAX*(1.05), label=str_wrap(group, 5)),
                colour = "black", alpha=1, size=tree_size, fontface="bold", inherit.aes = FALSE) +
      geom_text(data=base_data, aes(x = title, y = -SET*1.6, label=str_wrap(order, 5)),
                colour = basecol, alpha=1, size=5, fontface="bold", inherit.aes = FALSE) +
      geom_text(data=base_data, aes(x = title, y = -MAX*(0.9), label=str_wrap(project, 5)),
                colour = fontcol, alpha=1, size=project_size, fontface="bold", inherit.aes = FALSE)

  }

    if(status == "specific"){

        
project = LABEL[condition]

blacklist = c("regulation","positive","negative","a","the","of","gene","type","I","II","to","by","an")
blacklist = as.data.frame(blacklist)
colnames(blacklist) = "Description"

    DEG_list_path = paste0(filepath,"/",DEG_list_name)
    DEG_list_pathDF = as.data.frame(DEG_list_path)
    
    TXTbased = nrow(as.data.frame(DEG_list_pathDF[grepl("txt", DEG_list_pathDF[,1]),]))
    CSVbased = nrow(as.data.frame(DEG_list_pathDF[grepl("csv", DEG_list_pathDF[,1]),]))
    BEDbased = nrow(as.data.frame(DEG_list_pathDF[grepl("bed", DEG_list_pathDF[,1]),]))

    if(TXTbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(DEG_list_pathDF[,1], split = "_list.txt"))[1,1],"_GO_list.txt")}
    if(CSVbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(DEG_list_pathDF[,1], split = "_list.csv"))[1,1],"_GO_list.csv")}
    if(BEDbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(DEG_list_pathDF[,1], split = "_list.bed"))[1,1],"_GO_list.bed")}

    GO_list_path = OUTPUTpath

    if(TXTbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(DEG_list_pathDF[,1], split = "input_DEG_list.txt"))[1,1],"MANGO_PREPROCESSING_list.txt")}
    if(CSVbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(DEG_list_pathDF[,1], split = "input_DEG_list.csv"))[1,1],"MANGO_PREPROCESSING_list.csv")}
    if(BEDbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(DEG_list_pathDF[,1], split = "input_DEG_list.bed"))[1,1],"MANGO_PREPROCESSING_list.bed")}

    MANGO_PREPROCESSING_list_path = OUTPUTpath
    MANGO_PREPROCESSING_list_pathDF = as.data.frame(MANGO_PREPROCESSING_list_path)

    if(TXTbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.txt"))[1,1],"MANGO_SEPERATE.txt")}
    if(CSVbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.csv"))[1,1],"MANGO_SEPERATE.csv")}
    if(BEDbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.bed"))[1,1],"MANGO_SEPERATE.bed")}
    MANGO_SEPERATE_path = OUTPUTpath

    if(TXTbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.txt"))[1,1],"MANGO_SEPERATE_forMULTI_range.txt")}
    if(CSVbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.csv"))[1,1],"MANGO_SEPERATE_forMULTI_range.csv")}
    if(BEDbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.bed"))[1,1],"MANGO_SEPERATE_forMULTI_range.bed")}
    MANGO_SEPERATE_forMULTI_range_path = OUTPUTpath

    fileTABLE_path = MANGO_PREPROCESSING_list_path
    GOTABLE_path = GO_list_path

    fileTABLE = utils::read.table(fileTABLE_path, sep="\t", header=TRUE, fill=TRUE)

    MANGO_SEPERATE_INPUT = MANGO_SEPERATE_path
    SEPERATE_INPUT = utils::read.table(MANGO_SEPERATE_INPUT, sep="\t", header=TRUE, fill=TRUE)
    SEPERATE_INPUT = SEPERATE_INPUT[order(-SEPERATE_INPUT[,2]),]
    SEPERATE_INPUT = stats::na.omit(SEPERATE_INPUT)

    fileTABLE = utils::read.table(fileTABLE_path, sep="\t", header=TRUE, fill=TRUE)
    GOTABLE  = utils::read.table(GOTABLE_path,  sep="\t", header=TRUE, fill=TRUE)

    SEPERATE_INPUT = SEPERATE_INPUT[which(SEPERATE_INPUT[,ncol(SEPERATE_INPUT)] == condition),]
    SEPERATE_INPUT = SEPERATE_INPUT[order(-SEPERATE_INPUT[,(condition+1)]),]

 i = 1

    GO_RAW = utils::read.delim(GOTABLE[i,1], sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
    GO_DF  = GO_RAW[,c(2,5)]
    colnames(GO_DF)[2] = paste0("RF_",i)
    GO_DF[is.na(GO_DF)] = 0

    for(i in 2:nrow(GOTABLE)){
        GO_RAW = utils::read.delim(GOTABLE[i,1], sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
        GO_DF_frag  = GO_RAW[,c(2,5)]
        colnames(GO_DF_frag)[2] = paste0("RF_",i)
        GO_DF_frag[is.na(GO_DF_frag)] = 0
        
        GO_DF = full_join(GO_DF,GO_DF_frag,by = "Description")
        GO_DF[is.na(GO_DF)] = 0
        }
    

    i = 1

    SEPERATE_INPUT_DF = as.data.frame(strsplit(SEPERATE_INPUT[i,1], split = "//"))
    colnames(SEPERATE_INPUT_DF) = "Description"

    for(s in 1:nrow(SEPERATE_INPUT_DF)){
            SEPERATE_INPUT_DF[s,1] = gsub("_", " ", SEPERATE_INPUT_DF[s,1])
        }

    DATA_RAW_DF = inner_join(SEPERATE_INPUT_DF,GO_DF,by = "Description")

    termset = DATA_RAW_DF[,1]
    termset = as.data.frame(termset)

    termset_frag = as.data.frame(strsplit(termset[1,1], split = " "))
    colnames(termset_frag) = "Description"

    for(m in 2:nrow(termset)){
        termset_fragfrag = as.data.frame(strsplit(termset[m,1], split = " "))
        colnames(termset_fragfrag) = "Description"
        termset_frag = rbind(termset_frag,termset_fragfrag)
    }

    termset_frag = anti_join(termset_frag,blacklist,by="Description")
    termset_fraguniq = termset_frag
    termset_fraguniq$hit = 0

    termset_fraguniq = unique(termset_fraguniq)

    for(s in 1:nrow(termset_fraguniq)){
        termset_fraguniq[s,2] = nrow(as.data.frame(termset_frag[which(termset_frag[,1] == termset_fraguniq[s,1]),]))
        }

    CN = paste(termset_fraguniq[1,1],termset_fraguniq[2,1],sep=" ")

    SEPERATE_INPUT[i,1] = CN


    for(i in 2:nrow(SEPERATE_INPUT)){
    SEPERATE_INPUT_DF = as.data.frame(strsplit(SEPERATE_INPUT[i,1], split = "//"))
    colnames(SEPERATE_INPUT_DF) = "Description"

    for(s in 1:nrow(SEPERATE_INPUT_DF)){
            SEPERATE_INPUT_DF[s,1] = gsub("_", " ", SEPERATE_INPUT_DF[s,1])
        }

    DATA_RAW_DF = inner_join(SEPERATE_INPUT_DF,GO_DF,by = "Description")

    termset = DATA_RAW_DF[,1]
    termset = as.data.frame(termset)

    termset_frag = as.data.frame(strsplit(termset[1,1], split = " "))
    colnames(termset_frag) = "Description"

    for(m in 2:nrow(termset)){
        termset_fragfrag = as.data.frame(strsplit(termset[m,1], split = " "))
        colnames(termset_fragfrag) = "Description"
        termset_frag = rbind(termset_frag,termset_fragfrag)
    }

    termset_frag = anti_join(termset_frag,blacklist,by="Description")
    termset_fraguniq = termset_frag
    termset_fraguniq$hit = 0

    termset_fraguniq = unique(termset_fraguniq)

    for(s in 1:nrow(termset_fraguniq)){
        termset_fraguniq[s,2] = nrow(as.data.frame(termset_frag[which(termset_frag[,1] == termset_fraguniq[s,1]),]))
        }

    CN = paste(termset_fraguniq[1,1],termset_fraguniq[2,1],sep=" ")

    SEPERATE_INPUT[i,1] = CN
        }

    SEPERATE_INPUT_fit = SEPERATE_INPUT[,c(-1,-ncol(SEPERATE_INPUT))]

    DF_tmp = cbind("","",0)
    colnames(DF_tmp) = c("individual","group","value")
    DF_tmp = DF_tmp[-1,]
    DF_tmp = as.data.frame(DF_tmp)
    DF_tmp[,3] = as.numeric(DF_tmp[,3])

    for(i in 1:nrow(SEPERATE_INPUT)){
        for(s in 1:ncol(SEPERATE_INPUT_fit)){
            DF_tmp_frag = cbind(LABEL[s],SEPERATE_INPUT[i,1],SEPERATE_INPUT_fit[i,s])
            DF_tmp_frag = as.data.frame(DF_tmp_frag)
            DF_tmp_frag[,3] = as.numeric(DF_tmp_frag[,3])
            colnames(DF_tmp_frag) = c("individual","group","value")

            DF_tmp = rbind(DF_tmp,DF_tmp_frag)
            }

        DF_tmp_frag = cbind("",SEPERATE_INPUT[i,1],0)
        colnames(DF_tmp_frag) = c("individual","group","value")
        DF_tmp_frag = as.data.frame(DF_tmp_frag)
        DF_tmp_frag[,3] = as.numeric(DF_tmp_frag[,3])

        DF_tmp = rbind(DF_tmp,DF_tmp_frag)
        DF_tmp = rbind(DF_tmp,DF_tmp_frag)

        }

    DFcir = DF_tmp
    DFcir$id = c(1:nrow(DFcir))



    
DFcir <- DFcir %>%
  group_by(group) %>%
  mutate(
    sep = (individual == "" & value == 0),
    sep_n = cumsum(sep),
    boundary = sep & (sep_n %% 2 == 0),                 # 각 블록의 "두 번째 빈줄"에 boundary=TRUE
    block_id = cumsum(dplyr::lag(boundary, default=0)) + 1  # boundary는 '다음 줄'부터 블록 증가
  ) %>%
  ungroup() %>%
  group_by(group) %>%
  mutate(
    n_blocks = n_distinct(block_id),
    group = ifelse(n_blocks == 1, group, paste0(group, " tree", block_id))
  ) %>%
  ungroup() %>%
  select(-sep, -sep_n, -boundary, -block_id, -n_blocks)

label_data <- DFcir
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar  
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

empty_bar = 2
base_data <- DFcir %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))
 
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

base_data$case = ""
base_data[1,ncol(base_data)] = paste0(LABEL[condition])

MAX = ceiling(ceiling(max(DFcir[,3]))/10)*10
SET = MAX/2
MAX = SET + MAX

terms = DFcir[,2]
terms = as.data.frame(terms)
terms = unique(terms)
terms$num = 0


MANGO_SEPERATE_INPUT = MANGO_SEPERATE_path
    SEPERATE_INPUT = utils::read.table(MANGO_SEPERATE_INPUT, sep="\t", header=TRUE, fill=TRUE)
    SEPERATE_INPUT = SEPERATE_INPUT[order(-SEPERATE_INPUT[,2]),]
    SEPERATE_INPUT = stats::na.omit(SEPERATE_INPUT)
    SEPERATE_INPUT = SEPERATE_INPUT[which(SEPERATE_INPUT[,ncol(SEPERATE_INPUT)] == condition),]
    SEPERATE_INPUT = SEPERATE_INPUT[order(-SEPERATE_INPUT[,(condition+1)]),]

i = condition
    GO_RAW = utils::read.delim(GOTABLE[i,1], sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
    GO_DF  = GO_RAW[,c(2,5)]
    colnames(GO_DF)[2] = paste0("RF_",i)
    GO_DF[is.na(GO_DF)] = 0

    for(i in 1:nrow(terms)){
    SEPERATE_INPUT_DF = as.data.frame(strsplit(SEPERATE_INPUT[i,1], split = "//"))
    colnames(SEPERATE_INPUT_DF) = "Description"

    for(s in 1:nrow(SEPERATE_INPUT_DF)){
            SEPERATE_INPUT_DF[s,1] = gsub("_", " ", SEPERATE_INPUT_DF[s,1])
        }

    DATA_RAW_DF = inner_join(SEPERATE_INPUT_DF,GO_DF,by = "Description")
    terms[i,2] = nrow(DATA_RAW_DF)
        }

colnames(terms)[1] = "group"

DFcir = inner_join(DFcir,terms,by="group")

for(i in 1:nrow(DFcir)){
    if(DFcir[i,1] != LABEL[condition]){
        DFcir[i,ncol(DFcir)] = 0
        }
    }


DEG_list_DF = utils::read.delim(DEG_list_pathDF[1,1], sep = "\t", header = FALSE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
DEG = utils::read.delim(DEG_list_DF[condition,1], sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
DEG[,3] = abs(DEG[,3])

i = condition
    GO_RAW = utils::read.delim(GOTABLE[i,1], sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
    GO_DF  = GO_RAW[,c(2,5,11)]

terms$genenum = 0
terms$geneS = 0

i = 1

for(i in 1:nrow(SEPERATE_INPUT)){

SEPERATE_INPUT_DF = as.data.frame(strsplit(SEPERATE_INPUT[i,1], split = "//"))
    colnames(SEPERATE_INPUT_DF) = "Description"

    for(s in 1:nrow(SEPERATE_INPUT_DF)){
            SEPERATE_INPUT_DF[s,1] = gsub("_", " ", SEPERATE_INPUT_DF[s,1])
        }

    DATA_RAW_DF = inner_join(SEPERATE_INPUT_DF,GO_DF,by = "Description")

geneset = DATA_RAW_DF[,3]
geneset = as.data.frame(geneset)
colnames(geneset) = "gene"

geneset_s = "NOT"
geneset_s = as.data.frame(geneset_s)
colnames(geneset_s) = "gene"
geneset_s = geneset_s[-1,]

for(w in 1:nrow(geneset)){
    geneset_s_frag = as.data.frame(strsplit(geneset[w,1], split = "/"))
    colnames(geneset_s_frag) = "gene"
    geneset_s = rbind(geneset_s,geneset_s_frag)
    geneset_s = unique(geneset_s)
    }

geneset_s = inner_join(geneset_s,DEG,by="gene")

terms[i,3] = nrow(geneset_s)/nrow(DEG)*100
terms[i,4] = median(geneset_s[,3])
    }

ratio_max = ceiling(ceiling(max(terms[,3]))/10)*10

R_max = -SET*0.7
R_min = -SET*1.3

for(i in 1:nrow(terms)){
    terms[i,3] = ((terms[i,3]/10)*(R_max-R_min))+R_min
    }

base_data = inner_join(base_data,terms,by="group")


terms$order = c(1:nrow(terms))

base_data_add = terms[,c(1,ncol(terms))]

base_data = inner_join(base_data,base_data_add,by="group")
base_data$project = ""
base_data[1,ncol(base_data)]=project

    if(trends == "UP"){
        basecol = "#DD686CFF"
        fillcol = "pink2"
        fontcol = "darkred"
    }

    if(trends == "DOWN"){
        basecol = "#33E4FFFF"
        fillcol = "skyblue2"
        fontcol = "darkblue"
    }

    if(trends == "SIG"){
        basecol = "tan"
        fillcol = "tan2"
        fontcol = "brown"
    }

if(trends == "UP"){colorset = c("gray40", "#F5B78EFF", "#F19C7CFF", "#EA8171FF", "#DD686CFF", "#CA5268FF", "#B13F64FF")}
if(trends == "DOWN"){colorset = c("gray40", "#CCFDFFFF","#99F8FFFF","#66F0FFFF","#33E4FFFF","#00AACCFF","#007A99FF")}

plot_multi = ggplot(DFcir) +      
geom_bar(aes(x=as.factor(id), y=value, fill=num), stat="identity", alpha=1) +
geom_segment(data=grid_data, aes(x = end-1, y = SET*2, xend = 0, yend = SET*2), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
geom_segment(data=grid_data, aes(x = end-1, y = SET, xend = 0, yend = SET), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
geom_segment(data=grid_data, aes(x = end-1, y = 0, xend = 0, yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +

geom_segment(data=base_data, aes(x = title, y = 0, xend = title, yend = SET*2), linetype = "dashed", color = "gray40") + 
geom_segment(data=base_data, aes(x = title, y = -SET*1.3, xend = title, yend = -SET*0.7), linetype = "dashed", color = "gray40") +

geom_segment(data=grid_data, aes(x = end-1, y = -SET*1.3, xend = 0, yend = -SET*1.3), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
geom_segment(data=grid_data, aes(x = end-1, y = -SET*1, xend = 0, yend = -SET*1), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
geom_segment(data=grid_data, aes(x = end-1, y = -SET*0.7, xend = 0, yend = -SET*0.7), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +

geom_bar(aes(x=as.factor(id), y=value, fill=num), stat="identity", alpha=1) +

geom_segment(data=base_data, aes(x = start, y = 0, xend = end, yend = 0), colour = "black", alpha=1, size=1 , inherit.aes = FALSE ) +
geom_segment(data=base_data, aes(x = start, y = -SET*0.7, xend = end, yend = -SET*0.7), colour = "black", alpha=1, size=1 , inherit.aes = FALSE ) +
geom_segment(data=base_data, aes(x = start, y = -SET*1.3, xend = end, yend = -SET*1.3), colour = "black", alpha=1, size=1 , inherit.aes = FALSE ) +

geom_point(data=base_data, aes(x = title, y = genenum, size = geneS), color = fontcol,show.legend = TRUE) +

scale_size_continuous("FC of genes",range = c(1, 5)) +

annotate("text", x = rep(max(DFcir$id)+1,4), y = c(0,SET,SET*2,SET*(2.2)), label = c("0",SET ,SET*2,"HWES") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +

annotate("text", x = rep(max(DFcir$id)-1.5,4), y = c(-SET*1.3,-SET*1,-SET*0.7,-SET*0.55), label = c("0%",paste0(ratio_max/2,"%"),paste0(ratio_max,"%"),"DEG ratio") , color="grey", size=2 , angle=0, fontface="bold", hjust=0) +

ylim(-MAX*(0.9),MAX*(1.05)) +
theme_minimal() +
theme(#legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()) +
coord_polar() + 
scale_fill_gradientn("Number of terms", colours = colorset) +

geom_text(data=label_data, aes(x=id, y=-SET*(0.6), label=individual, hjust=hjust), color=fillcol, fontface="bold",alpha=1, size=label_size, angle= label_data$angle, inherit.aes = FALSE ) +

geom_text(data=base_data, aes(x = title, y = MAX*(1.05), label=str_wrap(group, 5)), colour = "black", alpha=1, size=tree_size, fontface="bold", inherit.aes = FALSE) +
geom_text(data=base_data, aes(x = title, y = -SET*1.6, label=str_wrap(order, 5)), colour = basecol, alpha=1, size=5, fontface="bold", inherit.aes = FALSE) +
geom_text(data=base_data, aes(x = title, y = -MAX*(0.9), label=str_wrap(project, 5)), colour = fontcol, alpha=1, size=project_size, fontface="bold", inherit.aes = FALSE)
        }

if(status == "DA"){
        

blacklist = c("regulation","positive","negative","a","the","of","gene","type","I","II","to","by","an")
blacklist = as.data.frame(blacklist)
colnames(blacklist) = "Description"

    DEG_list_path = paste0(filepath,"/",DEG_list_name)
    DEG_list_pathDF = as.data.frame(DEG_list_path)
    
    TXTbased = nrow(as.data.frame(DEG_list_pathDF[grepl("txt", DEG_list_pathDF[,1]),]))
    CSVbased = nrow(as.data.frame(DEG_list_pathDF[grepl("csv", DEG_list_pathDF[,1]),]))
    BEDbased = nrow(as.data.frame(DEG_list_pathDF[grepl("bed", DEG_list_pathDF[,1]),]))

    if(TXTbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(DEG_list_pathDF[,1], split = "_list.txt"))[1,1],"_GO_list.txt")}
    if(CSVbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(DEG_list_pathDF[,1], split = "_list.csv"))[1,1],"_GO_list.csv")}
    if(BEDbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(DEG_list_pathDF[,1], split = "_list.bed"))[1,1],"_GO_list.bed")}

    GO_list_path = OUTPUTpath

    if(TXTbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(DEG_list_pathDF[,1], split = "input_DEG_list.txt"))[1,1],"MANGO_PREPROCESSING_list.txt")}
    if(CSVbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(DEG_list_pathDF[,1], split = "input_DEG_list.csv"))[1,1],"MANGO_PREPROCESSING_list.csv")}
    if(BEDbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(DEG_list_pathDF[,1], split = "input_DEG_list.bed"))[1,1],"MANGO_PREPROCESSING_list.bed")}

    MANGO_PREPROCESSING_list_path = OUTPUTpath
    MANGO_PREPROCESSING_list_pathDF = as.data.frame(MANGO_PREPROCESSING_list_path)

    if(TXTbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.txt"))[1,1],"MANGO_TERMLISTING.txt")}
    if(CSVbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.csv"))[1,1],"MANGO_TERMLISTING.csv")}
    if(BEDbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.bed"))[1,1],"MANGO_TERMLISTING.bed")}
    MANGO_TERMLISTING_path = OUTPUTpath

    if(TXTbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.txt"))[1,1],"MANGO_COMPARE.txt")}
    if(CSVbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.csv"))[1,1],"MANGO_COMPARE.csv")}
    if(BEDbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.bed"))[1,1],"MANGO_COMPARE.bed")}
    MANGO_COMPARE_path = OUTPUTpath

    if(TXTbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.txt"))[1,1],"MANGO_SEPERATE.txt")}
    if(CSVbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.csv"))[1,1],"MANGO_SEPERATE.csv")}
    if(BEDbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.bed"))[1,1],"MANGO_SEPERATE.bed")}
    MANGO_SEPERATE_path = OUTPUTpath

    if(TXTbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.txt"))[1,1],"MANGO_SEPERATE_forMULTI_range.txt")}
    if(CSVbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.csv"))[1,1],"MANGO_SEPERATE_forMULTI_range.csv")}
    if(BEDbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.bed"))[1,1],"MANGO_SEPERATE_forMULTI_range.bed")}
    MANGO_SEPERATE_forMULTI_range_path = OUTPUTpath

    fileTABLE_path = MANGO_PREPROCESSING_list_path
    GOTABLE_path = GO_list_path

    fileTABLE = utils::read.table(fileTABLE_path, sep="\t", header=TRUE, fill=TRUE)

    MANGO_SEPERATE_INPUT = MANGO_SEPERATE_forMULTI_range_path
    SEPERATE_INPUT = utils::read.table(MANGO_SEPERATE_INPUT, sep="\t", header=TRUE, fill=TRUE)

    condition = SEPERATE_INPUT[1,ncol(SEPERATE_INPUT)]
    condition = eval(parse(text = condition))

    SEPERATE_INPUT = SEPERATE_INPUT[,-ncol(SEPERATE_INPUT)]
    SEPERATE_INPUT = stats::na.omit(SEPERATE_INPUT)
    SEPERATE_INPUT = SEPERATE_INPUT[order(-SEPERATE_INPUT[,(condition[1]+1)]),]

    fileTABLE = utils::read.table(fileTABLE_path, sep="\t", header=TRUE, fill=TRUE)
    GOTABLE  = utils::read.table(GOTABLE_path,  sep="\t", header=TRUE, fill=TRUE)

    i = 1

    GO_RAW = utils::read.delim(GOTABLE[i,1], sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
    GO_DF  = GO_RAW[,c(2,5)]
    colnames(GO_DF)[2] = paste0("RF_",i)
    GO_DF[is.na(GO_DF)] = 0

    for(i in 2:nrow(GOTABLE)){
        GO_RAW = utils::read.delim(GOTABLE[i,1], sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
        GO_DF_frag  = GO_RAW[,c(2,5)]
        colnames(GO_DF_frag)[2] = paste0("RF_",i)
        GO_DF_frag[is.na(GO_DF_frag)] = 0
        
        GO_DF = full_join(GO_DF,GO_DF_frag,by = "Description")
        GO_DF[is.na(GO_DF)] = 0
        }

    i = 1

    SEPERATE_INPUT_DF = as.data.frame(strsplit(SEPERATE_INPUT[i,1], split = "//"))
    colnames(SEPERATE_INPUT_DF) = "Description"

    for(s in 1:nrow(SEPERATE_INPUT_DF)){
            SEPERATE_INPUT_DF[s,1] = gsub("_", " ", SEPERATE_INPUT_DF[s,1])
        }

    DATA_RAW_DF = inner_join(SEPERATE_INPUT_DF,GO_DF,by = "Description")

    termset = DATA_RAW_DF[,1]
    termset = as.data.frame(termset)

    termset_frag = as.data.frame(strsplit(termset[1,1], split = " "))
    colnames(termset_frag) = "Description"

    for(m in 2:nrow(termset)){
        termset_fragfrag = as.data.frame(strsplit(termset[m,1], split = " "))
        colnames(termset_fragfrag) = "Description"
        termset_frag = rbind(termset_frag,termset_fragfrag)
    }

    termset_frag = anti_join(termset_frag,blacklist,by="Description")
    termset_fraguniq = termset_frag
    termset_fraguniq$hit = 0

    termset_fraguniq = unique(termset_fraguniq)

    for(s in 1:nrow(termset_fraguniq)){
        termset_fraguniq[s,2] = nrow(as.data.frame(termset_frag[which(termset_frag[,1] == termset_fraguniq[s,1]),]))
        }

    CN = paste(termset_fraguniq[1,1],termset_fraguniq[2,1],sep=" ")

    SEPERATE_INPUT[i,1] = CN

    if(nrow(SEPERATE_INPUT)>1){

    for(i in 2:nrow(SEPERATE_INPUT)){
    SEPERATE_INPUT_DF = as.data.frame(strsplit(SEPERATE_INPUT[i,1], split = "//"))
    colnames(SEPERATE_INPUT_DF) = "Description"

    for(s in 1:nrow(SEPERATE_INPUT_DF)){
            SEPERATE_INPUT_DF[s,1] = gsub("_", " ", SEPERATE_INPUT_DF[s,1])
        }

    DATA_RAW_DF = inner_join(SEPERATE_INPUT_DF,GO_DF,by = "Description")

    termset = DATA_RAW_DF[,1]
    termset = as.data.frame(termset)

    termset_frag = as.data.frame(strsplit(termset[1,1], split = " "))
    colnames(termset_frag) = "Description"

    for(m in 2:nrow(termset)){
        termset_fragfrag = as.data.frame(strsplit(termset[m,1], split = " "))
        colnames(termset_fragfrag) = "Description"
        termset_frag = rbind(termset_frag,termset_fragfrag)
    }

    termset_frag = anti_join(termset_frag,blacklist,by="Description")
    termset_fraguniq = termset_frag
    termset_fraguniq$hit = 0

    termset_fraguniq = unique(termset_fraguniq)

    for(s in 1:nrow(termset_fraguniq)){
        termset_fraguniq[s,2] = nrow(as.data.frame(termset_frag[which(termset_frag[,1] == termset_fraguniq[s,1]),]))
        }

    CN = paste(termset_fraguniq[1,1],termset_fraguniq[2,1],sep=" ")

    SEPERATE_INPUT[i,1] = CN
        }
        }

    SEPERATE_INPUT_fit = SEPERATE_INPUT[,-1]


    DF_tmp = cbind("","",0)
    colnames(DF_tmp) = c("individual","group","value")
    DF_tmp = DF_tmp[-1,]
    DF_tmp = as.data.frame(DF_tmp)
    DF_tmp[,3] = as.numeric(DF_tmp[,3])

    for(i in 1:nrow(SEPERATE_INPUT)){
        for(s in 1:ncol(SEPERATE_INPUT_fit)){
            DF_tmp_frag = cbind(LABEL[s],SEPERATE_INPUT[i,1],SEPERATE_INPUT_fit[i,s])
            DF_tmp_frag = as.data.frame(DF_tmp_frag)
            DF_tmp_frag[,3] = as.numeric(DF_tmp_frag[,3])
            colnames(DF_tmp_frag) = c("individual","group","value")

            DF_tmp = rbind(DF_tmp,DF_tmp_frag)
            }

        DF_tmp_frag = cbind("",SEPERATE_INPUT[i,1],0)
        colnames(DF_tmp_frag) = c("individual","group","value")
        DF_tmp_frag = as.data.frame(DF_tmp_frag)
        DF_tmp_frag[,3] = as.numeric(DF_tmp_frag[,3])

        DF_tmp = rbind(DF_tmp,DF_tmp_frag)
        DF_tmp = rbind(DF_tmp,DF_tmp_frag)

        }

    DFcir = DF_tmp
    DFcir$id = c(1:nrow(DFcir))

termlist = DFcir[,2]
termlist = as.data.frame(termlist)
termlist = unique(termlist)
colnames(termlist) = "group"
termlist$set = c(1:nrow(termlist))

SEPERATE_INPUT_set = SEPERATE_INPUT[,1]
SEPERATE_INPUT_set = as.data.frame(SEPERATE_INPUT_set)
SEPERATE_INPUT_set = unique(SEPERATE_INPUT_set)

for(q in 1:nrow(SEPERATE_INPUT_set)){
tt = as.data.frame(SEPERATE_INPUT[which(SEPERATE_INPUT[,1] == SEPERATE_INPUT_set[q,1]),])

if(nrow(tt)>1){
    
    for(w in 1:nrow(tt)){
        tt[w,1] = paste0(tt[w,1]," tree",w)
        }

    SEPERATE_INPUT = as.data.frame(SEPERATE_INPUT[which(SEPERATE_INPUT[,1] != SEPERATE_INPUT_set[q,1]),])
    SEPERATE_INPUT = rbind(SEPERATE_INPUT,tt)
    }
}

SEPERATE_INPUT = SEPERATE_INPUT[order(-SEPERATE_INPUT[,(condition[1]+1)]),]

label_data <- DFcir
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar  
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

empty_bar = 2
base_data <- DFcir %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))
 
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

base_data = inner_join(base_data,termlist,by="group")

base_data$case = ""

LABEL_DF = as.data.frame(LABEL)

LABEL_DF2 = LABEL_DF[condition,]
LABEL_DF2 = as.data.frame(LABEL_DF2)

NAME = LABEL_DF2[1,1]
for(i in 2:nrow(LABEL_DF2)){
    NAME = paste0(NAME,"\t",LABEL_DF2[i,1])
    }

base_data[1,ncol(base_data)] = NAME

MAX = ceiling(ceiling(max(DFcir[,3]))/10)*10
SET = MAX/2
MAX = SET + MAX

terms = DFcir[,2]
terms = as.data.frame(terms)
terms = unique(terms)

condition_DF = as.data.frame(condition)

for(i in 1:nrow(condition_DF)){
terms$num = 0
colnames(terms)[ncol(terms)] = paste0("RF_",condition_DF[i,1])
    }


    MANGO_SEPERATE_INPUT = MANGO_SEPERATE_forMULTI_range_path
    SEPERATE_INPUT = utils::read.table(MANGO_SEPERATE_INPUT, sep="\t", header=TRUE, fill=TRUE)

    condition = SEPERATE_INPUT[1,ncol(SEPERATE_INPUT)]
    condition = eval(parse(text = condition))

    SEPERATE_INPUT = SEPERATE_INPUT[,-ncol(SEPERATE_INPUT)]
    SEPERATE_INPUT = stats::na.omit(SEPERATE_INPUT)
    SEPERATE_INPUT = SEPERATE_INPUT[order(-SEPERATE_INPUT[,(condition[1]+1)]),]

for(i in 1:nrow(LABEL_DF)){
    GO_RAW = utils::read.delim(GOTABLE[i,1], sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
    GO_DF  = GO_RAW[,c(2,5)]
    colnames(GO_DF)[2] = paste0("RF_",i)
    GO_DF[is.na(GO_DF)] = 0

    for(z in 1:nrow(terms)){
    SEPERATE_INPUT_DF = as.data.frame(strsplit(SEPERATE_INPUT[z,1], split = "//"))
    colnames(SEPERATE_INPUT_DF) = "Description"

    for(s in 1:nrow(SEPERATE_INPUT_DF)){
            SEPERATE_INPUT_DF[s,1] = gsub("_", " ", SEPERATE_INPUT_DF[s,1])
        }

    DATA_RAW_DF = inner_join(SEPERATE_INPUT_DF,GO_DF,by = "Description")
    terms[z,(i+1)] = nrow(DATA_RAW_DF)
        }
    }

colnames(terms)[1] = "group"

LABEL_DF = as.data.frame(LABEL)
full_condition = c(1:nrow(LABEL_DF))
pos_condition = condition


full_condition_DF = as.data.frame(full_condition)
colnames(full_condition_DF) = "condition"
pos_condition_DF = as.data.frame(pos_condition)
colnames(pos_condition_DF) = "condition"
nega_condition_DF = anti_join(full_condition_DF,pos_condition_DF,by = "condition")
nega_condition = nega_condition_DF[,1]

for(i in 1:nrow(nega_condition_DF)){
    coln = nega_condition[i]
    terms[,coln+1] = 0
    }

DFcir2 = cbind(LABEL[1],terms[1,1],terms[1,2])
DFcir2 = as.data.frame(DFcir2)
colnames(DFcir2) = c("individual","group","num")
DFcir2 = DFcir2[-1,]

for(i in 1:nrow(SEPERATE_INPUT)){
    for(s in 1:nrow(LABEL_DF)){
        DFcir2_frag = cbind(LABEL[s],terms[i,1],terms[i,s+1])
        DFcir2_frag = as.data.frame(DFcir2_frag)
        colnames(DFcir2_frag) = c("individual","group","num")

        DFcir2 = rbind(DFcir2,DFcir2_frag)
        }
    }

DFcir = full_join(DFcir,DFcir2,by=c("individual","group"))
DFcir[is.na(DFcir)] = 0

DEG_list_DF = utils::read.delim(DEG_list_pathDF[1,1], sep = "\t", header = FALSE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)

for(i in 1:nrow(LABEL_DF)){
    DEG = utils::read.delim(DEG_list_DF[i,1], sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
    DEG = DEG[,c(1,3)]

    terms[,(i+1)] = terms[,(i+1)]/nrow(DEG)*100

    }

terms_geneS = terms

    MANGO_SEPERATE_INPUT = MANGO_SEPERATE_forMULTI_range_path
    SEPERATE_INPUT = utils::read.table(MANGO_SEPERATE_INPUT, sep="\t", header=TRUE, fill=TRUE)

    condition = SEPERATE_INPUT[1,ncol(SEPERATE_INPUT)]
    condition = eval(parse(text = condition))

    SEPERATE_INPUT = SEPERATE_INPUT[,-ncol(SEPERATE_INPUT)]
    SEPERATE_INPUT = stats::na.omit(SEPERATE_INPUT)
    SEPERATE_INPUT = SEPERATE_INPUT[order(-SEPERATE_INPUT[,(condition[1]+1)]),]


        for(z in 1:nrow(terms)){
            SEPERATE_INPUT_DF = as.data.frame(strsplit(SEPERATE_INPUT[z,1], split = "//"))
            colnames(SEPERATE_INPUT_DF) = "Description"
        
        for(s in 1:nrow(SEPERATE_INPUT_DF)){
            SEPERATE_INPUT_DF[s,1] = gsub("_", " ", SEPERATE_INPUT_DF[s,1])
            }

        for(q in 1:nrow(GOTABLE)){
            GO_RAW = utils::read.delim(GOTABLE[q,1], sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
            GO_DF  = GO_RAW[,c(2,11)]
            colnames(GO_DF)[2] = paste0("RF_",q)
            GO_DF[is.na(GO_DF)] = 0

            DEG = utils::read.delim(DEG_list_DF[q,1], sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
            DEG = DEG[,c(1,3)]
            
            DATA_RAW_DF = inner_join(SEPERATE_INPUT_DF,GO_DF,by = "Description")
            DATA_RAW_gene = DATA_RAW_DF[,2]
            DATA_RAW_gene = as.data.frame(DATA_RAW_gene)
            
            DF_gene = "gene"
            DF_gene = as.data.frame(DF_gene)
            colnames(DF_gene) = "gene"
            DF_gene = DF_gene[-1,]

        for(s in 1:nrow(SEPERATE_INPUT_DF)){
            DF_gene_frag = as.data.frame(strsplit(DATA_RAW_gene[s,1], split = "/"))
            colnames(DF_gene_frag) = "gene"
            DF_gene = rbind(DF_gene,DF_gene_frag)
            DF_gene = unique(DF_gene)
            DF_gene = stats::na.omit(DF_gene)
            }

            DF_gene = inner_join(DF_gene,DEG,by="gene")
            terms_geneS[z,(q+1)] = median(DF_gene[,2])
            }
            }

terms_geneS[is.na(terms_geneS)] = 0

for(i in 1:nrow(nega_condition_DF)){
    coln = nega_condition[i]
    terms_geneS[,coln+1] = 0
    }

terms_tmp = terms[,-1]
terms$genenum = 0

for(i in 1:nrow(terms_tmp)){
    terms[i,ncol(terms)] = sum(terms_tmp[i,])/nrow(condition_DF)
    }

terms_geneS_tmp = terms_geneS[,-1]
terms_geneS$geneS = 0

for(i in 1:nrow(terms_geneS_tmp)){
    terms_geneS[i,ncol(terms_geneS)] = sum(terms_geneS_tmp[i,])/nrow(condition_DF)
    }

terms = terms[,c(1,ncol(terms))]
terms_geneS = terms_geneS[,c(1,ncol(terms_geneS))]

terms = inner_join(terms,terms_geneS,by="group")

ratio_max = ceiling(ceiling(max(terms$genenum))/5)*5

R_max = -SET*0.7
R_min = -SET*1.3

for(i in 1:nrow(terms)){
    terms[i,2] = ((terms[i,2]/5)*(R_max-R_min))+R_min
    }

base_data = inner_join(base_data,terms,by="group")
DFcir[,5] = as.numeric(DFcir[,5])

terms_fix = terms[,1]
terms_fix = as.data.frame(terms_fix)
colnames(terms_fix) = "group"

DFcir = inner_join(DFcir,terms_fix,by = "group")


DFcir_group = DFcir[,2]
DFcir_group = as.data.frame(DFcir_group)
colnames(DFcir_group) = "group"
DFcir_group$sub = "NOT"

DFcir_group_set = unique(DFcir_group)

DFcir_group_set$mark = 0

for(i in 1:nrow(DFcir_group_set)){
    DFcir_group_set[i,3] = nrow(DFcir_group[which(DFcir_group[,1] == DFcir_group_set[i,1]),])
    }

DFcir_group_set[,3] = DFcir_group_set[,3]/(nrow(LABEL_DF)+2)

DFcir_group_dupset = DFcir_group_set[which(DFcir_group_set[,3] > 1),]

if(nrow(DFcir_group_dupset)>0){

for(i in 1:nrow(DFcir_group_dupset)){
    max = DFcir_group_dupset[i,3]
    
    for(s in 2:max){
        pass_set = c(1:(nrow(LABEL_DF)+2))
        pass_set = pass_set + ((nrow(LABEL_DF)+2)*(s-1))
        codepath = DFcir[which(DFcir[,2] == DFcir_group_dupset[i,1]),][pass_set,4]

        DFcir[codepath,2] = rep(paste0(DFcir_group_dupset[i,1]," Tree",s), times = (nrow(LABEL_DF)+2)) 
    }
}
    }

termlist = DFcir[,2]
termlist = as.data.frame(termlist)
termlist = unique(termlist)
colnames(termlist) = "group"
termlist$set = c(1:nrow(termlist))


label_data <- DFcir
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar  
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

empty_bar = 2
base_data <- DFcir %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))
 
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

base_data = inner_join(base_data,termlist,by="group")

base_data$case = ""
LABEL_DF2 = LABEL_DF[condition,]
LABEL_DF2 = as.data.frame(LABEL_DF2)

NAME = LABEL_DF2[1,1]
for(i in 2:nrow(LABEL_DF2)){
    NAME = paste0(NAME,"\t",LABEL_DF2[i,1])
    }

base_data[1,ncol(base_data)] = NAME

MAX = ceiling(ceiling(max(DFcir[,3]))/10)*10
SET = MAX/2
MAX = SET + MAX

terms = DFcir[,2]
terms = as.data.frame(terms)
terms = unique(terms)

condition_DF = as.data.frame(condition)

for(i in 1:nrow(condition_DF)){
terms$num = 0
colnames(terms)[ncol(terms)] = paste0("RF_",condition_DF[i,1])
    }

    MANGO_SEPERATE_INPUT = MANGO_SEPERATE_forMULTI_range_path
    SEPERATE_INPUT = utils::read.table(MANGO_SEPERATE_INPUT, sep="\t", header=TRUE, fill=TRUE)

    condition = SEPERATE_INPUT[1,ncol(SEPERATE_INPUT)]
    condition = eval(parse(text = condition))

    SEPERATE_INPUT = SEPERATE_INPUT[,-ncol(SEPERATE_INPUT)]
    SEPERATE_INPUT = stats::na.omit(SEPERATE_INPUT)
    SEPERATE_INPUT = SEPERATE_INPUT[order(-SEPERATE_INPUT[,(condition[1]+1)]),]

LABEL_DF = as.data.frame(LABEL)

for(i in 1:nrow(LABEL_DF)){
    GO_RAW = utils::read.delim(GOTABLE[i,1], sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
    GO_DF  = GO_RAW[,c(2,5)]
    colnames(GO_DF)[2] = paste0("RF_",i)
    GO_DF[is.na(GO_DF)] = 0

    for(z in 1:nrow(terms)){
    SEPERATE_INPUT_DF = as.data.frame(strsplit(SEPERATE_INPUT[z,1], split = "//"))
    colnames(SEPERATE_INPUT_DF) = "Description"

    for(s in 1:nrow(SEPERATE_INPUT_DF)){
            SEPERATE_INPUT_DF[s,1] = gsub("_", " ", SEPERATE_INPUT_DF[s,1])
        }

    DATA_RAW_DF = inner_join(SEPERATE_INPUT_DF,GO_DF,by = "Description")
    terms[z,(i+1)] = nrow(DATA_RAW_DF)
        }
    }

colnames(terms)[1] = "group"

LABEL_DF = as.data.frame(LABEL)
full_condition = c(1:nrow(LABEL_DF))
pos_condition = condition


full_condition_DF = as.data.frame(full_condition)
colnames(full_condition_DF) = "condition"
pos_condition_DF = as.data.frame(pos_condition)
colnames(pos_condition_DF) = "condition"
nega_condition_DF = anti_join(full_condition_DF,pos_condition_DF,by = "condition")
nega_condition = nega_condition_DF[,1]

for(i in 1:nrow(nega_condition_DF)){
    coln = nega_condition[i]
    terms[,coln+1] = 0
    }

DEG_list_DF = utils::read.delim(DEG_list_pathDF[1,1], sep = "\t", header = FALSE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)

for(i in 1:nrow(LABEL_DF)){
    DEG = utils::read.delim(DEG_list_DF[i,1], sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
    DEG = DEG[,c(1,3)]

    terms[,(i+1)] = terms[,(i+1)]/nrow(DEG)*100

    }

terms_geneS = terms

    MANGO_SEPERATE_INPUT = MANGO_SEPERATE_forMULTI_range_path
    SEPERATE_INPUT = utils::read.table(MANGO_SEPERATE_INPUT, sep="\t", header=TRUE, fill=TRUE)

    condition = SEPERATE_INPUT[1,ncol(SEPERATE_INPUT)]
    condition = eval(parse(text = condition))

    SEPERATE_INPUT = SEPERATE_INPUT[,-ncol(SEPERATE_INPUT)]
    SEPERATE_INPUT = stats::na.omit(SEPERATE_INPUT)
    SEPERATE_INPUT = SEPERATE_INPUT[order(-SEPERATE_INPUT[,(condition[1]+1)]),]


        for(z in 1:nrow(terms)){
            SEPERATE_INPUT_DF = as.data.frame(strsplit(SEPERATE_INPUT[z,1], split = "//"))
            colnames(SEPERATE_INPUT_DF) = "Description"
        
        for(s in 1:nrow(SEPERATE_INPUT_DF)){
            SEPERATE_INPUT_DF[s,1] = gsub("_", " ", SEPERATE_INPUT_DF[s,1])
            }

        for(q in 1:nrow(GOTABLE)){
            GO_RAW = utils::read.delim(GOTABLE[q,1], sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
            GO_DF  = GO_RAW[,c(2,11)]
            colnames(GO_DF)[2] = paste0("RF_",q)
            GO_DF[is.na(GO_DF)] = 0

            DEG = utils::read.delim(DEG_list_DF[q,1], sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
            DEG = DEG[,c(1,3)]
            
            DATA_RAW_DF = inner_join(SEPERATE_INPUT_DF,GO_DF,by = "Description")
            DATA_RAW_gene = DATA_RAW_DF[,2]
            DATA_RAW_gene = as.data.frame(DATA_RAW_gene)
            
            DF_gene = "gene"
            DF_gene = as.data.frame(DF_gene)
            colnames(DF_gene) = "gene"
            DF_gene = DF_gene[-1,]

        for(s in 1:nrow(SEPERATE_INPUT_DF)){
            DF_gene_frag = as.data.frame(strsplit(DATA_RAW_gene[s,1], split = "/"))
            colnames(DF_gene_frag) = "gene"
            DF_gene = rbind(DF_gene,DF_gene_frag)
            DF_gene = unique(DF_gene)
            DF_gene = stats::na.omit(DF_gene)
            }

            DF_gene = inner_join(DF_gene,DEG,by="gene")
            terms_geneS[z,(q+1)] = median(DF_gene[,2])
            }
            }

terms_geneS[is.na(terms_geneS)] = 0

for(i in 1:nrow(nega_condition_DF)){
    coln = nega_condition[i]
    terms_geneS[,coln+1] = 0
    }

terms_tmp = terms[,-1]
terms$genenum = 0

for(i in 1:nrow(terms_tmp)){
    terms[i,ncol(terms)] = sum(terms_tmp[i,])/nrow(condition_DF)
    }

terms_geneS_tmp = terms_geneS[,-1]
terms_geneS$geneS = 0

for(i in 1:nrow(terms_geneS_tmp)){
    terms_geneS[i,ncol(terms_geneS)] = sum(terms_geneS_tmp[i,])/nrow(condition_DF)
    }

terms = terms[,c(1,ncol(terms))]
terms_geneS = terms_geneS[,c(1,ncol(terms_geneS))]

terms = inner_join(terms,terms_geneS,by="group")

ratio_max = ceiling(ceiling(max(terms$genenum))/2)*2

R_max = -SET*0.7
R_min = -SET*1.3

for(i in 1:nrow(terms)){
    terms[i,2] = ((terms[i,2]/2)*(R_max-R_min))+R_min
    }

base_data = inner_join(base_data,terms,by="group")
DFcir[,5] = as.numeric(DFcir[,5])

    if(trends == "UP"){
        basecol = "#DD686CFF"
        fillcol = "pink2"
        fontcol = "darkred"
    }

    if(trends == "DOWN"){
        basecol = "#33E4FFFF"
        fillcol = "skyblue2"
        fontcol = "darkblue"
    }

    if(trends == "SIG"){
        basecol = "tan"
        fillcol = "tan2"
        fontcol = "brown"
    }

if(trends == "UP"){colorset = c("gray40", "#F5B78EFF", "#F19C7CFF", "#EA8171FF", "#DD686CFF", "#CA5268FF", "#B13F64FF")}
if(trends == "DOWN"){colorset = c("gray40", "#CCFDFFFF","#99F8FFFF","#66F0FFFF","#33E4FFFF","#00AACCFF","#007A99FF")}

plot_multi = ggplot(DFcir) +      
geom_bar(aes(x=as.factor(id), y=value, fill=num), stat="identity", alpha=1) +
geom_segment(data=grid_data, aes(x = end-1, y = SET*2, xend = 0, yend = SET*2), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
geom_segment(data=grid_data, aes(x = end-1, y = SET, xend = 0, yend = SET), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
geom_segment(data=grid_data, aes(x = end-1, y = 0, xend = 0, yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +


geom_segment(data=base_data, aes(x = title, y = 0, xend = title, yend = SET*2), linetype = "dashed", color = "gray40") + 
geom_segment(data=base_data, aes(x = title, y = -SET*1.3, xend = title, yend = -SET*0.7), linetype = "dashed", color = "gray40") + 

geom_segment(data=grid_data, aes(x = end-1, y = -SET*1.3, xend = 0, yend = -SET*1.3), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
geom_segment(data=grid_data, aes(x = end-1, y = -SET*1, xend = 0, yend = -SET*1), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
geom_segment(data=grid_data, aes(x = end-1, y = -SET*0.7, xend = 0, yend = -SET*0.7), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +



geom_bar(aes(x=as.factor(id), y=value, fill=num), stat="identity", alpha=1) +

geom_segment(data=base_data, aes(x = start, y = 0, xend = end, yend = 0), colour = "black", alpha=1, size=1 , inherit.aes = FALSE ) +
geom_segment(data=base_data, aes(x = start, y = -SET*0.7, xend = end, yend = -SET*0.7), colour = "black", alpha=1, size=1 , inherit.aes = FALSE ) +
geom_segment(data=base_data, aes(x = start, y = -SET*1.3, xend = end, yend = -SET*1.3), colour = "black", alpha=1, size=1 , inherit.aes = FALSE ) +

geom_point(data=base_data, aes(x = title, y = genenum, size = geneS), color = fontcol,show.legend = TRUE) +
scale_size_continuous(range = c(0.5, 3)) +

annotate("text", x = rep(max(DFcir$id)+1,4), y = c(0,SET,SET*2,SET*(2.2)), label = c("0",SET ,SET*2,"HWES") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +

annotate("text", x = rep(max(DFcir$id)-1.5,4), y = c(-SET*1.3,-SET*1,-SET*0.7,-SET*0.55), label = c("0%",paste0(ratio_max/2,"%"),paste0(ratio_max,"%"),"DEG ratio") , color="grey", size=2 , angle=0, fontface="bold", hjust=0) +


ylim(-MAX*(0.9),MAX*(1.05)) +
theme_minimal() +
theme(#legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()) +
coord_polar() + 
scale_fill_gradientn("Number of terms", colours = colorset) +



geom_text(data=label_data, aes(x=id, y=-SET*(0.6), label=individual, hjust=hjust), color=fillcol, fontface="bold",alpha=1, size=label_size, angle= label_data$angle, inherit.aes = FALSE ) +


geom_text(data=base_data, aes(x = title, y = MAX*(1.05), label=str_wrap(group, 5)), colour = "black", alpha=1, size=tree_size, fontface="bold", inherit.aes = FALSE) +
geom_text(data=base_data, aes(x = title, y = -SET*1.6, label=str_wrap(set, 5)), colour = basecol, alpha=1, size=5, fontface="bold", inherit.aes = FALSE) +
geom_text(data=base_data, aes(x = title, y = -MAX*(0.9), label=str_wrap(case, 5)), colour = fontcol, alpha=1, size=project_size, fontface="bold", inherit.aes = FALSE)
    }


    options(repr.plot.width = width, repr.plot.height = height, repr.plot.res = 300, repr.plot.pointsize = 10)
    plot_multi
    }

