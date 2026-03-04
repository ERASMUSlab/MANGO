#' MANGO multi-case term bar plot (term-level)
#'
#' Draw a bar plot of GO terms within a selected MANGO tree for multi-case results.
#' The function supports three modes via \code{status}: \code{"common"}, \code{"specific"}, and \code{"DA"}.
#' Bars represent the term score (tscore-like), dashed segments/points encode gene distribution ratio
#' and median fold-change of genes in each term (as implemented in the original logic).
#'
#' @param filepath Directory path containing input files.
#' @param DEG_list_name DEG list table filename under \code{filepath}.
#' @param status One of \code{"common"}, \code{"specific"}, \code{"DA"}.
#' @param trends One of \code{"UP"}, \code{"DOWN"}, \code{"SIG"}.
#' @param LABEL Character vector of case labels (used for titles).
#' @param condition Condition index (used to pick the case/contrast).
#' @param clusternum Tree(cluster) index to display.
#' @param width Plot width (used by repr.plot.width option in notebooks).
#' @param height Plot height (used by repr.plot.height option in notebooks).
#'
#' @return A ggplot object.
#' @export
#'
#' @importFrom dplyr inner_join
#' @importFrom ggplot2 ggplot aes geom_bar geom_segment geom_point scale_fill_gradientn
#'   scale_size_continuous labs ylab xlab facet_grid theme element_text element_rect
#' @importFrom forcats fct_reorder
#' @importFrom utils read.table read.delim
#' @importFrom stats na.omit

MANGO_TERM_barPLOT_forMULTI = function(filepath,
                                       DEG_list_name,
                                       status,
                                       trends,
                                       LABEL,
                                       condition,
                                       clusternum,
                                       width,
                                       height){

    if(status == "common"){
        
        LABEL_DF = as.data.frame(LABEL)

blacklist = c("regulation","positive","negative","a","the","of","gene","type","I","II","to","by")
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

    fileTABLE = utils::read.table(fileTABLE_path, sep="\t", header=T, fill=TRUE)

    MANGO_SEPERATE_INPUT = MANGO_SEPERATE_path
    SEPERATE_INPUT = utils::read.table(MANGO_SEPERATE_INPUT, sep="\t", header=T, fill=TRUE)
    SEPERATE_INPUT = SEPERATE_INPUT[order(-SEPERATE_INPUT[,2]),]
    SEPERATE_INPUT = stats::na.omit(SEPERATE_INPUT)
    SEPERATE_INPUT = SEPERATE_INPUT[which(SEPERATE_INPUT[,ncol(SEPERATE_INPUT)] == (nrow(LABEL_DF) +1)),]

    fileTABLE = utils::read.table(fileTABLE_path, sep="\t", header=T, fill=TRUE)
    GOTABLE  = utils::read.table(GOTABLE_path,  sep="\t", header=T, fill=TRUE)

    i = condition
    GO_RAW = utils::read.delim(GOTABLE[i,1], sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
    GO_DF  = GO_RAW[,c(2,5,11)]
    colnames(GO_DF)[2] = paste0("RF_",i)

    GO_DF[is.na(GO_DF)] = 0

    i = 1
    SEPERATE_INPUT_DF = as.data.frame(strsplit(SEPERATE_INPUT[i,1], split = "//"))
    colnames(SEPERATE_INPUT_DF) = "Description"

    for(s in 1:nrow(SEPERATE_INPUT_DF)){
            SEPERATE_INPUT_DF[s,1] = gsub("_", " ", SEPERATE_INPUT_DF[s,1])
        }

    DATA_RAW_DF = inner_join(SEPERATE_INPUT_DF,GO_DF,by = "Description") 
    DATA_RAW_DF = DATA_RAW_DF[order(-DATA_RAW_DF[,2]),]
    DATA_RAW_DF = cbind(DATA_RAW_DF,i)

    colnames(DATA_RAW_DF) = c("Description","RF","gene","cluster")

    for(i in 2:nrow(SEPERATE_INPUT)){
    SEPERATE_INPUT_DF = as.data.frame(strsplit(SEPERATE_INPUT[i,1], split = "//"))
    colnames(SEPERATE_INPUT_DF) = "Description"

    for(s in 1:nrow(SEPERATE_INPUT_DF)){
            SEPERATE_INPUT_DF[s,1] = gsub("_", " ", SEPERATE_INPUT_DF[s,1])
        }

    DATA_RAW_DF_frag = inner_join(SEPERATE_INPUT_DF,GO_DF,by = "Description")        
    DATA_RAW_DF_frag = DATA_RAW_DF_frag[order(-DATA_RAW_DF_frag[,2]),]
    DATA_RAW_DF_frag = cbind(DATA_RAW_DF_frag,i)

    colnames(DATA_RAW_DF_frag) = c("Description","RF","gene","cluster")

    DATA_RAW_DF = rbind(DATA_RAW_DF,DATA_RAW_DF_frag)
    }

    DATA_RAW_DF[,5] = 0

    colnames(DATA_RAW_DF)[5] = "score"



for(i in 1:nrow(DATA_RAW_DF)){
gene = as.data.frame(strsplit(DATA_RAW_DF[i,3], split = "/"))
colnames(gene) = "gene"

DEG = utils::read.table(DEG_list_pathDF[1,1], sep="\t", header=F, fill=TRUE)
DEG = utils::read.table(DEG[1,1], sep="\t", header=T, fill=TRUE)
DEG[,3] = 2^abs(DEG[,3])
DEG = inner_join(DEG,gene,by="gene")

DATA_RAW_DF[i,5] = 1/(log10(DATA_RAW_DF[i,2])*(-1))*median(DEG[,3])
    }


DATA_RAW_DF = DATA_RAW_DF[which(DATA_RAW_DF[,4] == clusternum),]
DATA_RAW_DF = DATA_RAW_DF[order(DATA_RAW_DF[,5]),]

DATA_RAW_DF$generatio = 0

DEG = utils::read.table(DEG_list_pathDF[1,1], sep="\t", header=F, fill=TRUE)
DEG = utils::read.table(DEG[condition,1], sep="\t", header=T, fill=TRUE)

for(i in 1:nrow(DATA_RAW_DF)){
    DATA_RAW_DF[i,6] = (nrow(as.data.frame(strsplit(DATA_RAW_DF[i,3], split = "/")))/nrow(DEG)*10)
    }

DEG = utils::read.table(DEG_list_pathDF[1,1], sep="\t", header=F, fill=TRUE)
DEG = utils::read.table(DEG[condition,1], sep="\t", header=T, fill=TRUE)
DEG[,3] = 2^abs(DEG[,3])
DATA_RAW_DF$geneS = 0

for(i in 1:nrow(DATA_RAW_DF)){
geneSset = as.data.frame(strsplit(DATA_RAW_DF[i,3], split = "/"))
colnames(geneSset) = "gene"
geneSset = inner_join(geneSset,DEG,by="gene")
DATA_RAW_DF[i,7] = median(geneSset[,3])
    }

DATA_DF = DATA_RAW_DF[,c(1,5,6,7)]
DATA_DF$cluster = clusternum

GO_RAW = utils::read.delim(GOTABLE[condition,1], sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
    GO_DF  = GO_RAW[,c(2,9)]

DATA_DF = inner_join(DATA_DF,GO_DF,by="Description")
DATA_DF$P = log10(DATA_DF[,6])*(-1)



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

if(trends == "UP"){colorset = c( "#F5B78EFF", "#F19C7CFF", "#EA8171FF", "#DD686CFF", "#CA5268FF", "#B13F64FF")}
if(trends == "DOWN"){colorset = c( "#CCFDFFFF","#99F8FFFF","#66F0FFFF","#33E4FFFF","#00AACCFF","#007A99FF")}

plot_term = ggplot(DATA_DF, aes(score, fct_reorder(Description, score), fill=P)) +
    scale_fill_gradientn("Normalized\nAdjust P-value",colours = colorset) +
    geom_bar(stat="identity",
             position=position_dodge(width = 0.8),
             colour="black",
             width=0.85,
             linewidth=0.1)+
    geom_segment(aes(x = 0, y = Description, xend = generatio, yend = Description), linetype = "dashed", color = "gray12") + 
    geom_point(aes(x = generatio, y = Description, size = geneS), color = "gray12",show.legend = TRUE) +
    scale_size_continuous("FC of genes",range = c(1, 5)) +
    theme_dose(12) +
    labs(title=paste0(LABEL[condition]," result"),
         subtitle = paste0("Tree ",clusternum))+
    ylab("\n\n\n") +
    xlab("Tscore & Gene distribution ratio") +
    facet_grid(cluster ~ ., scales = "free_y", space = "free_y") +
    theme(strip.text.y = element_text(face="bold", size=14, color = "darkblue"),
          strip.background = element_rect(fill = fillcol, color = fontcol, linewidth = 0.6)) +
    theme(plot.title=element_text(face="bold",hjust=1,size=25),
          plot.subtitle=element_text(face="bold",hjust=1,size=20),
          axis.text.x=element_text(face="bold",size=7),
          axis.text.y=element_text(face="bold",size=15),
          axis.title.x = element_text(face="bold",size = 15),
          legend.title=element_text(face="bold",size=10), 
          legend.text=element_text(face="bold",size=8))

    }

    if(status == "specific"){

        blacklist = c("regulation","positive","negative","a","the","of","gene","type","I","II","to","by")
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

    fileTABLE = utils::read.table(fileTABLE_path, sep="\t", header=T, fill=TRUE)

    MANGO_SEPERATE_INPUT = MANGO_SEPERATE_path
    SEPERATE_INPUT = utils::read.table(MANGO_SEPERATE_INPUT, sep="\t", header=T, fill=TRUE)
    SEPERATE_INPUT = SEPERATE_INPUT[order(-SEPERATE_INPUT[,(condition+1)]),]
    SEPERATE_INPUT = stats::na.omit(SEPERATE_INPUT)
    SEPERATE_INPUT = SEPERATE_INPUT[which(SEPERATE_INPUT[,ncol(SEPERATE_INPUT)] == condition),]


    fileTABLE = utils::read.table(fileTABLE_path, sep="\t", header=T, fill=TRUE)
    GOTABLE  = utils::read.table(GOTABLE_path,  sep="\t", header=T, fill=TRUE)

    i = condition
    GO_RAW = utils::read.delim(GOTABLE[i,1], sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
    GO_DF  = GO_RAW[,c(2,5,11)]
    colnames(GO_DF)[2] = paste0("RF_",i)

    GO_DF[is.na(GO_DF)] = 0

     i = 1
    SEPERATE_INPUT_DF = as.data.frame(strsplit(SEPERATE_INPUT[i,1], split = "//"))
    colnames(SEPERATE_INPUT_DF) = "Description"

    for(s in 1:nrow(SEPERATE_INPUT_DF)){
            SEPERATE_INPUT_DF[s,1] = gsub("_", " ", SEPERATE_INPUT_DF[s,1])
        }

    DATA_RAW_DF = inner_join(SEPERATE_INPUT_DF,GO_DF,by = "Description") 
    DATA_RAW_DF = DATA_RAW_DF[order(-DATA_RAW_DF[,2]),]
    DATA_RAW_DF = cbind(DATA_RAW_DF,i)

    colnames(DATA_RAW_DF) = c("Description","RF","gene","cluster")

    for(i in 2:nrow(SEPERATE_INPUT)){
    SEPERATE_INPUT_DF = as.data.frame(strsplit(SEPERATE_INPUT[i,1], split = "//"))
    colnames(SEPERATE_INPUT_DF) = "Description"

    for(s in 1:nrow(SEPERATE_INPUT_DF)){
            SEPERATE_INPUT_DF[s,1] = gsub("_", " ", SEPERATE_INPUT_DF[s,1])
        }

    DATA_RAW_DF_frag = inner_join(SEPERATE_INPUT_DF,GO_DF,by = "Description")        
    DATA_RAW_DF_frag = DATA_RAW_DF_frag[order(-DATA_RAW_DF_frag[,2]),]
    DATA_RAW_DF_frag = cbind(DATA_RAW_DF_frag,i)

    colnames(DATA_RAW_DF_frag) = c("Description","RF","gene","cluster")

    DATA_RAW_DF = rbind(DATA_RAW_DF,DATA_RAW_DF_frag)
    }

    DATA_RAW_DF[,5] = 0
    colnames(DATA_RAW_DF)[5] = "score"

    for(i in 1:nrow(DATA_RAW_DF)){
        gene = as.data.frame(strsplit(DATA_RAW_DF[i,3], split = "/"))
        colnames(gene) = "gene"
        
        DEG = utils::read.table(DEG_list_pathDF[1,1], sep="\t", header=F, fill=TRUE)
        DEG = utils::read.table(DEG[condition,1], sep="\t", header=T, fill=TRUE)
        DEG[,3] = 2^abs(DEG[,3])
        DEG = inner_join(DEG,gene,by="gene")
        
        DATA_RAW_DF[i,5] = 1/(log10(DATA_RAW_DF[i,2])*(-1))*median(DEG[,3])
    }


DATA_RAW_DF = DATA_RAW_DF[which(DATA_RAW_DF[,4] == clusternum),]
DATA_RAW_DF = DATA_RAW_DF[order(DATA_RAW_DF[,5]),]

DATA_RAW_DF$generatio = 0

DEG = utils::read.table(DEG_list_pathDF[1,1], sep="\t", header=F, fill=TRUE)
DEG = utils::read.table(DEG[condition,1], sep="\t", header=T, fill=TRUE)

for(i in 1:nrow(DATA_RAW_DF)){
    DATA_RAW_DF[i,6] = (nrow(as.data.frame(strsplit(DATA_RAW_DF[i,3], split = "/")))/nrow(DEG)*10)
    }

DEG = utils::read.table(DEG_list_pathDF[1,1], sep="\t", header=F, fill=TRUE)
DEG = utils::read.table(DEG[condition,1], sep="\t", header=T, fill=TRUE)
DEG[,3] = 2^abs(DEG[,3])
DATA_RAW_DF$geneS = 0

for(i in 1:nrow(DATA_RAW_DF)){
geneSset = as.data.frame(strsplit(DATA_RAW_DF[i,3], split = "/"))
colnames(geneSset) = "gene"
geneSset = inner_join(geneSset,DEG,by="gene")
DATA_RAW_DF[i,7] = median(geneSset[,3])
    }

DATA_DF = DATA_RAW_DF[,c(1,5,6,7)]
DATA_DF$cluster = clusternum

GO_RAW = utils::read.delim(GOTABLE[condition,1], sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
    GO_DF  = GO_RAW[,c(2,9)]

DATA_DF = inner_join(DATA_DF,GO_DF,by="Description")
DATA_DF$P = log10(DATA_DF[,6])*(-1)



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

if(trends == "UP"){colorset = c( "#F5B78EFF", "#F19C7CFF", "#EA8171FF", "#DD686CFF", "#CA5268FF", "#B13F64FF")}
if(trends == "DOWN"){colorset = c( "#CCFDFFFF","#99F8FFFF","#66F0FFFF","#33E4FFFF","#00AACCFF","#007A99FF")}

plot_term = ggplot(DATA_DF, aes(score, fct_reorder(Description, score), fill=P)) +
    scale_fill_gradientn("Normalized\nAdjust P-value",colours = colorset) +
    geom_bar(stat="identity",
             position=position_dodge(width = 0.8),
             colour="black",
             width=0.85,
             linewidth=0.1)+
    labs(title=paste0(LABEL),
         subtitle = paste0("Tree ",clusternum))+
    geom_segment(aes(x = 0, y = Description, xend = generatio, yend = Description), linetype = "dashed", color = "gray12") + 
    geom_point(aes(x = generatio, y = Description, size = geneS), color = "gray12",show.legend = TRUE) +
    scale_size_continuous("FC of genes",range = c(1, 5)) +
    theme_dose(12) +
    labs(title=paste0(LABEL[condition]," result"),
         subtitle = paste0("Tree ",clusternum))+
    ylab("\n\n\n") +
    xlab("Tscore & Gene distribution ratio") +
    facet_grid(cluster ~ ., scales = "free_y", space = "free_y") +
    theme(strip.text.y = element_text(face="bold", size=14, color = "darkblue"),
          strip.background = element_rect(fill = fillcol, color = fontcol, linewidth = 0.6)) +
    theme(plot.title=element_text(face="bold",hjust=1,size=25),
          plot.subtitle=element_text(face="bold",hjust=1,size=20),
          axis.text.x=element_text(face="bold",size=7),
          axis.text.y=element_text(face="bold",size=15),
          axis.title.x = element_text(face="bold",size = 15),
          legend.title=element_text(face="bold",size=10), 
          legend.text=element_text(face="bold",size=8))
    }


    if(status == "DA"){

        blacklist = c("regulation","positive","negative","a","the","of","gene","type","I","II","to","by")
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

    fileTABLE = utils::read.table(fileTABLE_path, sep="\t", header=T, fill=TRUE)

    MANGO_SEPERATE_INPUT = MANGO_SEPERATE_forMULTI_range_path
    SEPERATE_INPUT = utils::read.table(MANGO_SEPERATE_INPUT, sep="\t", header=T, fill=TRUE)
    SEPERATE_INPUT = stats::na.omit(SEPERATE_INPUT)
    SEPERATE_INPUT = SEPERATE_INPUT[,-ncol(SEPERATE_INPUT)]


    fileTABLE = utils::read.table(fileTABLE_path, sep="\t", header=T, fill=TRUE)
    GOTABLE  = utils::read.table(GOTABLE_path,  sep="\t", header=T, fill=TRUE)

    i = condition
    GO_RAW = utils::read.delim(GOTABLE[i,1], sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
    GO_DF  = GO_RAW[,c(2,5,11)]
    colnames(GO_DF)[2] = paste0("RF_",i)

    GO_DF[is.na(GO_DF)] = 0

     i = 1
    SEPERATE_INPUT_DF = as.data.frame(strsplit(SEPERATE_INPUT[i,1], split = "//"))
    colnames(SEPERATE_INPUT_DF) = "Description"

    for(s in 1:nrow(SEPERATE_INPUT_DF)){
            SEPERATE_INPUT_DF[s,1] = gsub("_", " ", SEPERATE_INPUT_DF[s,1])
        }

    DATA_RAW_DF = inner_join(SEPERATE_INPUT_DF,GO_DF,by = "Description") 
    DATA_RAW_DF = DATA_RAW_DF[order(-DATA_RAW_DF[,2]),]
    DATA_RAW_DF = cbind(DATA_RAW_DF,i)

    colnames(DATA_RAW_DF) = c("Description","RF","gene","cluster")

    if(nrow(SEPERATE_INPUT)>1){

    for(i in 2:nrow(SEPERATE_INPUT)){
    SEPERATE_INPUT_DF = as.data.frame(strsplit(SEPERATE_INPUT[i,1], split = "//"))
    colnames(SEPERATE_INPUT_DF) = "Description"

    for(s in 1:nrow(SEPERATE_INPUT_DF)){
            SEPERATE_INPUT_DF[s,1] = gsub("_", " ", SEPERATE_INPUT_DF[s,1])
        }

    DATA_RAW_DF_frag = inner_join(SEPERATE_INPUT_DF,GO_DF,by = "Description")        
    DATA_RAW_DF_frag = DATA_RAW_DF_frag[order(-DATA_RAW_DF_frag[,2]),]
    DATA_RAW_DF_frag = cbind(DATA_RAW_DF_frag,i)

    colnames(DATA_RAW_DF_frag) = c("Description","RF","gene","cluster")

    DATA_RAW_DF = rbind(DATA_RAW_DF,DATA_RAW_DF_frag)
    }
        }

    DATA_RAW_DF[,5] = 0
    colnames(DATA_RAW_DF)[5] = "score"

    for(i in 1:nrow(DATA_RAW_DF)){
        gene = as.data.frame(strsplit(DATA_RAW_DF[i,3], split = "/"))
        colnames(gene) = "gene"
        
        DEG = utils::read.table(DEG_list_pathDF[1,1], sep="\t", header=F, fill=TRUE)
        DEG = utils::read.table(DEG[condition,1], sep="\t", header=T, fill=TRUE)
        DEG[,3] = 2^abs(DEG[,3])
        DEG = inner_join(DEG,gene,by="gene")
        
        DATA_RAW_DF[i,5] = 1/(log10(DATA_RAW_DF[i,2])*(-1))*median(DEG[,3])
    }


DATA_RAW_DF = DATA_RAW_DF[which(DATA_RAW_DF[,4] == clusternum),]
DATA_RAW_DF = DATA_RAW_DF[order(DATA_RAW_DF[,5]),]

DATA_RAW_DF$generatio = 0

DEG = utils::read.table(DEG_list_pathDF[1,1], sep="\t", header=F, fill=TRUE)
DEG = utils::read.table(DEG[condition,1], sep="\t", header=T, fill=TRUE)

for(i in 1:nrow(DATA_RAW_DF)){
    DATA_RAW_DF[i,6] = (nrow(as.data.frame(strsplit(DATA_RAW_DF[i,3], split = "/")))/nrow(DEG)*10)
    }

DEG = utils::read.table(DEG_list_pathDF[1,1], sep="\t", header=F, fill=TRUE)
DEG = utils::read.table(DEG[condition,1], sep="\t", header=T, fill=TRUE)
DEG[,3] = 2^abs(DEG[,3])
DATA_RAW_DF$geneS = 0

for(i in 1:nrow(DATA_RAW_DF)){
geneSset = as.data.frame(strsplit(DATA_RAW_DF[i,3], split = "/"))
colnames(geneSset) = "gene"
geneSset = inner_join(geneSset,DEG,by="gene")
DATA_RAW_DF[i,7] = median(geneSset[,3])
    }

DATA_DF = DATA_RAW_DF[,c(1,5,6,7)]
DATA_DF$cluster = clusternum

GO_RAW = utils::read.delim(GOTABLE[condition,1], sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
    GO_DF  = GO_RAW[,c(2,9)]

DATA_DF = inner_join(DATA_DF,GO_DF,by="Description")
DATA_DF$P = log10(DATA_DF[,6])*(-1)



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

if(trends == "UP"){colorset = c( "#F5B78EFF", "#F19C7CFF", "#EA8171FF", "#DD686CFF", "#CA5268FF", "#B13F64FF")}
if(trends == "DOWN"){colorset = c( "#CCFDFFFF","#99F8FFFF","#66F0FFFF","#33E4FFFF","#00AACCFF","#007A99FF")}

options(repr.plot.width = width, repr.plot.height = height, repr.plot.res = 500, repr.plot.pointsize = 10)

plot_term = ggplot(DATA_DF, aes(score, fct_reorder(Description, score), fill=P)) +
    scale_fill_gradientn("Normalized\nAdjust P-value",colours = colorset) +
    geom_bar(stat="identity",
             position=position_dodge(width = 0.8),
             colour="black",
             width=0.85,
             linewidth=0.1)+
    geom_segment(aes(x = 0, y = Description, xend = generatio, yend = Description), linetype = "dashed", color = "gray12") + 
    geom_point(aes(x = generatio, y = Description, size = geneS), color = "gray12",show.legend = TRUE) +
    scale_size_continuous("FC of genes",range = c(1, 5)) +
    theme_dose(12) +
    labs(title=paste0(LABEL[condition]," result"),
         subtitle = paste0("Tree ",clusternum))+
    ylab("\n\n\n") +
    xlab("Tscore & Gene distribution ratio") +
    facet_grid(cluster ~ ., scales = "free_y", space = "free_y") +
    theme(strip.text.y = element_text(face="bold", size=14, color = "darkblue"),
          strip.background = element_rect(fill = fillcol, color = fontcol, linewidth = 0.6)) +
    theme(plot.title=element_text(face="bold",hjust=1,size=25),
          plot.subtitle=element_text(face="bold",hjust=1,size=20),
          axis.text.x=element_text(face="bold",size=7),
          axis.text.y=element_text(face="bold",size=15),
          axis.title.x = element_text(face="bold",size = 15),
          legend.title=element_text(face="bold",size=10), 
          legend.text=element_text(face="bold",size=8))
    }

    options(repr.plot.width = width, repr.plot.height = height, repr.plot.res = 500, repr.plot.pointsize = 10)
plot_term
    }
