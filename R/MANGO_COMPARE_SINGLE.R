#' Compute MANGO tree scores for a single case
#'
#' Single-condition version of \code{MANGO_COMPARE}. Produces a one-column
#' comparison table for the single case.
#'
#' @param input_TERMLISTING Path to the tree/subtree table produced by term listing.
#' @param filepath Base directory containing \code{DEG_list_name}.
#' @param fileTABLE_path Path to a table listing a single MANGO preprocessing file.
#' @param DEG_list_name DEG list filename under \code{filepath}.
#' @param GOTABLE_path Path to the GO list table written by preprocessing.
#' @param PASSED_NUM Threshold on passed dependency count (column 3).
#' @param outputPATH Output path to write the comparison table.
#'
#' @return Writes \code{outputPATH}. Returns \code{NULL} invisibly.
#'
#' @examples
#' stopifnot(is.function(MANGO_COMPARE_SINGLE))
#' @export


MANGO_COMPARE_SINGLE = function(input_TERMLISTING,
                                filepath,
                                fileTABLE_path,
                                DEG_list_name,
                                GOTABLE_path,
                                PASSED_NUM,
                                outputPATH){
    
    DEG_list_path = paste0(filepath,"/",DEG_list_name)
    DEG_list_pathDF = as.data.frame(DEG_list_path)

    MANGO_TERMLISTING = utils::read.table(input_TERMLISTING, sep="\t", header=T, fill=TRUE)
    MANGO_TERMLISTING = unique(MANGO_TERMLISTING)

    fileTABLE = utils::read.table(fileTABLE_path, sep="\t", header=T, fill=TRUE)
    GOTABLE = utils::read.table(GOTABLE_path, sep="\t", header=T, fill=TRUE)

    i = 1

    GO_RAW = utils::read.delim(GOTABLE[i,1], sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
    GO_DF = GO_RAW[,c(2,5)]
    colnames(GO_DF)[2] = paste0("RF_",i)

    GO_DF[is.na(GO_DF)] = 0

    MANGO_TERM = MANGO_TERMLISTING[,1]
    MANGO_TERM = as.data.frame(MANGO_TERM)
    MANGO_TERM = unique(MANGO_TERM)
    MANGO_TERM = as.data.frame(MANGO_TERM)
    colnames(MANGO_TERM) = "MANGO_TERM"

    MANGO_TERM_EX = MANGO_TERMLISTING[,2]
    MANGO_TERM_EX = as.data.frame(MANGO_TERM_EX)
    MANGO_TERM_EX = unique(MANGO_TERM_EX)
    MANGO_TERM_EX = as.data.frame(MANGO_TERM_EX)
    colnames(MANGO_TERM_EX) = "MANGO_TERM"

    MANGO_TERM = dplyr::anti_join(MANGO_TERM,MANGO_TERM_EX,by="MANGO_TERM")
    MANGO_TERM = unique(MANGO_TERM)

    MANGO_TERM = cbind(MANGO_TERM,0)
    colnames(MANGO_TERM)[2] = paste0("RF_",1)

    MANGO_TERM$label = "PASS"

    for(i in 1:nrow(MANGO_TERM)){
        
        TREEset = as.data.frame(strsplit(MANGO_TERM[i,1], split = "//"))
        colnames(TREEset) = "Description"

    for(s in 1:nrow(TREEset)){
            TREEset[s,1] = gsub("_", " ", TREEset[s,1])
        }

    ratio_mom = nrow(TREEset)

    TREEset = dplyr::left_join(TREEset,GO_DF,by="Description")
    TREEset[is.na(TREEset)] = 0

    TREEset_pos = TREEset[which(TREEset[,2] > 0),]

    if(nrow(TREEset_pos)<PASSED_NUM){MANGO_TERM[i,ncol(MANGO_TERM)] = "FALSE"}

    ratio_son = nrow(TREEset)
    ratio_hit = ratio_son/ratio_mom

        e = 1

    GO_RAW = utils::read.delim(GOTABLE[e,1], sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
    DEG_list = utils::read.table(DEG_list_pathDF[1,1], sep="\t", header=F, fill=TRUE)
    DEG_RAW = utils::read.table(DEG_list[e,1], sep="\t", header=T, fill=TRUE)

    for(w in 1:nrow(TREEset)){

    Description_genelist = as.data.frame(strsplit(GO_RAW[which(GO_RAW$Description == TREEset[w,1]),11], split = "/"))

    if(nrow(Description_genelist)>0){
        colnames(Description_genelist) = "gene"
        Description_genelist = dplyr::inner_join(Description_genelist,DEG_RAW,by="gene")
        TREEset[w,(e+1)] = 1/(log10(TREEset[w,(e+1)])*(-1))*median(abs(Description_genelist[,3]))
        }
        }
            MANGO_TERM[i,2] = mean(TREEset[[2]][TREEset[[2]] != 0])
            MANGO_TERM[i,2] = MANGO_TERM[i,2]*ratio_hit
                                 }

    MANGO_TERM[is.na(MANGO_TERM)] = 0
    MANGO_TERM = MANGO_TERM[which(MANGO_TERM[,3] == "PASS"),1:2]

    utils::write.table(MANGO_TERM,
                file = outputPATH, col.names=T, row.names=F, quote=F,sep="\t")
    }
