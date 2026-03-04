#' Build the MANGO term tree list across multiple cases
#'
#' Aggregates \code{total_Dependency} strings across preprocessing outputs,
#' computes pairwise overlap ratio, and creates a tree/subtree mapping based on
#' a similarity threshold.
#'
#' @param PASSED_RATIO Threshold on \code{Ratio_of_passed_Dependency} (column 4).
#' @param PASSED_NUM Threshold on passed dependency count (column 3).
#' @param similarity Similarity cutoff (%) used to define similar trees.
#' @param fileTABLE_path Path to a table listing MANGO preprocessing files (one per row).
#' @param outputPATH Output path to write the final tree/subtree table.
#'
#' @return Writes \code{outputPATH}. Returns \code{NULL} invisibly.
#'
#' @examples
#' stopifnot(is.function(MANGO_TERMLISTING))
#' @export


MANGO_TERMLISTING = function(PASSED_RATIO,
                             PASSED_NUM,
                             similarity,
                             fileTABLE_path,
                             outputPATH){

    fileTABLE = utils::read.table(fileTABLE_path, sep="\t", header=T, fill=TRUE)

    MANGO_PREPROCESSING = utils::read.table(fileTABLE[1,1], sep="\t", header=T, fill=TRUE)
    MANGO_PREPROCESSING = MANGO_PREPROCESSING[which(MANGO_PREPROCESSING[,4] >= PASSED_RATIO),]
    MANGO_PREPROCESSING = MANGO_PREPROCESSING[which(MANGO_PREPROCESSING[,3] >= PASSED_NUM),]
    MANGO_PREPROCESSING = MANGO_PREPROCESSING[order(-MANGO_PREPROCESSING$Ratio_of_passed_Dependency),]

    MANGO_TREE = MANGO_PREPROCESSING[,6]
    MANGO_TREE = as.data.frame(MANGO_TREE)
    colnames(MANGO_TREE) = "total_Dependency"

    for(i in 2:nrow(fileTABLE)){
        MANGO_PREPROCESSING = utils::read.table(fileTABLE[i,1], sep="\t", header=T, fill=TRUE)
        MANGO_PREPROCESSING = MANGO_PREPROCESSING[which(MANGO_PREPROCESSING[,4] >= PASSED_RATIO),]
        MANGO_PREPROCESSING = MANGO_PREPROCESSING[which(MANGO_PREPROCESSING[,3] >= PASSED_NUM),] 
        MANGO_PREPROCESSING = MANGO_PREPROCESSING[order(-MANGO_PREPROCESSING$Ratio_of_passed_Dependency),]
    
        MANGO_TREE_frag = MANGO_PREPROCESSING[,6]
        MANGO_TREE_frag = as.data.frame(MANGO_TREE_frag)
        colnames(MANGO_TREE_frag) = "total_Dependency"

        MANGO_TREE = rbind(MANGO_TREE,MANGO_TREE_frag)
    }

    MANGO_TREE = unique(MANGO_TREE)

    MANGO_TREE_BB = MANGO_TREE

    MANGO_TREE_COMPARE = cbind(MANGO_TREE_BB[1,1],MANGO_TREE_BB)
    colnames(MANGO_TREE_COMPARE) = c("T1","T2")

    for(i in 2:nrow(MANGO_TREE_BB)){
        MANGO_TREE_COMPARE_frag = cbind(MANGO_TREE_BB[i,1],MANGO_TREE_BB)
        colnames(MANGO_TREE_COMPARE_frag) = c("T1","T2")
        MANGO_TREE_COMPARE = rbind(MANGO_TREE_COMPARE,MANGO_TREE_COMPARE_frag)
        }

    MANGO_TREE_COMPARE$RATIO = 0

    for(i in 1:nrow(MANGO_TREE_COMPARE)){
        T1 = as.data.frame(strsplit(MANGO_TREE_COMPARE[i,1], split = "//"))
        T2 = as.data.frame(strsplit(MANGO_TREE_COMPARE[i,2], split = "//"))
    
        colnames(T1) = "Term"
        colnames(T2) = "Term"
    
        inter = nrow(dplyr::inner_join(T1,T2,by = "Term"))
        nT1 = nrow(T1)
        nT2 = nrow(T2)
    
        if(nT1>=nT2){mom = nT1}
        if(nT1<nT2){mom = nT2}
    
    MANGO_TREE_COMPARE[i,3] = round((inter/mom)*100,3)
    }

    MANGO_TREE_COMPARE_SIM = MANGO_TREE_COMPARE[which(MANGO_TREE_COMPARE[,3] >= similarity),]
    MANGO_TREE_COMPARE_DIFF = MANGO_TREE_COMPARE[which(MANGO_TREE_COMPARE[,3] < similarity),]

    MANGO_TREE_COMPARE_SIM_NOT100 = MANGO_TREE_COMPARE_SIM[which(MANGO_TREE_COMPARE_SIM[,3] < 100),]

    if(nrow(MANGO_TREE_COMPARE_SIM_NOT100) == 0){
        MANGO_TREE_simlist = as.data.frame(matrix(ncol=2, nrow=0))
        colnames(MANGO_TREE_simlist) = c("tree","subtree")
    
        MANGO_TREE_simlist_frag = as.data.frame(matrix(ncol=2, nrow=0))
        colnames(MANGO_TREE_simlist_frag) = c("tree","subtree")
    } 
    else {

        i = 1
        T1 = nrow(as.data.frame(strsplit(MANGO_TREE_COMPARE_SIM_NOT100[i,1], split = "//")))
        T2 = nrow(as.data.frame(strsplit(MANGO_TREE_COMPARE_SIM_NOT100[i,2], split = "//")))

        if(T1>=T2){
            MANGO_TREE_simlist = cbind(MANGO_TREE_COMPARE_SIM_NOT100[i,1],MANGO_TREE_COMPARE_SIM_NOT100[i,2])
            MANGO_TREE_simlist = as.data.frame(MANGO_TREE_simlist)
            colnames(MANGO_TREE_simlist) = c("tree","subtree")
            }

        if(T1<T2){
            MANGO_TREE_simlist = cbind(MANGO_TREE_COMPARE_SIM_NOT100[i,2],MANGO_TREE_COMPARE_SIM_NOT100[i,1])
            MANGO_TREE_simlist = as.data.frame(MANGO_TREE_simlist)
            colnames(MANGO_TREE_simlist) = c("tree","subtree")
            }

        for(i in 2:nrow(MANGO_TREE_COMPARE_SIM_NOT100)){
            T1 = nrow(as.data.frame(strsplit(MANGO_TREE_COMPARE_SIM_NOT100[i,1], split = "//")))
            T2 = nrow(as.data.frame(strsplit(MANGO_TREE_COMPARE_SIM_NOT100[i,2], split = "//")))
        
            if(T1>=T2){
                MANGO_TREE_simlist_frag = cbind(MANGO_TREE_COMPARE_SIM_NOT100[i,1],MANGO_TREE_COMPARE_SIM_NOT100[i,2])
                MANGO_TREE_simlist_frag = as.data.frame(MANGO_TREE_simlist_frag)
                colnames(MANGO_TREE_simlist_frag) = c("tree","subtree")
                }

            if(T1<T2){
                MANGO_TREE_simlist_frag = cbind(MANGO_TREE_COMPARE_SIM_NOT100[i,2],MANGO_TREE_COMPARE_SIM_NOT100[i,1])
                MANGO_TREE_simlist_frag = as.data.frame(MANGO_TREE_simlist_frag)
                colnames(MANGO_TREE_simlist_frag) = c("tree","subtree")
                }

            MANGO_TREE_simlist = rbind(MANGO_TREE_simlist,MANGO_TREE_simlist_frag)
            MANGO_TREE_simlist = unique(MANGO_TREE_simlist)
            }
        } 

    MANGO_TREE_simlist_tree_1 = MANGO_TREE_simlist[,1]
    MANGO_TREE_simlist_tree_1 = as.data.frame(MANGO_TREE_simlist_tree_1)
    MANGO_TREE_simlist_tree_1 = unique(MANGO_TREE_simlist_tree_1)
    MANGO_TREE_simlist_tree_1 = as.data.frame(MANGO_TREE_simlist_tree_1)
    colnames(MANGO_TREE_simlist_tree_1) = "tree"

    MANGO_TREE_simlist_tree_2 = MANGO_TREE_simlist[,2]
    MANGO_TREE_simlist_tree_2 = as.data.frame(MANGO_TREE_simlist_tree_2)
    MANGO_TREE_simlist_tree_2 = unique(MANGO_TREE_simlist_tree_2)
    MANGO_TREE_simlist_tree_2 = as.data.frame(MANGO_TREE_simlist_tree_2)
    colnames(MANGO_TREE_simlist_tree_2) = "tree"

    MANGO_TREE_simlist_1 = rbind(MANGO_TREE_simlist_tree_1,MANGO_TREE_simlist_tree_2)
    MANGO_TREE_simlist_1 = unique(MANGO_TREE_simlist_1)
    MANGO_TREE_simlist_1 = as.data.frame(MANGO_TREE_simlist_1)
    colnames(MANGO_TREE_simlist_1) = "tree"

    MANGO_TREE_COMPARE_SIM_100 = MANGO_TREE_COMPARE_SIM[which(MANGO_TREE_COMPARE_SIM[,3] == 100),]

    MANGO_TREE_simlist_2 = as.data.frame(MANGO_TREE_COMPARE_SIM_100[,1])
    MANGO_TREE_simlist_2 = unique(MANGO_TREE_simlist_2)
    MANGO_TREE_simlist_2 = as.data.frame(MANGO_TREE_simlist_2)
    colnames(MANGO_TREE_simlist_2) = "tree"

    MANGO_TREE_simlist_tree = rbind(MANGO_TREE_simlist_tree_1,MANGO_TREE_simlist_tree_2)

    MANGO_TREE_simlist_tree = unique(MANGO_TREE_simlist_tree)
    MANGO_TREE_simlist_tree = as.data.frame(MANGO_TREE_simlist_tree)
    colnames(MANGO_TREE_simlist_tree) = "tree"

    MANGO_TREE_simlist_tree$tree <- as.character(MANGO_TREE_simlist_tree$tree)

    MANGO_TREE_simlist_2 = as.data.frame(MANGO_TREE_COMPARE_SIM_100[,1])
    MANGO_TREE_simlist_2 = unique(MANGO_TREE_simlist_2)
    MANGO_TREE_simlist_2 = as.data.frame(MANGO_TREE_simlist_2)
    colnames(MANGO_TREE_simlist_2) = "tree"
    MANGO_TREE_simlist_2$subtree = "NOT"

    MANGO_TREE_simlist = rbind(MANGO_TREE_simlist,MANGO_TREE_simlist_2)
    MANGO_TREE_simlist = unique(MANGO_TREE_simlist)

    MANGO_TREE_simlist_REST = as.data.frame(MANGO_TREE_COMPARE_DIFF[,1])
    MANGO_TREE_simlist_REST = unique(MANGO_TREE_simlist_REST)
    MANGO_TREE_simlist_REST = as.data.frame(MANGO_TREE_simlist_REST)
    colnames(MANGO_TREE_simlist_REST) = "tree"
    MANGO_TREE_simlist_REST = dplyr::anti_join(MANGO_TREE_simlist_REST,MANGO_TREE_simlist_tree,by = "tree")
    MANGO_TREE_simlist_REST$subtree = "NOT"

    MANGO_TREE = rbind(MANGO_TREE_simlist,MANGO_TREE_simlist_REST)

    utils::write.table(MANGO_TREE,
                file = outputPATH, col.names=T, row.names=F, quote=F,sep="\t")
    }
