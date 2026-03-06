#' Separate MANGO terms by condition specificity (with PASS/FAIL filtering)
#'
#' This function separates MANGO terms by condition using the comparison matrix,
#' then applies additional filtering based on the number of matched GO terms
#' per tree (>= PASSED_NUM) for each condition and for non-separated terms.
#'
#' @param input_COMPARE Path to MANGO_COMPARE output file.
#' @param fileTABLE_path Path to MANGO_PREPROCESSING_list file (condition table).
#' @param GOTABLE_path Path to GO result list table (GO_list_path).
#' @param FC Fold-change ratio cutoff used for separating condition-specific terms.
#' @param PASSED_NUM Minimum number of matched GO terms required to keep a tree.
#' @param outputPATH Output path to write the separated table.
#'
#' @return Writes a tab-delimited file to \code{outputPATH}. Returns NULL invisibly.
#' @export
#'
#' @examples
#' stopifnot(is.function(MANGO_SEPERATE))
#' \donttest{
#' # MANGO_SEPERATE(
#' #   input_COMPARE = "MANGO_COMPARE.txt",
#' #   fileTABLE_path = "MANGO_PREPROCESSING_list.txt",
#' #   GOTABLE_path = "input_DEG_GO_list.txt",
#' #   FC = 2,
#' #   PASSED_NUM = 4,
#' #   outputPATH = "MANGO_SEPERATE.txt"
#' # )
#' }

MANGO_SEPERATE = function(input_COMPARE,
                          fileTABLE_path,
                          GOTABLE_path,
                          FC,
                          PASSED_NUM,
                          outputPATH){

fileTABLE = read.table(fileTABLE_path, sep="\t", header=TRUE, fill=TRUE)    
GOTABLE  = read.table(GOTABLE_path,  sep="\t", header=TRUE, fill=TRUE)

    DOWN_MANGO_COMPARE = read.table(input_COMPARE, sep="\t", header=TRUE, fill=TRUE)
    DOWN_MANGO_COMPARE = unique(DOWN_MANGO_COMPARE)

    DOWN_MANGO_SEPERATE = DOWN_MANGO_COMPARE[1,]
    DOWN_MANGO_SEPERATE$condition = 0
    DOWN_MANGO_SEPERATE = DOWN_MANGO_SEPERATE[-1,]

    for(i in 1:nrow(fileTABLE)){
        DOWN_MANGO_SEPERATE_frag = DOWN_MANGO_COMPARE[which(DOWN_MANGO_COMPARE[,(i+1)] > 0),]
        DOWN_MANGO_SEPERATE_frag_1 = DOWN_MANGO_SEPERATE_frag[,c(1,(i+1))]

        DOWN_MANGO_SEPERATE_frag = DOWN_MANGO_COMPARE[which(DOWN_MANGO_COMPARE[,(i+1)] > 0),]
        DOWN_MANGO_SEPERATE_frag = DOWN_MANGO_SEPERATE_frag[,-(i+1)]
        rownames(DOWN_MANGO_SEPERATE_frag) = DOWN_MANGO_SEPERATE_frag[,1]
        DOWN_MANGO_SEPERATE_frag$avoidV = 0
        DOWN_MANGO_SEPERATE_frag = DOWN_MANGO_SEPERATE_frag[,-1]

        DOWN_MANGO_SEPERATE_frag_2 = DOWN_MANGO_SEPERATE_frag

        if(nrow(DOWN_MANGO_SEPERATE_frag) == 0) {next}

        if(nrow(fileTABLE)==2){
            DOWN_MANGO_SEPERATE_frag_2[,1] = DOWN_MANGO_SEPERATE_frag_1[,2]/(DOWN_MANGO_SEPERATE_frag[,1]+0.0000000001)
        }

        if(nrow(fileTABLE)>2){
            for(r in 1:(nrow(fileTABLE)-1)){
                DOWN_MANGO_SEPERATE_frag_2[,r] = DOWN_MANGO_SEPERATE_frag_1[,2]/(DOWN_MANGO_SEPERATE_frag[,r]+0.0000000001)
            }
        }

        if(nrow(fileTABLE)==2){
            DOWN_MANGO_SEPERATE_frag_2 = DOWN_MANGO_SEPERATE_frag_2[which(DOWN_MANGO_SEPERATE_frag_2[,1] >= FC),]
        }

        if(nrow(fileTABLE)>2){
            for(r in 1:(nrow(fileTABLE)-1)){
                DOWN_MANGO_SEPERATE_frag_2 = DOWN_MANGO_SEPERATE_frag_2[which(DOWN_MANGO_SEPERATE_frag_2[,r] >= FC),]
            }
        }

        if(nrow(DOWN_MANGO_SEPERATE_frag_2) == 0) {next}

        DOWN_MANGO_SEPERATE_frag_2 = rownames(DOWN_MANGO_SEPERATE_frag_2)
        DOWN_MANGO_SEPERATE_frag_2 = as.data.frame(DOWN_MANGO_SEPERATE_frag_2)
        colnames(DOWN_MANGO_SEPERATE_frag_2) = "MANGO_TERM"

        DOWN_MANGO_SEPERATE_frag = dplyr::inner_join(DOWN_MANGO_SEPERATE_frag_2, DOWN_MANGO_COMPARE, by= "MANGO_TERM")
        DOWN_MANGO_SEPERATE_frag$condition = i

        DOWN_MANGO_SEPERATE = rbind(DOWN_MANGO_SEPERATE, DOWN_MANGO_SEPERATE_frag)
    }

    DOWN_MANGO_SEPERATE_term = DOWN_MANGO_SEPERATE[,1]
    DOWN_MANGO_SEPERATE_term = as.data.frame(DOWN_MANGO_SEPERATE_term)
    colnames(DOWN_MANGO_SEPERATE_term) = "MANGO_TERM"

    DOWN_MANGO_NOT_SEPERATE = dplyr::anti_join(DOWN_MANGO_COMPARE, DOWN_MANGO_SEPERATE_term, by = "MANGO_TERM")

    NOT_NUM = 1 + nrow(fileTABLE)

    if(nrow(DOWN_MANGO_NOT_SEPERATE)>0){
        DOWN_MANGO_NOT_SEPERATE$condition = NOT_NUM
    }

    if(nrow(DOWN_MANGO_SEPERATE)>0){DOWN_MANGO_SEPERATE$label = "PASS"}
    if(nrow(DOWN_MANGO_SEPERATE)==0){DOWN_MANGO_SEPERATE$label = character(0)}

z = 1
    DOWN_MANGO_SEPERATE_fin = DOWN_MANGO_SEPERATE[which(DOWN_MANGO_SEPERATE$condition == z),]

    if(nrow(DOWN_MANGO_SEPERATE_fin)>0){
        
        for(i in 1:nrow(DOWN_MANGO_SEPERATE_fin)){
            DOWN_MANGO_SEPERATE_fin_DF = as.data.frame(strsplit(DOWN_MANGO_SEPERATE_fin[i,1], split = "//"))
            colnames(DOWN_MANGO_SEPERATE_fin_DF) = "Description"
            
            for(w in 1:nrow(DOWN_MANGO_SEPERATE_fin_DF)){
                DOWN_MANGO_SEPERATE_fin_DF[w,1] = gsub("_", " ", DOWN_MANGO_SEPERATE_fin_DF[w,1])
            }

            GO_RAW = read.delim(GOTABLE[z,1], sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
            GO_DF  = GO_RAW[,c(2,5,11)]
            colnames(GO_DF)[2] = paste0("RF_", i)
            GO_DF[is.na(GO_DF)] = 0
            DOWN_MANGO_SEPERATE_fin_DF_TF = dplyr::inner_join(DOWN_MANGO_SEPERATE_fin_DF, GO_DF, by="Description")

            if(nrow(DOWN_MANGO_SEPERATE_fin_DF_TF) < PASSED_NUM){
                DOWN_MANGO_SEPERATE_fin[i, ncol(DOWN_MANGO_SEPERATE_fin)] = "FAIL"
                }
            }
        }

    for(z in 2:nrow(fileTABLE)){
        DOWN_MANGO_SEPERATE_fin_frag = DOWN_MANGO_SEPERATE[which(DOWN_MANGO_SEPERATE$condition == z),]

        if(nrow(DOWN_MANGO_SEPERATE_fin_frag)>0){

            for(i in 1:nrow(DOWN_MANGO_SEPERATE_fin_frag)){
                DOWN_MANGO_SEPERATE_fin_DF_frag = as.data.frame(strsplit(DOWN_MANGO_SEPERATE_fin_frag[i,1], split = "//"))
                colnames(DOWN_MANGO_SEPERATE_fin_DF_frag) = "Description"

                for(w in 1:nrow(DOWN_MANGO_SEPERATE_fin_DF_frag)){
                    DOWN_MANGO_SEPERATE_fin_DF_frag[w,1] = gsub("_", " ", DOWN_MANGO_SEPERATE_fin_DF_frag[w,1])
                }

                GO_RAW = read.delim(GOTABLE[z,1], sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
                GO_DF  = GO_RAW[,c(2,5,11)]
                colnames(GO_DF)[2] = paste0("RF_", i)
                GO_DF[is.na(GO_DF)] = 0
                DOWN_MANGO_SEPERATE_fin_DF_TF_frag = dplyr::inner_join(DOWN_MANGO_SEPERATE_fin_DF_frag, GO_DF, by="Description")

                if(nrow(DOWN_MANGO_SEPERATE_fin_DF_TF_frag) < PASSED_NUM){
                    DOWN_MANGO_SEPERATE_fin_frag[i, ncol(DOWN_MANGO_SEPERATE_fin_frag)] = "FAIL"
                    }
                }
            }

        DOWN_MANGO_SEPERATE_fin = rbind(DOWN_MANGO_SEPERATE_fin, DOWN_MANGO_SEPERATE_fin_frag)
    }

    if(nrow(DOWN_MANGO_SEPERATE_fin)>0){DOWN_MANGO_SEPERATE_fin = DOWN_MANGO_SEPERATE_fin[which(DOWN_MANGO_SEPERATE_fin[,ncol(DOWN_MANGO_SEPERATE_fin)] == "PASS"),]}

    if(nrow(DOWN_MANGO_NOT_SEPERATE) > 0){
        DOWN_MANGO_NOT_SEPERATE$label = "PASS"
        
        for(i in 1:nrow(DOWN_MANGO_NOT_SEPERATE)){
            DOWN_MANGO_NOT_SEPERATE_DF = as.data.frame(strsplit(DOWN_MANGO_NOT_SEPERATE[i,1], split = "//"))
            colnames(DOWN_MANGO_NOT_SEPERATE_DF) = "Description"

            for(w in 1:nrow(DOWN_MANGO_NOT_SEPERATE_DF)){
                DOWN_MANGO_NOT_SEPERATE_DF[w,1] = gsub("_", " ", DOWN_MANGO_NOT_SEPERATE_DF[w,1])
            }

            for(s in 1:nrow(GOTABLE)){
                GO_RAW = read.delim(GOTABLE[s,1], sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
                GO_DF  = GO_RAW[,c(2,5,11)]
                colnames(GO_DF)[2] = paste0("RF_", i)
                GO_DF[is.na(GO_DF)] = 0
                DOWN_MANGO_NOT_SEPERATE_DF_TF = dplyr::inner_join(DOWN_MANGO_NOT_SEPERATE_DF, GO_DF, by="Description")

                if(nrow(DOWN_MANGO_NOT_SEPERATE_DF_TF) < PASSED_NUM){
                    DOWN_MANGO_NOT_SEPERATE[i, ncol(DOWN_MANGO_NOT_SEPERATE)] = "FAIL"
                }
            }
        }

        DOWN_MANGO_NOT_SEPERATE = DOWN_MANGO_NOT_SEPERATE[which(DOWN_MANGO_NOT_SEPERATE[, ncol(DOWN_MANGO_NOT_SEPERATE)] == "PASS"),]
    }

if(nrow(fileTABLE)>2){
if(nrow(DOWN_MANGO_NOT_SEPERATE)>0){
    for(i in 1:nrow(DOWN_MANGO_NOT_SEPERATE)){
        DOWN_MANGO_NOT_SEPERATE_frag = DOWN_MANGO_NOT_SEPERATE[i,]
        RFset = t(DOWN_MANGO_NOT_SEPERATE_frag[,2:(nrow(fileTABLE)+1)])

        UPPER = as.data.frame(RFset[order(-RFset[,1]),])[round(nrow(RFset)/4),]
        LOWER = as.data.frame(RFset[order(-RFset[,1]),])[round((nrow(RFset)/4)*3),]

        UPPER = UPPER*FC
        LOWER = LOWER/FC

        for(s in 2:(nrow(fileTABLE)+1)){
    
            if(DOWN_MANGO_NOT_SEPERATE[i,ncol(DOWN_MANGO_NOT_SEPERATE)]=="PASS"){
                if(DOWN_MANGO_NOT_SEPERATE_frag[,s]<LOWER){DOWN_MANGO_NOT_SEPERATE[i,ncol(DOWN_MANGO_NOT_SEPERATE)]="FALSE"}
                if(DOWN_MANGO_NOT_SEPERATE_frag[,s]>UPPER){DOWN_MANGO_NOT_SEPERATE[i,ncol(DOWN_MANGO_NOT_SEPERATE)]="FALSE"}
        
                }
            }
        }

    DOWN_MANGO_NOT_SEPERATE = DOWN_MANGO_NOT_SEPERATE[which(DOWN_MANGO_NOT_SEPERATE[, ncol(DOWN_MANGO_NOT_SEPERATE)] == "PASS"),]
    }
    }


    print("Number of specific tree ls")
    print(nrow(DOWN_MANGO_SEPERATE_fin))

    print("Number of common tree ls")
    print(nrow(DOWN_MANGO_NOT_SEPERATE))

    ttt = " If you want analysis more deeper, Do dynamic analysis"

    if(nrow(DOWN_MANGO_NOT_SEPERATE)==0){
        if(nrow(DOWN_MANGO_SEPERATE_fin)==0){
            ttt = "Do dynamic analysis"
            }
        }
    print(ttt)

    DOWN_MANGO_SEPERATE = rbind(DOWN_MANGO_SEPERATE_fin, DOWN_MANGO_NOT_SEPERATE)
    DOWN_MANGO_SEPERATE = DOWN_MANGO_SEPERATE[,-ncol(DOWN_MANGO_SEPERATE)]

     write.table(DOWN_MANGO_SEPERATE,
                file = outputPATH, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
}
