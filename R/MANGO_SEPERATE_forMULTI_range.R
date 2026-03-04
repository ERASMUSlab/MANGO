#' Select significant MANGO terms within a condition range
#'
#' Filters MANGO comparison results to terms present in all selected conditions
#' (\code{condition_A:condition_B}) and exceeding an \code{FC} ratio against
#' the remaining conditions.
#'
#' @param input_COMPARE Path to the MANGO comparison table.
#' @param fileTABLE_path Path to a table listing conditions (used to determine number of columns).
#' @param FC Fold-change ratio cutoff used for filtering.
#' @param outputPATH Output path to write the filtered table.
#' @param condition_A Start condition index (1-based, matching your code convention).
#' @param condition_B End condition index (1-based, matching your code convention).
#'
#' @return Writes \code{outputPATH}. Returns \code{NULL} invisibly.
#'
#' @examples
#' stopifnot(is.function(MANGO_SEPERATE_forMULTI_range))
#' @export


MANGO_SEPERATE_forMULTI_range = function(input_COMPARE,
                                         fileTABLE_path,
                                         FC,
                                         outputPATH,
                                         condition_A,
                                         condition_B){

    DOWN_MANGO_COMPARE = utils::read.table(input_COMPARE, sep="\t", header=T, fill=TRUE)
    DOWN_MANGO_COMPARE = unique(DOWN_MANGO_COMPARE)

    fileTABLE = utils::read.table(fileTABLE_path, sep="\t", header=T, fill=TRUE)

    DOWN_MANGO_SEPERATE_tmp = DOWN_MANGO_COMPARE

    for(i in condition_A:condition_B){
        DOWN_MANGO_SEPERATE_tmp = DOWN_MANGO_SEPERATE_tmp[which(DOWN_MANGO_SEPERATE_tmp[,(i+1)] > 0),]
    }

    real_condition_A_pos = condition_A+1
    real_condition_B_pos = condition_B+1

    DOWN_MANGO_SEPERATE_tmp1 = DOWN_MANGO_SEPERATE_tmp[,c(1,real_condition_A_pos:real_condition_B_pos)]
    DOWN_MANGO_SEPERATE_tmp2 = DOWN_MANGO_SEPERATE_tmp[,c(-real_condition_A_pos:-real_condition_B_pos)]

    DOWN_MANGO_SEPERATE_sig_term = DOWN_MANGO_COMPARE[,1]
    DOWN_MANGO_SEPERATE_sig_term = as.data.frame(DOWN_MANGO_SEPERATE_sig_term)
    colnames(DOWN_MANGO_SEPERATE_sig_term) = "MANGO_TERM"

    finnum1 = ncol(DOWN_MANGO_SEPERATE_tmp1)-1

    DOWN_MANGO_SEPERATE_tmp2_tmp = DOWN_MANGO_SEPERATE_tmp2
    
    finnum2 = ncol(DOWN_MANGO_SEPERATE_tmp2_tmp)-1

    for(m in 1:finnum1){

        DOWN_MANGO_SEPERATE_tmp2_tmp = DOWN_MANGO_SEPERATE_tmp2
    
        for(n in 1:finnum2){
            DOWN_MANGO_SEPERATE_tmp2_tmp[,(n+1)] = DOWN_MANGO_SEPERATE_tmp1[,(m+1)]/(DOWN_MANGO_SEPERATE_tmp2_tmp[,(n+1)]+0.0000000001)
        }

        for(n in 1:finnum2){
            DOWN_MANGO_SEPERATE_tmp2_tmp = DOWN_MANGO_SEPERATE_tmp2_tmp[which(DOWN_MANGO_SEPERATE_tmp2_tmp[,(n+1)] >= FC),]
        }

    DOWN_MANGO_SEPERATE_sig_term_frag = DOWN_MANGO_SEPERATE_tmp2_tmp[,1]
    DOWN_MANGO_SEPERATE_sig_term_frag = as.data.frame(DOWN_MANGO_SEPERATE_sig_term_frag)
    colnames(DOWN_MANGO_SEPERATE_sig_term_frag) = "MANGO_TERM"

    DOWN_MANGO_SEPERATE_sig_term = dplyr::inner_join(DOWN_MANGO_SEPERATE_sig_term,
                                              DOWN_MANGO_SEPERATE_sig_term_frag,
                                              by = "MANGO_TERM")
        }

    DOWN_MANGO_SEPERATE_sig = dplyr::inner_join(DOWN_MANGO_SEPERATE_sig_term,DOWN_MANGO_COMPARE,by = "MANGO_TERM")

    utils::write.table(DOWN_MANGO_SEPERATE_sig,
                file = outputPATH, col.names=T, row.names=F, quote=F,sep="\t")
    }

