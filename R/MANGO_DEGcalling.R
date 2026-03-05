#' Call DEGs with DESeq2 and write per-contrast DEG files for MANGO
#'
#' Runs DESeq2 for multiple contrasts defined by \code{mango_design}, writes
#' normalized counts and DEG result tables, and appends selected DEG file paths
#' into \code{input_DEG_list.txt} under \code{filepath}.
#'
#' @param filepath Directory path containing input files (expects \code{count/RAWcount.txt}).
#' @param type DEG calling mode. One of \code{"standard"} or \code{"broad"}.
#'   \code{"standard"} uses \code{adjp=0.05} and \code{FC=2}; \code{"broad"} uses
#'   \code{adjp=0.05} and \code{FC=1.5}.
#' @param full_condition Character vector of condition names (must match sample condition labels).
#' @param number_of_rep Integer vector of replicate counts for each condition in \code{full_condition}.
#' @param mango_design Character vector specifying contrasts and direction to export.
#'   Each element should be formatted like \code{"<objective>_<control>_<trend>"},
#'   where \code{trend} is one of \code{"UP"}, \code{"DOWN"}, \code{"SIG"}.
#'
#' @return Invisibly returns \code{NULL}. Results are written to files under \code{filepath}.
#' @export
#'
#' @importFrom utils read.table write.table
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results counts
#' @importFrom SummarizedExperiment assay
MANGO_DEGcalling = function(filepath,
                            type = "standard",
                            full_condition,
                            number_of_rep,
                            mango_design){

    MANGO_CMD = paste0(": > ", filepath, "/input_DEG_list.txt")
    base::system(MANGO_CMD, intern = TRUE)

    if(type == "standard"){
        adjp = 0.05
        FC = 2
    }

    if(type == "broad"){
        adjp = 0.05
        FC = 1.5
    }

    mango_design_DF = as.data.frame(mango_design)

    s = 1
    mango_design_DF_split = as.data.frame(strsplit(mango_design_DF[s,1], split = "_"))

    for(i in 2:nrow(mango_design_DF)){
        mango_design_DF_split_frag = as.data.frame(strsplit(mango_design_DF[i,1], split = "_"))
        mango_design_DF_split = cbind(mango_design_DF_split, mango_design_DF_split_frag)
    }

    for(i in 1:nrow(mango_design_DF)){
        colnames(mango_design_DF_split)[i] = paste0("case", i)
    }

    for(e in 1:nrow(mango_design_DF)){
        control_of_sample = mango_design_DF_split[2, e]
        objective_of_sample = mango_design_DF_split[1, e]
        trends = mango_design_DF_split[3, e]

        input_matrix_path = paste(filepath, "/count/RAWcount.txt", sep = "")
        input_matrix = utils::read.table(input_matrix_path, sep = "\t", header = TRUE)

        analysis_name = paste("condition_", objective_of_sample, "_vs_", control_of_sample, sep = "")
        numcer_of_rep = number_of_rep
        numcer_of_rep_forcal = as.data.frame(numcer_of_rep)

        input_precondition = rep(full_condition[1], numcer_of_rep[1])
        input_precondition = as.data.frame(input_precondition)
        colnames(input_precondition) = "precondition1"

        for(i in 2:nrow(numcer_of_rep_forcal)){
            input_precondition_frag = rep(full_condition[i], numcer_of_rep[i])
            input_precondition_frag = as.data.frame(input_precondition_frag)
            colnames(input_precondition_frag) = "precondition1"
            input_precondition = rbind(input_precondition, input_precondition_frag)
        }

        parametere_adjp = adjp
        parametere_LFC = log2(FC)

        raw = input_matrix
        rawdata = raw[, -1]
        rawN = as.vector(raw[, 1])
        rownames(rawdata) = rawN
        rawdata = as.matrix(rawdata)

        precondition1 = input_precondition[, 1]

        condition = as.factor(precondition1)
        condition = as.data.frame(condition)
        metadata = condition

        mata_name = colnames(raw[, -1])
        rownames(metadata) = mata_name

        dds <- suppressMessages(
            DESeq2::DESeqDataSetFromMatrix(countData = rawdata,
                                          colData = metadata,
                                          design = ~ condition)
        )

        keep <- rowSums(SummarizedExperiment::assay(dds) >= 0) >= (ncol(rawdata) / 2) + 1
        dds <- dds[keep, ]
        dds <- suppressMessages(DESeq2::DESeq(dds))

        normalized_counts_tmp <- DESeq2::counts(dds, normalized = TRUE)
        gene = rownames(normalized_counts_tmp)
        gene = as.data.frame(gene)
        DESeq2_Norm = cbind(gene, normalized_counts_tmp)
        rownames(DESeq2_Norm) = NULL

        dds$condition <- stats::relevel(dds$condition, ref = control_of_sample)
        dds <- suppressMessages(DESeq2::DESeq(dds))
        res <- suppressMessages(DESeq2::results(dds, name = analysis_name, pAdjustMethod = "BH"))

        resOrdered <- res[order(res$padj), ]

        final = as.data.frame(resOrdered)
        gene = rownames(resOrdered)
        gene = as.data.frame(gene)
        Full_data = cbind(gene, resOrdered)
        rownames(Full_data) = NULL

        resSig <- subset(resOrdered, padj < parametere_adjp)
        final <- subset(resSig, abs(log2FoldChange) > parametere_LFC)
        final = as.data.frame(final)
        gene = rownames(final)
        gene = as.data.frame(gene)
        sig_data = cbind(gene, final)
        rownames(sig_data) = NULL

        DEG_UP = sig_data[1, ]

        for(i in 1:nrow(sig_data)){
            RNA_LFC = sig_data[i, 3]

            if(RNA_LFC > 0){
                DEG_UP_row = sig_data[i, ]
                DEG_UP_row = as.data.frame(DEG_UP_row)
                DEG_UP = rbind(DEG_UP, DEG_UP_row)
            }
        }

        DEG_UP = DEG_UP[-1, ]

        DEG_DOWN = sig_data[1, ]

        for(i in 1:nrow(sig_data)){
            RNA_LFC = sig_data[i, 3]

            if(RNA_LFC < 0){
                DEG_DOWN_row = sig_data[i, ]
                DEG_DOWN_row = as.data.frame(DEG_DOWN_row)
                DEG_DOWN = rbind(DEG_DOWN, DEG_DOWN_row)
            }
        }

        DEG_DOWN = DEG_DOWN[-1, ]

        make_DIR_command = paste("mkdir ", filepath, "/DEG_result", sep = "")
        base::system(make_DIR_command, intern = TRUE)

        make_DIR_command = paste("mkdir ", filepath, "/DEG_result/", analysis_name, "_", adjp, "_", FC, sep = "")
        base::system(make_DIR_command, intern = TRUE)

        DESeq2_Norm_path = paste(filepath, "/count/DESeq2_normalized.txt", sep = "")
        Full_data_path = paste(filepath, "/DEG_result/", analysis_name, "_", adjp, "_", FC, "/DEG_FULL.txt", sep = "")
        sig_data_path = paste(filepath, "/DEG_result/", analysis_name, "_", adjp, "_", FC, "/DEG_SIG.txt", sep = "")
        DEG_UP_path = paste(filepath, "/DEG_result/", analysis_name, "_", adjp, "_", FC, "/DEG_UP.txt", sep = "")
        DEG_DOWN_path = paste(filepath, "/DEG_result/", analysis_name, "_", adjp, "_", FC, "/DEG_DOWN.txt", sep = "")

        if(trends == "UP"){ mainPATH = DEG_UP_path }
        if(trends == "DOWN"){ mainPATH = DEG_DOWN_path }
        if(trends == "SIG"){ mainPATH = sig_data_path }

        MANGO_CMD = paste0("echo ", mainPATH, " >>", filepath, "/input_DEG_list.txt")
        base::system(MANGO_CMD, intern = TRUE)

        utils::write.table(DESeq2_Norm,
                           file = DESeq2_Norm_path,
                           col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

        utils::write.table(Full_data,
                           file = Full_data_path,
                           col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

        utils::write.table(sig_data,
                           file = sig_data_path,
                           col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

        utils::write.table(DEG_UP,
                           file = DEG_UP_path,
                           col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

        utils::write.table(DEG_DOWN,
                           file = DEG_DOWN_path,
                           col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
    }

    invisible(NULL)
}
