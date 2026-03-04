#' Draw a heatmap of MANGO results
#'
#' Create a heatmap from MANGO output tables (MANGO_SEPERATE or
#' MANGO_SEPERATE_forMULTI_range).
#'
#' @param filepath Working directory path.
#' @param DEG_list_name DEG list filename under `filepath`.
#' @param trends One of "UP", "DOWN", "SIG" (controls color palette).
#' @param condition "ALL", "SIG", or a specific condition index/value used in the
#'   last column of MANGO_SEPERATE output.
#' @param width Plot width (for notebook display options).
#' @param height Plot height (for notebook display options).
#' @param max Maximum value for heatmap color scale (used in breaks).
#' @param dynamic_analyisis "F" for single/multi heatmap using MANGO_SEPERATE,
#'   "T" for dynamic analysis heatmap using MANGO_SEPERATE_forMULTI_range.
#'
#' @return A `pheatmap` object.
#' @export
#'
#' @examples
#' stopifnot(is.function(MANGO_HeatMap))
#' \donttest{
#' # MANGO_HeatMap(
#' #   filepath = ".",
#' #   DEG_list_name = "input_DEG_list.txt",
#' #   trends = "UP",
#' #   condition = "ALL",
#' #   width = 10,
#' #   height = 6,
#' #   max = 5,
#' #   dynamic_analyisis = "F"
#' # )
#' }
MANGO_HeatMap <- function(filepath,
                          DEG_list_name,
                          trends,
                          condition = "ALL",
                          width,
                          height,
                          max,
                          dynamic_analyisis = "F") {

    ppp <- NULL

    DEG_list_path <- paste0(filepath, "/", DEG_list_name)
    DEG_list_pathDF <- as.data.frame(DEG_list_path)

    TXTbased <- nrow(as.data.frame(DEG_list_pathDF[grepl("txt", DEG_list_pathDF[, 1]), ]))
    CSVbased <- nrow(as.data.frame(DEG_list_pathDF[grepl("csv", DEG_list_pathDF[, 1]), ]))
    BEDbased <- nrow(as.data.frame(DEG_list_pathDF[grepl("bed", DEG_list_pathDF[, 1]), ]))

    if (TXTbased > 0) {
        OUTPUTpath <- paste0(as.data.frame(strsplit(DEG_list_pathDF[, 1], split = "input_DEG_list.txt"))[1, 1],
                             "MANGO_PREPROCESSING_list.txt")
    }
    if (CSVbased > 0) {
        OUTPUTpath <- paste0(as.data.frame(strsplit(DEG_list_pathDF[, 1], split = "input_DEG_list.csv"))[1, 1],
                             "MANGO_PREPROCESSING_list.csv")
    }
    if (BEDbased > 0) {
        OUTPUTpath <- paste0(as.data.frame(strsplit(DEG_list_pathDF[, 1], split = "input_DEG_list.bed"))[1, 1],
                             "MANGO_PREPROCESSING_list.bed")
    }

    MANGO_PREPROCESSING_list_path <- OUTPUTpath
    MANGO_PREPROCESSING_list_pathDF <- as.data.frame(MANGO_PREPROCESSING_list_path)

    if (TXTbased > 0) {
        OUTPUTpath <- paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[, 1], split = "MANGO_PREPROCESSING_list.txt"))[1, 1],
                             "MANGO_SEPERATE.txt")
    }
    if (CSVbased > 0) {
        OUTPUTpath <- paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[, 1], split = "MANGO_PREPROCESSING_list.csv"))[1, 1],
                             "MANGO_SEPERATE.csv")
    }
    if (BEDbased > 0) {
        OUTPUTpath <- paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[, 1], split = "MANGO_PREPROCESSING_list.bed"))[1, 1],
                             "MANGO_SEPERATE.bed")
    }
    MANGO_SEPERATE_path <- OUTPUTpath

    if (TXTbased > 0) {
        OUTPUTpath <- paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[, 1], split = "MANGO_PREPROCESSING_list.txt"))[1, 1],
                             "MANGO_SEPERATE_forMULTI_range.txt")
    }
    if (CSVbased > 0) {
        OUTPUTpath <- paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[, 1], split = "MANGO_PREPROCESSING_list.csv"))[1, 1],
                             "MANGO_SEPERATE_forMULTI_range.csv")
    }
    if (BEDbased > 0) {
        OUTPUTpath <- paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[, 1], split = "MANGO_PREPROCESSING_list.bed"))[1, 1],
                             "MANGO_SEPERATE_forMULTI_range.bed")
    }
    MANGO_SEPERATE_forMULTI_range_path <- OUTPUTpath

    fileTABLE_path <- MANGO_PREPROCESSING_list_path

    if (dynamic_analyisis == "F") {

        input_SEPERATE <- MANGO_SEPERATE_path

        options(repr.plot.width = width, repr.plot.height = height,
                repr.plot.res = 400, repr.plot.pointsize = 10)

        INPUT <- utils::read.table(input_SEPERATE, sep = "\t", header = TRUE, fill = TRUE)
        fileTABLE <- utils::read.table(fileTABLE_path, sep = "\t", header = TRUE, fill = TRUE)
        INPUTforNUM <- nrow(fileTABLE)

        INPUT_BB <- INPUT[1, -ncol(INPUT)]

        if (condition == "ALL") {
            INPUT <- INPUT[, -ncol(INPUT)]
            INPUT <- unique(INPUT)
        }

        if (condition == "SIG") {
            for (i in seq_len(INPUTforNUM)) {
                INPUT_BB_frag <- INPUT[which(INPUT$condition == i), -ncol(INPUT)]
                INPUT_BB_frag <- unique(INPUT_BB_frag)
                INPUT_BB <- rbind(INPUT_BB, INPUT_BB_frag)
                INPUT_BB <- unique(INPUT_BB)
            }
            INPUT <- INPUT_BB
        }

        if (condition != "ALL" && condition != "SIG") {
            INPUT <- INPUT[which(INPUT[, ncol(INPUT)] == condition), -ncol(INPUT)]
            INPUT <- unique(INPUT)
        }

        if (trends == "UP") {
            col_set <- grDevices::colorRampPalette(
                c("white", "#F5B78EFF", "#F19C7CFF", "#EA8171FF", "#DD686CFF", "#CA5268FF", "#B13F64FF")
            )(100)
        }
        if (trends == "DOWN") {
            col_set <- grDevices::colorRampPalette(
                c("white", "#CCFDFFFF", "#99F8FFFF", "#66F0FFFF", "#33E4FFFF", "#00AACCFF", "#007A99FF")
            )(100)
        }
        if (trends == "SIG") {
            col_set <- grDevices::colorRampPalette(
                c("white", "tan1", "tan2", "tan3", "orange", "brown")
            )(100)
        }

        namen <- INPUT[, 1]
        conte <- INPUT[, -1]
        rownames(conte) <- as.vector(namen)

        PCA <- conte
        ppp <- pheatmap::pheatmap(
            PCA,
            cluster_rows = FALSE,
            cluster_cols = FALSE,
            show_rownames = FALSE,
            show_colnames = FALSE,
            legend = FALSE,
            border_color = "grey80",
            breaks = seq(from = 0, to = max, length.out = 100),
            color = col_set
        )
    }

    if (dynamic_analyisis == "T") {

        input_SEPERATE <- MANGO_SEPERATE_forMULTI_range_path

        options(repr.plot.width = width, repr.plot.height = height,
                repr.plot.res = 400, repr.plot.pointsize = 10)

        INPUT <- utils::read.table(input_SEPERATE, sep = "\t", header = TRUE, fill = TRUE)
        INPUT <- INPUT[, -ncol(INPUT)]
        INPUT <- unique(INPUT)

        if (trends == "UP") {
            col_set <- grDevices::colorRampPalette(
                c("grey80", "#F5B78EFF", "#F19C7CFF", "#EA8171FF", "#DD686CFF", "#CA5268FF", "#B13F64FF")
            )(100)
        }
        if (trends == "DOWN") {
            col_set <- grDevices::colorRampPalette(
                c("grey80", "#CCFDFFFF", "#99F8FFFF", "#66F0FFFF", "#33E4FFFF", "#00AACCFF", "#007A99FF")
            )(100)
        }
        if (trends == "SIG") {
            col_set <- grDevices::colorRampPalette(
                c("grey80", "tan1", "tan2", "tan3", "orange", "brown")
            )(100)
        }

        namen <- INPUT[, 1]
        conte <- INPUT[, -1]
        rownames(conte) <- as.vector(namen)

        PCA <- conte
        ppp <- pheatmap::pheatmap(
            PCA,
            cluster_rows = FALSE,
            cluster_cols = FALSE,
            show_rownames = FALSE,
            show_colnames = FALSE,
            legend = FALSE,
            border_color = "grey80",
            breaks = seq(from = 0, to = max, length.out = 100),
            color = col_set
        )
    }

    ppp
}
