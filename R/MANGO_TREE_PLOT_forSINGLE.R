#' Plot MANGO trees for a single condition
#'
#' Visualize tree-level summaries from `MANGO_SEPERATE` output for a single-case
#' analysis. Supports bar plot and circular (polar) plot styles.
#'
#' @param filepath Working directory path.
#' @param trends One of "UP", "DOWN", "SIG" (controls color palette).
#' @param width Plot width (for notebook display options).
#' @param height Plot height (for notebook display options).
#' @param plot Plot type: "bar" or "cir".
#' @param DEG_list_name DEG list filename under `filepath`.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' stopifnot(is.function(MANGO_TREE_PLOT_forSINGLE))
#' \donttest{
#' # MANGO_TREE_PLOT_forSINGLE(
#' #   filepath = ".",
#' #   trends = "UP",
#' #   width = 10,
#' #   height = 6,
#' #   plot = "bar",
#' #   DEG_list_name = "input_DEG_list.txt"
#' # )
#' }
MANGO_TREE_PLOT_forSINGLE <- function(filepath,
                                      trends,
                                      width,
                                      height,
                                      plot,
                                      DEG_list_name) {

    condition <- 1L

    blacklist <- c("regulation", "positive", "negative", "a", "the", "of",
                   "gene", "type", "I", "II", "to", "by")
    blacklist <- as.data.frame(blacklist)
    colnames(blacklist) <- "Description"

    DEG_list_path <- paste0(filepath, "/", DEG_list_name)
    DEG_list_pathDF <- as.data.frame(DEG_list_path)

    TXTbased <- nrow(as.data.frame(DEG_list_pathDF[grepl("txt", DEG_list_pathDF[, 1]), ]))
    CSVbased <- nrow(as.data.frame(DEG_list_pathDF[grepl("csv", DEG_list_pathDF[, 1]), ]))
    BEDbased <- nrow(as.data.frame(DEG_list_pathDF[grepl("bed", DEG_list_pathDF[, 1]), ]))

    if (TXTbased > 0) {
        OUTPUTpath <- paste0(as.data.frame(strsplit(DEG_list_pathDF[, 1], split = "_list.txt"))[1, 1], "_GO_list.txt")
    }
    if (CSVbased > 0) {
        OUTPUTpath <- paste0(as.data.frame(strsplit(DEG_list_pathDF[, 1], split = "_list.csv"))[1, 1], "_GO_list.csv")
    }
    if (BEDbased > 0) {
        OUTPUTpath <- paste0(as.data.frame(strsplit(DEG_list_pathDF[, 1], split = "_list.bed"))[1, 1], "_GO_list.bed")
    }
    GO_list_path <- OUTPUTpath

    if (TXTbased > 0) {
        OUTPUTpath <- paste0(as.data.frame(strsplit(DEG_list_pathDF[, 1], split = "input_DEG_list.txt"))[1, 1], "MANGO_PREPROCESSING_list.txt")
    }
    if (CSVbased > 0) {
        OUTPUTpath <- paste0(as.data.frame(strsplit(DEG_list_pathDF[, 1], split = "input_DEG_list.csv"))[1, 1], "MANGO_PREPROCESSING_list.csv")
    }
    if (BEDbased > 0) {
        OUTPUTpath <- paste0(as.data.frame(strsplit(DEG_list_pathDF[, 1], split = "input_DEG_list.bed"))[1, 1], "MANGO_PREPROCESSING_list.bed")
    }
    MANGO_PREPROCESSING_list_path <- OUTPUTpath
    MANGO_PREPROCESSING_list_pathDF <- as.data.frame(MANGO_PREPROCESSING_list_path)

    if (TXTbased > 0) {
        OUTPUTpath <- paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[, 1], split = "MANGO_PREPROCESSING_list.txt"))[1, 1], "MANGO_SEPERATE.txt")
    }
    if (CSVbased > 0) {
        OUTPUTpath <- paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[, 1], split = "MANGO_PREPROCESSING_list.csv"))[1, 1], "MANGO_SEPERATE.csv")
    }
    if (BEDbased > 0) {
        OUTPUTpath <- paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[, 1], split = "MANGO_PREPROCESSING_list.bed"))[1, 1], "MANGO_SEPERATE.bed")
    }
    MANGO_SEPERATE_path <- OUTPUTpath

    if (TXTbased > 0) {
        OUTPUTpath <- paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[, 1], split = "MANGO_PREPROCESSING_list.txt"))[1, 1], "MANGO_SEPERATE_forMULTI_range.txt")
    }
    if (CSVbased > 0) {
        OUTPUTpath <- paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[, 1], split = "MANGO_PREPROCESSING_list.csv"))[1, 1], "MANGO_SEPERATE_forMULTI_range.csv")
    }
    if (BEDbased > 0) {
        OUTPUTpath <- paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[, 1], split = "MANGO_PREPROCESSING_list.bed"))[1, 1], "MANGO_SEPERATE_forMULTI_range.bed")
    }
    MANGO_SEPERATE_forMULTI_range_path <- OUTPUTpath

    fileTABLE_path <- MANGO_PREPROCESSING_list_path
    GOTABLE_path <- GO_list_path

    fileTABLE <- utils::read.table(fileTABLE_path, sep = "\t", header = TRUE, fill = TRUE)

    MANGO_SEPERATE_INPUT <- MANGO_SEPERATE_path
    SEPERATE_INPUT <- utils::read.table(MANGO_SEPERATE_INPUT, sep = "\t", header = TRUE, fill = TRUE)
    SEPERATE_INPUT <- SEPERATE_INPUT[order(-SEPERATE_INPUT[, 2]), ]
    SEPERATE_INPUT <- stats::na.omit(SEPERATE_INPUT)

    fileTABLE <- utils::read.table(fileTABLE_path, sep = "\t", header = TRUE, fill = TRUE)
    GOTABLE <- utils::read.table(GOTABLE_path, sep = "\t", header = TRUE, fill = TRUE)

    i <- 1L
    GO_RAW <- utils::read.delim(GOTABLE[i, 1], sep = "\t", header = TRUE,
                               quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
    GO_DF <- GO_RAW[, c(2, 5, 11)]
    colnames(GO_DF)[2] <- paste0("RF_", i)
    GO_DF[is.na(GO_DF)] <- 0

    i <- 1L
    SEPERATE_INPUT_DF <- as.data.frame(strsplit(SEPERATE_INPUT[i, 1], split = "//"))
    colnames(SEPERATE_INPUT_DF) <- "Description"
    for (s in seq_len(nrow(SEPERATE_INPUT_DF))) {
        SEPERATE_INPUT_DF[s, 1] <- gsub("_", " ", SEPERATE_INPUT_DF[s, 1])
    }

    DATA_RAW_DF <- dplyr::inner_join(SEPERATE_INPUT_DF, GO_DF, by = "Description")
    DATA_RAW_DF <- DATA_RAW_DF[order(-DATA_RAW_DF[, (condition + 1)]), ]
    DATA_RAW_DF <- cbind(DATA_RAW_DF, i)
    colnames(DATA_RAW_DF) <- c("Description", "RF", "gene", "cluster")

    if (nrow(SEPERATE_INPUT) >= 2) {
        for (i in 2:nrow(SEPERATE_INPUT)) {
            SEPERATE_INPUT_DF <- as.data.frame(strsplit(SEPERATE_INPUT[i, 1], split = "//"))
            colnames(SEPERATE_INPUT_DF) <- "Description"
            for (s in seq_len(nrow(SEPERATE_INPUT_DF))) {
                SEPERATE_INPUT_DF[s, 1] <- gsub("_", " ", SEPERATE_INPUT_DF[s, 1])
            }

            DATA_RAW_DF_frag <- dplyr::inner_join(SEPERATE_INPUT_DF, GO_DF, by = "Description")
            DATA_RAW_DF_frag <- DATA_RAW_DF_frag[order(-DATA_RAW_DF_frag[, (condition + 1)]), ]
            DATA_RAW_DF_frag <- cbind(DATA_RAW_DF_frag, i)
            colnames(DATA_RAW_DF_frag) <- c("Description", "RF", "gene", "cluster")
            DATA_RAW_DF <- rbind(DATA_RAW_DF, DATA_RAW_DF_frag)
        }
    }

    DATA_RAW_DF[, 5] <- 0
    colnames(DATA_RAW_DF)[5] <- "score"

    for (i in seq_len(nrow(DATA_RAW_DF))) {
        gene <- as.data.frame(strsplit(DATA_RAW_DF[i, 3], split = "/"))
        colnames(gene) <- "gene"

        DEG <- utils::read.table(DEG_list_pathDF[1, 1], sep = "\t", header = FALSE, fill = TRUE)
        DEG <- utils::read.table(DEG[1, 1], sep = "\t", header = TRUE, fill = TRUE)

        DEG <- dplyr::inner_join(DEG, gene, by = "gene")

        DATA_RAW_DF[i, 5] <- 1 / (log10(DATA_RAW_DF[i, 2]) * (-1)) * stats::median(DEG[, 3])
    }

    i <- 1L
    subset <- DATA_RAW_DF[which(DATA_RAW_DF[, 4] == i), ]
    termset <- as.data.frame(subset[, 1])

    termset_frag <- as.data.frame(strsplit(termset[1, 1], split = " "))
    colnames(termset_frag) <- "Description"
    if (nrow(termset) >= 2) {
        for (m in 2:nrow(termset)) {
            termset_fragfrag <- as.data.frame(strsplit(termset[m, 1], split = " "))
            colnames(termset_fragfrag) <- "Description"
            termset_frag <- rbind(termset_frag, termset_fragfrag)
        }
    }

    termset_frag <- dplyr::anti_join(termset_frag, blacklist, by = "Description")
    termset_fraguniq <- termset_frag
    termset_fraguniq$hit <- 0
    termset_fraguniq <- unique(termset_fraguniq)

    for (s in seq_len(nrow(termset_fraguniq))) {
        termset_fraguniq[s, 2] <- nrow(as.data.frame(termset_frag[which(termset_frag[, 1] == termset_fraguniq[s, 1]), ]))
    }

    CN <- paste(termset_fraguniq[1, 1], termset_fraguniq[2, 1], sep = " ")
    TreeScore <- SEPERATE_INPUT[i, 2]

    geneset <- as.data.frame(subset[, 3])
    geneset_frag <- as.data.frame(strsplit(geneset[1, 1], split = "/"))
    colnames(geneset_frag) <- "Description"
    if (nrow(geneset) >= 2) {
        for (m in 2:nrow(geneset)) {
            geneset_fragfrag <- as.data.frame(strsplit(geneset[m, 1], split = "/"))
            colnames(geneset_fragfrag) <- "Description"
            geneset_frag <- rbind(geneset_frag, geneset_fragfrag)
        }
    }

    geneset_frag <- unique(geneset_frag)
    genenum <- nrow(geneset_frag)

    DEG <- utils::read.table(DEG_list_pathDF[1, 1], sep = "\t", header = FALSE, fill = TRUE)
    DEG <- utils::read.table(DEG[1, 1], sep = "\t", header = TRUE, fill = TRUE)
    colnames(DEG)[1] <- "Description"
    DEG <- dplyr::inner_join(DEG, geneset_frag, by = "Description")
    geneS <- stats::median(DEG[, 3])

    termnum <- nrow(subset)

    DATA_DATAFRAME <- cbind(CN, TreeScore, genenum, termnum, geneS)
    colnames(DATA_DATAFRAME) <- c("region", "sum_length", "mean_gain", "n", "geneS")
    DATA_DATAFRAME <- as.data.frame(DATA_DATAFRAME)

    MAX <- max(DATA_RAW_DF[, 4])

    if (MAX >= 2) {
        for (i in 2:MAX) {
            subset <- DATA_RAW_DF[which(DATA_RAW_DF[, 4] == i), ]
            termset <- as.data.frame(subset[, 1])

            termset_frag <- as.data.frame(strsplit(termset[1, 1], split = " "))
            colnames(termset_frag) <- "Description"

            if (nrow(termset) >= 2) {
                for (m in 2:nrow(termset)) {
                    termset_fragfrag <- as.data.frame(strsplit(termset[m, 1], split = " "))
                    colnames(termset_fragfrag) <- "Description"
                    termset_frag <- rbind(termset_frag, termset_fragfrag)
                }
            }

            termset_frag <- dplyr::anti_join(termset_frag, blacklist, by = "Description")
            termset_fraguniq <- termset_frag
            termset_fraguniq$hit <- 0
            termset_fraguniq <- unique(termset_fraguniq)

            for (s in seq_len(nrow(termset_fraguniq))) {
                termset_fraguniq[s, 2] <- nrow(as.data.frame(termset_frag[which(termset_frag[, 1] == termset_fraguniq[s, 1]), ]))
            }

            CN <- paste(termset_fraguniq[1, 1], termset_fraguniq[2, 1], sep = " ")
            TreeScore <- SEPERATE_INPUT[i, 2]

            geneset <- as.data.frame(subset[, 3])
            geneset_frag <- as.data.frame(strsplit(geneset[1, 1], split = "/"))
            colnames(geneset_frag) <- "Description"

            if (nrow(geneset) >= 2) {
                for (m in 2:nrow(geneset)) {
                    geneset_fragfrag <- as.data.frame(strsplit(geneset[m, 1], split = "/"))
                    colnames(geneset_fragfrag) <- "Description"
                    geneset_frag <- rbind(geneset_frag, geneset_fragfrag)
                }
            }

            geneset_frag <- unique(geneset_frag)
            genenum <- nrow(geneset_frag)

            DEG <- utils::read.table(DEG_list_pathDF[1, 1], sep = "\t", header = FALSE, fill = TRUE)
            DEG <- utils::read.table(DEG[1, 1], sep = "\t", header = TRUE, fill = TRUE)
            colnames(DEG)[1] <- "Description"
            DEG <- dplyr::inner_join(DEG, geneset_frag, by = "Description")
            geneS <- stats::median(DEG[, 3])

            termnum <- nrow(subset)

            DATA_DATAFRAME_frag <- cbind(CN, TreeScore, genenum, termnum, geneS)
            colnames(DATA_DATAFRAME_frag) <- c("region", "sum_length", "mean_gain", "n", "geneS")
            DATA_DATAFRAME_frag <- as.data.frame(DATA_DATAFRAME_frag)

            DATA_DATAFRAME <- rbind(DATA_DATAFRAME, DATA_DATAFRAME_frag)
        }
    }

    plot_df <- DATA_DATAFRAME
    plot_df[, 2] <- as.numeric(plot_df[, 2])
    plot_df[, 3] <- as.numeric(plot_df[, 3])
    plot_df[, 4] <- as.numeric(plot_df[, 4])

    totgene <- as.data.frame(strsplit(DATA_RAW_DF[1, 3], split = "/"))
    colnames(totgene) <- "gene"

    if (nrow(DATA_RAW_DF) >= 2) {
        for (i in 2:nrow(DATA_RAW_DF)) {
            totgene_frag <- as.data.frame(strsplit(DATA_RAW_DF[i, 3], split = "/"))
            colnames(totgene_frag) <- "gene"
            totgene <- rbind(totgene, totgene_frag)
        }
    }
    totgene <- unique(totgene)

    DEG <- utils::read.table(DEG_list_pathDF[1, 1], sep = "\t", header = FALSE, fill = TRUE)
    DEG <- utils::read.table(DEG[1, 1], sep = "\t", header = TRUE, fill = TRUE)

    plot_df[, 3] <- (plot_df[, 3] / nrow(DEG) * 10)

    MAXpara_1 <- max(plot_df[, 2])
    MAXpara_2 <- max(plot_df[, 3])

    plot_df_set <- unique(as.data.frame(plot_df[, 1]))

    for (q in seq_len(nrow(plot_df_set))) {
        tt <- as.data.frame(plot_df[which(plot_df[, 1] == plot_df_set[q, 1]), ])
        if (nrow(tt) > 1) {
            for (w in seq_len(nrow(tt))) {
                tt[w, 1] <- paste0(tt[w, 1], " tree", w)
            }
            plot_df <- as.data.frame(plot_df[which(plot_df[, 1] != plot_df_set[q, 1]), ])
            plot_df <- rbind(plot_df, tt)
        }
    }

    if (MAXpara_1 > MAXpara_2) MAXpara <- MAXpara_1
    if (MAXpara_2 > MAXpara_1) MAXpara <- MAXpara_2

    plot_df$tree <- seq_len(nrow(plot_df))
    plot_df[, 5] <- as.numeric(plot_df[, 5])

    plot_bar <- plot_df
    plot_cir <- plot_df

    plot_cir <- plot_cir %>%
        dplyr::mutate(xlab = reorder(stringr::str_wrap(region, 5), -sum_length))

    for (s in seq_len(nrow(plot_bar))) {
        plot_bar[s, 1] <- gsub(" ", " & ", plot_bar[s, 1])
    }

    plot_bar <- plot_bar[order(-plot_bar$sum_length), ]
    plot_bar$tree <- seq_len(nrow(plot_bar))

    plot_cir <- plot_cir[order(-plot_cir$sum_length), ]
    plot_cir$tree <- seq_len(nrow(plot_cir))

    MAXpara <- ceiling(MAXpara)
    set <- ceiling(MAXpara / 3)

    if (trends == "UP") {
        basecol <- "#DD686CFF"
        fillcol <- "pink2"
        fontcol <- "darkred"
    }

    if (trends == "DOWN") {
        basecol <- "#33E4FFFF"
        fillcol <- "skyblue2"
        fontcol <- "darkblue"
    }

    if (trends == "SIG") {
        basecol <- "tan"
        fillcol <- "tan2"
        fontcol <- "brown"
    }

    if (trends == "UP") colorset <- c("#F5B78EFF", "#F19C7CFF", "#EA8171FF", "#DD686CFF", "#CA5268FF", "#B13F64FF")
    if (trends == "DOWN") colorset <- c("#CCFDFFFF", "#99F8FFFF", "#66F0FFFF", "#33E4FFFF", "#00AACCFF", "#007A99FF")

    if (plot == "bar") {

        ppp <- ggplot2::ggplot(plot_bar, ggplot2::aes(sum_length, forcats::fct_reorder(region, sum_length), fill = n)) +
            ggplot2::scale_fill_gradientn("Number of terms", colours = colorset) +
            ggplot2::geom_bar(stat = "identity",
                              position = ggplot2::position_dodge(width = 0.8),
                              colour = "black",
                              width = 0.85,
                              linewidth = 0.1) +
            ggplot2::geom_segment(ggplot2::aes(x = 0, y = region, xend = mean_gain, yend = region),
                                  linetype = "dashed", color = "gray12") +
            ggplot2::geom_point(ggplot2::aes(x = mean_gain, y = region, size = geneS),
                                color = "gray12", show.legend = TRUE) +
            theme_dose(12) +
            ggplot2::ylab("\n\n\n") +
            ggplot2::xlab("HWES & Gene distribution ratio") +
            ggplot2::scale_size_continuous("FC of genes", range = c(1, 5)) +
            ggplot2::facet_grid(tree ~ ., scales = "free_y", space = "free_y") +
            ggplot2::theme(
                strip.text.y = ggplot2::element_text(face = "bold", size = 14, color = "darkblue"),
                strip.background = ggplot2::element_rect(fill = fillcol, color = fontcol, linewidth = 0.6)
            ) +
            ggplot2::theme(
                plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 25, color = "darkblue"),
                axis.text.x = ggplot2::element_text(face = "bold", size = 7),
                axis.text.y = ggplot2::element_text(face = "bold", size = 15),
                axis.title.x = ggplot2::element_text(face = "bold", size = 15),
                legend.title = ggplot2::element_text(face = "bold", size = 10),
                legend.text = ggplot2::element_text(face = "bold", size = 8)
            )
    }

    if (plot == "cir") {

        ppp <- ggplot2::ggplot(plot_cir) +
            ggplot2::geom_hline(ggplot2::aes(yintercept = y),
                                data = data.frame(y = seq_len(4) - 1L) * set,
                                color = "lightgrey") +
            ggplot2::geom_col(
                ggplot2::aes(x = reorder(stringr::str_wrap(region, 5), -sum_length),
                             y = sum_length, fill = n),
                position = "dodge2", show.legend = TRUE, alpha = 0.9
            ) +
            ggplot2::geom_point(
                ggplot2::aes(x = reorder(stringr::str_wrap(region, 5), -sum_length),
                             y = mean_gain, size = geneS),
                color = "gray12", show.legend = TRUE
            ) +
            ggplot2::geom_segment(
                ggplot2::aes(x = reorder(stringr::str_wrap(region, 5), -sum_length),
                             y = 0,
                             xend = reorder(stringr::str_wrap(region, 5), sum_length),
                             yend = (set * 3) + 0.00000001),
                linetype = "dashed", color = "gray12"
            ) +
            ggplot2::coord_polar(clip = "off") +
            ggplot2::annotate(x = 0, y = (set * 1), label = set, geom = "text", color = "gray12", size = 4, family = "Bell MT") +
            ggplot2::annotate(x = 0, y = (set * 2), label = set * 2, geom = "text", color = "gray12", size = 4, family = "Bell MT") +
            ggplot2::annotate(x = 0, y = (set * 3), label = set * 3, geom = "text", color = "gray12", size = 4, family = "Bell MT") +
            ggplot2::scale_y_continuous(limits = c(-set * 0.8, (set * 3) + 1),
                                        expand = c(0, 0), breaks = NULL) +
            ggplot2::scale_fill_gradientn("Number of terms", colours = colorset) +
            ggplot2::guides(fill = ggplot2::guide_colorsteps(barwidth = 15, barheight = 0.5,
                                                             title.position = "top", title.hjust = 0.5)) +
            ggplot2::geom_text(
                ggplot2::aes(x = reorder(stringr::str_wrap(region, 5), -sum_length),
                             y = -set * 0.2, label = tree),
                colour = fontcol, alpha = 1, size = 3, fontface = "bold", inherit.aes = FALSE
            ) +
            ggplot2::theme(
                axis.title = ggplot2::element_blank(),
                axis.ticks = ggplot2::element_blank(),
                axis.text.y = ggplot2::element_blank(),
                axis.ticks.y = ggplot2::element_blank(),
                axis.text.x = ggplot2::element_text(color = "gray12", size = 15, face = "bold", family = "Bell MT"),
                text = ggplot2::element_text(color = "gray12", family = "Bell MT"),
                panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                panel.grid = ggplot2::element_blank(),
                panel.grid.major.x = ggplot2::element_blank(),
                legend.position = "none"
            )
    }

    options(repr.plot.width = width, repr.plot.height = height,
            repr.plot.res = 300, repr.plot.pointsize = 10)

    ppp
}
