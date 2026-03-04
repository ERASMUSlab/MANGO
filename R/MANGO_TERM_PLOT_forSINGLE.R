#' Plot MANGO terms for a single tree (bar or chord)
#'
#' Generate a bar plot (term scores) or a chord diagram (term–gene links)
#' for a selected tree from MANGO single-case outputs.
#'
#' @param filepath Working directory path.
#' @param DEG_list_name DEG list filename under `filepath`.
#' @param plot Plot type: "bar" or "cir".
#' @param trends One of "UP", "DOWN", "SIG".
#' @param clusternum Tree (cluster) index to plot.
#' @param width Plot width (for notebook display options).
#' @param height Plot height (for notebook display options).
#' @param LABEL Title label for the plot.
#' @param lfc_leng Max abs LFC used for chord color scaling.
#'
#' @return For `plot="bar"`, returns a ggplot object (and prints it).
#'   For `plot="cir"`, returns `NULL` invisibly (circos draws to device).
#' @export
MANGO_TERM_PLOT_forSINGLE <- function(filepath,
                                      DEG_list_name,
                                      plot,
                                      trends,
                                      clusternum,
                                      width,
                                      height,
                                      LABEL = "text for result",
                                      lfc_leng = 2) {

  plot_type <- as.character(plot)[1]

  if (plot_type == "bar") {

    blacklist <- c("regulation","positive","negative","a","the","of","gene","type","I","II","to","by")
    blacklist <- as.data.frame(blacklist)
    colnames(blacklist) <- "Description"

    DEG_list_path <- paste0(filepath, "/", DEG_list_name)
    DEG_list_pathDF <- as.data.frame(DEG_list_path)

    TXTbased <- nrow(as.data.frame(DEG_list_pathDF[grepl("txt", DEG_list_pathDF[,1]),]))
    CSVbased <- nrow(as.data.frame(DEG_list_pathDF[grepl("csv", DEG_list_pathDF[,1]),]))
    BEDbased <- nrow(as.data.frame(DEG_list_pathDF[grepl("bed", DEG_list_pathDF[,1]),]))

    if (TXTbased > 0) OUTPUTpath <- paste0(as.data.frame(strsplit(DEG_list_pathDF[,1], split = "_list.txt"))[1,1], "_GO_list.txt")
    if (CSVbased > 0) OUTPUTpath <- paste0(as.data.frame(strsplit(DEG_list_pathDF[,1], split = "_list.csv"))[1,1], "_GO_list.csv")
    if (BEDbased > 0) OUTPUTpath <- paste0(as.data.frame(strsplit(DEG_list_pathDF[,1], split = "_list.bed"))[1,1], "_GO_list.bed")
    GO_list_path <- OUTPUTpath

    if (TXTbased > 0) OUTPUTpath <- paste0(as.data.frame(strsplit(DEG_list_pathDF[,1], split = "input_DEG_list.txt"))[1,1], "MANGO_PREPROCESSING_list.txt")
    if (CSVbased > 0) OUTPUTpath <- paste0(as.data.frame(strsplit(DEG_list_pathDF[,1], split = "input_DEG_list.csv"))[1,1], "MANGO_PREPROCESSING_list.csv")
    if (BEDbased > 0) OUTPUTpath <- paste0(as.data.frame(strsplit(DEG_list_pathDF[,1], split = "input_DEG_list.bed"))[1,1], "MANGO_PREPROCESSING_list.bed")
    MANGO_PREPROCESSING_list_path <- OUTPUTpath
    MANGO_PREPROCESSING_list_pathDF <- as.data.frame(MANGO_PREPROCESSING_list_path)

    if (TXTbased > 0) OUTPUTpath <- paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.txt"))[1,1], "MANGO_SEPERATE.txt")
    if (CSVbased > 0) OUTPUTpath <- paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.csv"))[1,1], "MANGO_SEPERATE.csv")
    if (BEDbased > 0) OUTPUTpath <- paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.bed"))[1,1], "MANGO_SEPERATE.bed")
    MANGO_SEPERATE_path <- OUTPUTpath

    if (TXTbased > 0) OUTPUTpath <- paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.txt"))[1,1], "MANGO_SEPERATE_forMULTI_range.txt")
    if (CSVbased > 0) OUTPUTpath <- paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.csv"))[1,1], "MANGO_SEPERATE_forMULTI_range.csv")
    if (BEDbased > 0) OUTPUTpath <- paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.bed"))[1,1], "MANGO_SEPERATE_forMULTI_range.bed")
    MANGO_SEPERATE_forMULTI_range_path <- OUTPUTpath

    fileTABLE_path <- MANGO_PREPROCESSING_list_path
    GOTABLE_path <- GO_list_path

    fileTABLE <- utils::read.table(fileTABLE_path, sep="\t", header=TRUE, fill=TRUE)

    MANGO_SEPERATE_INPUT <- MANGO_SEPERATE_path
    SEPERATE_INPUT <- utils::read.table(MANGO_SEPERATE_INPUT, sep="\t", header=TRUE, fill=TRUE)
    SEPERATE_INPUT <- SEPERATE_INPUT[order(-SEPERATE_INPUT[,2]),]
    SEPERATE_INPUT <- stats::na.omit(SEPERATE_INPUT)

    fileTABLE <- utils::read.table(fileTABLE_path, sep="\t", header=TRUE, fill=TRUE)
    GOTABLE  <- utils::read.table(GOTABLE_path,  sep="\t", header=TRUE, fill=TRUE)

    i <- 1
    GO_RAW <- utils::read.delim(GOTABLE[i,1], sep="\t", header=TRUE, quote="", comment.char="",
                               stringsAsFactors=FALSE, fill=TRUE)
    GO_DF <- GO_RAW[,c(2,5,11)]
    colnames(GO_DF)[2] <- paste0("RF_", i)
    GO_DF[is.na(GO_DF)] <- 0

    i <- 1
    SEPERATE_INPUT_DF <- as.data.frame(strsplit(SEPERATE_INPUT[i,1], split="//"))
    colnames(SEPERATE_INPUT_DF) <- "Description"
    for (s in seq_len(nrow(SEPERATE_INPUT_DF))) {
      SEPERATE_INPUT_DF[s,1] <- gsub("_", " ", SEPERATE_INPUT_DF[s,1])
    }

    DATA_RAW_DF <- dplyr::inner_join(SEPERATE_INPUT_DF, GO_DF, by="Description")
    DATA_RAW_DF <- DATA_RAW_DF[order(-DATA_RAW_DF[,(2)]),]
    DATA_RAW_DF <- cbind(DATA_RAW_DF, i)
    colnames(DATA_RAW_DF) <- c("Description","RF","gene","cluster")

    if (nrow(SEPERATE_INPUT) >= 2) {
      for (i in 2:nrow(SEPERATE_INPUT)) {
        SEPERATE_INPUT_DF <- as.data.frame(strsplit(SEPERATE_INPUT[i,1], split="//"))
        colnames(SEPERATE_INPUT_DF) <- "Description"
        for (s in seq_len(nrow(SEPERATE_INPUT_DF))) {
          SEPERATE_INPUT_DF[s,1] <- gsub("_", " ", SEPERATE_INPUT_DF[s,1])
        }
        DATA_RAW_DF_frag <- dplyr::inner_join(SEPERATE_INPUT_DF, GO_DF, by="Description")
        DATA_RAW_DF_frag <- DATA_RAW_DF_frag[order(-DATA_RAW_DF_frag[,(2)]),]
        DATA_RAW_DF_frag <- cbind(DATA_RAW_DF_frag, i)
        colnames(DATA_RAW_DF_frag) <- c("Description","RF","gene","cluster")
        DATA_RAW_DF <- rbind(DATA_RAW_DF, DATA_RAW_DF_frag)
      }
    }

    DATA_RAW_DF[,5] <- 0
    colnames(DATA_RAW_DF)[5] <- "score"

    for (i in seq_len(nrow(DATA_RAW_DF))) {
      gene <- as.data.frame(strsplit(DATA_RAW_DF[i,3], split="/"))
      colnames(gene) <- "gene"

      DEG <- utils::read.table(DEG_list_pathDF[1,1], sep="\t", header=FALSE, fill=TRUE)
      DEG <- utils::read.table(DEG[1,1], sep="\t", header=TRUE, fill=TRUE)
      DEG[,3] <- 2^abs(DEG[,3])
      DEG <- dplyr::inner_join(DEG, gene, by="gene")

      DATA_RAW_DF[i,5] <- 1/(log10(DATA_RAW_DF[i,2])*(-1))*stats::median(DEG[,3])
    }

    DATA_RAW_DF <- DATA_RAW_DF[which(DATA_RAW_DF[,4] == clusternum),]
    DATA_RAW_DF <- DATA_RAW_DF[order(DATA_RAW_DF[,5]),]

    DATA_RAW_DF$generatio <- 0

    DEG <- utils::read.table(DEG_list_pathDF[1,1], sep="\t", header=FALSE, fill=TRUE)
    DEG <- utils::read.table(DEG[1,1], sep="\t", header=TRUE, fill=TRUE)

    for (i in seq_len(nrow(DATA_RAW_DF))) {
      DATA_RAW_DF[i,6] <- (nrow(as.data.frame(strsplit(DATA_RAW_DF[i,3], split="/")))/nrow(DEG)*10)
    }

    DEG <- utils::read.table(DEG_list_pathDF[1,1], sep="\t", header=FALSE, fill=TRUE)
    DEG <- utils::read.table(DEG[1,1], sep="\t", header=TRUE, fill=TRUE)
    DEG[,3] <- 2^abs(DEG[,3])
    DATA_RAW_DF$geneS <- 0

    for (i in seq_len(nrow(DATA_RAW_DF))) {
      geneSset <- as.data.frame(strsplit(DATA_RAW_DF[i,3], split="/"))
      colnames(geneSset) <- "gene"
      geneSset <- dplyr::inner_join(geneSset, DEG, by="gene")
      DATA_RAW_DF[i,7] <- stats::median(geneSset[,3])
    }

    DATA_DF <- DATA_RAW_DF[,c(1,5,6,7)]
    DATA_DF$cluster <- clusternum

    i <- 1
    GO_RAW <- utils::read.delim(GOTABLE[i,1], sep="\t", header=TRUE, quote="", comment.char="",
                               stringsAsFactors=FALSE, fill=TRUE)
    GO_DF <- GO_RAW[,c(2,9)]

    DATA_DF <- dplyr::inner_join(DATA_DF, GO_DF, by="Description")
    DATA_DF$P <- log10(DATA_DF[,6])*(-1)

    if (trends == "UP")  { fillcol <- "pink2";    fontcol <- "darkred"  }
    if (trends == "DOWN"){ fillcol <- "skyblue2"; fontcol <- "darkblue" }
    if (trends == "SIG") { fillcol <- "tan2";     fontcol <- "brown"    }

    if (trends == "UP")   colorset <- c("#F5B78EFF","#F19C7CFF","#EA8171FF","#DD686CFF","#CA5268FF","#B13F64FF")
    if (trends == "DOWN") colorset <- c("#CCFDFFFF","#99F8FFFF","#66F0FFFF","#33E4FFFF","#00AACCFF","#007A99FF")
    if (trends == "SIG")  colorset <- c("tan1","tan2","tan3","orange","brown")

    options(repr.plot.width = width, repr.plot.height = height, repr.plot.res = 500, repr.plot.pointsize = 10)

    plot_bar <- ggplot2::ggplot(DATA_DF, ggplot2::aes(score, forcats::fct_reorder(Description, score), fill = P)) +
      ggplot2::scale_fill_gradientn("Normalized\nAdjust P-value", colours = colorset) +
      ggplot2::geom_bar(stat="identity",
                        position=ggplot2::position_dodge(width = 0.8),
                        colour="black",
                        width=0.85,
                        linewidth=0.1) +
      ggplot2::geom_segment(ggplot2::aes(x = 0, y = Description, xend = generatio, yend = Description),
                            linetype = "dashed", color = "gray12") +
      ggplot2::geom_point(ggplot2::aes(x = generatio, y = Description, size = geneS),
                          color = "gray12", show.legend = TRUE) +
      ggplot2::scale_size_continuous("FC of genes", range = c(1, 5)) +
      theme_dose(12) +
      ggplot2::labs(title = paste0(LABEL),
                    subtitle = paste0("Tree ", clusternum)) +
      ggplot2::ylab("\n\n\n") +
      ggplot2::xlab("Tscore & Gene distribution ratio") +
      ggplot2::facet_grid(cluster ~ ., scales = "free_y", space = "free_y") +
      ggplot2::theme(strip.text.y = ggplot2::element_text(face="bold", size=14, color="darkblue"),
                     strip.background = ggplot2::element_rect(fill = fillcol, color = fontcol, linewidth = 0.6)) +
      ggplot2::theme(plot.title=ggplot2::element_text(face="bold", hjust=1, size=25),
                     plot.subtitle=ggplot2::element_text(face="bold", hjust=1, size=20),
                     axis.text.x=ggplot2::element_text(face="bold", size=7),
                     axis.text.y=ggplot2::element_text(face="bold", size=15),
                     axis.title.x=ggplot2::element_text(face="bold", size=15),
                     legend.title=ggplot2::element_text(face="bold", size=10),
                     legend.text=ggplot2::element_text(face="bold", size=8))

    print(plot_bar)
  }

  if (plot_type == "cir") {

    condition <- 1
    blacklist <- c("regulation","positive","negative","a","the","of","gene","type","I","II","to","by")
    blacklist <- as.data.frame(blacklist)
    colnames(blacklist) <- "Description"

    DEG_list_path <- paste0(filepath,"/",DEG_list_name)
    DEG_list_pathDF <- as.data.frame(DEG_list_path)

    TXTbased <- nrow(as.data.frame(DEG_list_pathDF[grepl("txt", DEG_list_pathDF[,1]),]))
    CSVbased <- nrow(as.data.frame(DEG_list_pathDF[grepl("csv", DEG_list_pathDF[,1]),]))
    BEDbased <- nrow(as.data.frame(DEG_list_pathDF[grepl("bed", DEG_list_pathDF[,1]),]))

    if (TXTbased > 0) OUTPUTpath <- paste0(as.data.frame(strsplit(DEG_list_pathDF[,1], split = "_list.txt"))[1,1],"_GO_list.txt")
    if (CSVbased > 0) OUTPUTpath <- paste0(as.data.frame(strsplit(DEG_list_pathDF[,1], split = "_list.csv"))[1,1],"_GO_list.csv")
    if (BEDbased > 0) OUTPUTpath <- paste0(as.data.frame(strsplit(DEG_list_pathDF[,1], split = "_list.bed"))[1,1],"_GO_list.bed")
    GO_list_path <- OUTPUTpath

    if (TXTbased > 0) OUTPUTpath <- paste0(as.data.frame(strsplit(DEG_list_pathDF[,1], split = "input_DEG_list.txt"))[1,1],"MANGO_PREPROCESSING_list.txt")
    if (CSVbased > 0) OUTPUTpath <- paste0(as.data.frame(strsplit(DEG_list_pathDF[,1], split = "input_DEG_list.csv"))[1,1],"MANGO_PREPROCESSING_list.csv")
    if (BEDbased > 0) OUTPUTpath <- paste0(as.data.frame(strsplit(DEG_list_pathDF[,1], split = "input_DEG_list.bed"))[1,1],"MANGO_PREPROCESSING_list.bed")
    MANGO_PREPROCESSING_list_path <- OUTPUTpath
    MANGO_PREPROCESSING_list_pathDF <- as.data.frame(MANGO_PREPROCESSING_list_path)

    if (TXTbased > 0) OUTPUTpath <- paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.txt"))[1,1],"MANGO_SEPERATE.txt")
    if (CSVbased > 0) OUTPUTpath <- paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.csv"))[1,1],"MANGO_SEPERATE.csv")
    if (BEDbased > 0) OUTPUTpath <- paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.bed"))[1,1],"MANGO_SEPERATE.bed")
    MANGO_SEPERATE_path <- OUTPUTpath

    fileTABLE_path <- MANGO_PREPROCESSING_list_path
    GOTABLE_path <- GO_list_path

    MANGO_SEPERATE_INPUT <- MANGO_SEPERATE_path
    SEPERATE_INPUT <- utils::read.table(MANGO_SEPERATE_INPUT, sep="\t", header=TRUE, fill=TRUE)
    SEPERATE_INPUT <- SEPERATE_INPUT[order(-SEPERATE_INPUT[,2]),]
    SEPERATE_INPUT <- stats::na.omit(SEPERATE_INPUT)

    GOTABLE <- utils::read.table(GOTABLE_path, sep="\t", header=TRUE, fill=TRUE)

    i <- 1
    GO_RAW <- utils::read.delim(GOTABLE[i,1], sep="\t", header=TRUE, quote="", comment.char="",
                               stringsAsFactors=FALSE, fill=TRUE)
    GO_DF <- GO_RAW[,c(2,5,11)]
    colnames(GO_DF)[2] <- paste0("RF_", i)
    GO_DF[is.na(GO_DF)] <- 0

    i <- 1
    SEPERATE_INPUT_DF <- as.data.frame(strsplit(SEPERATE_INPUT[i,1], split="//"))
    colnames(SEPERATE_INPUT_DF) <- "Description"
    for (s in seq_len(nrow(SEPERATE_INPUT_DF))) {
      SEPERATE_INPUT_DF[s,1] <- gsub("_", " ", SEPERATE_INPUT_DF[s,1])
    }

    DATA_RAW_DF <- dplyr::inner_join(SEPERATE_INPUT_DF, GO_DF, by="Description")
    DATA_RAW_DF <- DATA_RAW_DF[order(-DATA_RAW_DF[,(condition+1)]),]
    DATA_RAW_DF <- cbind(DATA_RAW_DF, i)
    colnames(DATA_RAW_DF) <- c("Description","RF","gene","cluster")

    if (nrow(SEPERATE_INPUT) >= 2) {
      for (i in 2:nrow(SEPERATE_INPUT)) {
        SEPERATE_INPUT_DF <- as.data.frame(strsplit(SEPERATE_INPUT[i,1], split="//"))
        colnames(SEPERATE_INPUT_DF) <- "Description"
        for (s in seq_len(nrow(SEPERATE_INPUT_DF))) {
          SEPERATE_INPUT_DF[s,1] <- gsub("_", " ", SEPERATE_INPUT_DF[s,1])
        }
        DATA_RAW_DF_frag <- dplyr::inner_join(SEPERATE_INPUT_DF, GO_DF, by="Description")
        DATA_RAW_DF_frag <- DATA_RAW_DF_frag[order(-DATA_RAW_DF_frag[,(condition+1)]),]
        DATA_RAW_DF_frag <- cbind(DATA_RAW_DF_frag, i)
        colnames(DATA_RAW_DF_frag) <- c("Description","RF","gene","cluster")
        DATA_RAW_DF <- rbind(DATA_RAW_DF, DATA_RAW_DF_frag)
      }
    }

    i <- 1
    GO_RAW <- utils::read.delim(GOTABLE[i,1], sep="\t", header=TRUE, quote="", comment.char="",
                               stringsAsFactors=FALSE, fill=TRUE)
    GO_DF2 <- GO_RAW[,c(2,9)]
    DATA_RAW_DF <- dplyr::inner_join(DATA_RAW_DF, GO_DF2, by="Description")
    DATA_RAW_DF <- DATA_RAW_DF[order(DATA_RAW_DF[,ncol(DATA_RAW_DF)]),]
    DATA_RAW_DF <- DATA_RAW_DF[,-ncol(DATA_RAW_DF)]

    DATA_RAW_DF[,5] <- 0
    colnames(DATA_RAW_DF)[5] <- "score"

    for (i in seq_len(nrow(DATA_RAW_DF))) {
      gene <- as.data.frame(strsplit(DATA_RAW_DF[i,3], split="/"))
      colnames(gene) <- "gene"

      DEG <- utils::read.table(DEG_list_pathDF[1,1], sep="\t", header=FALSE, fill=TRUE)
      DEG <- utils::read.table(DEG[1,1], sep="\t", header=TRUE, fill=TRUE)

      DEG <- dplyr::inner_join(DEG, gene, by="gene")
      DATA_RAW_DF[i,5] <- 1/(log10(DATA_RAW_DF[i,2])*(-1))*stats::median(DEG[,3])
    }

    DATA_RAW_DF <- DATA_RAW_DF[which(DATA_RAW_DF[,4] == clusternum),]
    DATA_RAW_DF <- DATA_RAW_DF[order(DATA_RAW_DF[,5]),]

    treeterm <- as.data.frame(DATA_RAW_DF[,1])
    treegene_tmp <- as.data.frame(DATA_RAW_DF[,3])

    treegene <- data.frame(gene = character())
    for (i in seq_len(nrow(treegene_tmp))) {
      treegene_frag <- as.data.frame(strsplit(treegene_tmp[i,1], split="/"))
      colnames(treegene_frag) <- "gene"
      treegene <- rbind(treegene, treegene_frag)
      treegene <- unique(treegene)
    }

    chordDF <- treeterm
    colnames(chordDF) <- "Description"

    for (i in seq_len(nrow(treegene))) {
      block <- DATA_RAW_DF[grepl(treegene[i,1], DATA_RAW_DF[,3]), 1]
      block <- as.data.frame(block)
      colnames(block) <- "Description"
      block$label <- 1
      colnames(block)[2] <- treegene[i,1]

      chordDF <- dplyr::full_join(chordDF, block, by="Description")
      chordDF[is.na(chordDF)] <- 0
    }

    chordDF <- t(chordDF)
    chordDF <- as.data.frame(chordDF)
    colnames(chordDF) <- chordDF[1,]
    chordDF <- chordDF[-1,]

    for (i in seq_len(ncol(chordDF))) chordDF[,i] <- as.numeric(chordDF[,i])

    chordDF$hit <- 0
    for (i in seq_len(nrow(chordDF))) chordDF[i,ncol(chordDF)] <- sum(as.data.frame(t(chordDF[i,]))[,1])

    chordDF <- chordDF[order(chordDF$hit),]
    chordDF <- chordDF[,-ncol(chordDF)]

    chordDF <- cbind(rownames(chordDF), chordDF)
    colnames(chordDF)[1] <- "gene"

    DEG <- utils::read.table(DEG_list_pathDF[1,1], sep="\t", header=FALSE, fill=TRUE)
    DEG <- utils::read.table(DEG[1,1], sep="\t", header=TRUE, fill=TRUE)[,c(1,3)]

    chordDF <- dplyr::inner_join(chordDF, DEG, by="gene")
    chordDF <- as.matrix(chordDF)

    rownames(chordDF) <- chordDF[,1]
    chordDF <- chordDF[,-1]

    chordDF_num <- apply(chordDF, 2, as.numeric)
    rownames(chordDF_num) <- rownames(chordDF)
    colnames(chordDF_num) <- colnames(chordDF)

    if (trends == "UP") {
      mycol <- c("#F5B78EFF","#F3AE88FF","#F1A482FF","#EF9A7CFF","#ED9077FF",
                 "#EA8672FF","#E67C70FF","#E2726DFF","#DE696CFF","#D9606BFF",
                 "#D4576AFF","#CF4F69FF","#C94768FF","#C33F67FF","#BC3866FF",
                 "#B53165FF","#AE2B64FF","#A72563FF","#A02062FF","#991B61FF")
    }
    if (trends == "DOWN") {
      mycol <- c("#CCFDFFFF","#BFFBFFFF","#A5F8FFFF","#8AF3FFFF","#70EEFFFF",
                 "#56E7FFFF","#3DDAF9FF","#26CAE6FF","#10B9D3FF","#009FBCFF","#007A99FF")
    }

    options(repr.plot.width = width, repr.plot.height = height, repr.plot.res = 500, repr.plot.pointsize = 10)

    ch <- as.data.frame(chordDF_num)
    go_terms <- setdiff(colnames(ch), c("log2FoldChange","logFC","LFC","lfc"))
    genes <- rownames(ch)
    mat <- as.matrix(ch[, go_terms, drop = FALSE])
    mat[is.na(mat)] <- 0

    go_cols <- setNames(rep(mycol, length.out = length(go_terms)), go_terms)
    lfc <- ch$log2FoldChange
    names(lfc) <- genes
    lfc_max <- lfc_leng
    lfc2 <- pmax(pmin(lfc, lfc_max), -lfc_max)

    if (trends == "UP")   col_fun <- circlize::colorRamp2(c(0, lfc_max), c("white","pink3"))
    if (trends == "DOWN") col_fun <- circlize::colorRamp2(c(0, lfc_max), c("white","skyblue3"))

    gene_cols <- setNames(col_fun(lfc2), genes)

    grid_col <- c(go_cols, gene_cols)
    link_cols <- matrix(NA, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    for (g in rownames(mat)) {
      for (t in colnames(mat)) {
        if (mat[g, t] > 0) link_cols[g, t] <- go_cols[[t]]
      }
    }

    circlize::circos.clear()
    circlize::circos.par(start.degree = 270)

    circlize::chordDiagram(x = mat,
                           grid.col = grid_col,
                           col = link_cols,
                           transparency = 0.1,
                           annotationTrack = "grid",
                           preAllocateTracks = list(track.height = 0.14))

    circlize::circos.trackPlotRegion(track.index = 1,
                                     bg.border = NA,
                                     panel.fun = function(x, y) {
                                       sector <- circlize::get.cell.meta.data("sector.index")
                                       xlim <- circlize::get.cell.meta.data("xlim")
                                       ylim <- circlize::get.cell.meta.data("ylim")

                                       if (sector %in% go_terms) {
                                         circlize::circos.text(x = mean(xlim),
                                                               y = ylim[1] + circlize::mm_y(2),
                                                               labels = sector,
                                                               facing = "downward",
                                                               niceFacing = TRUE,
                                                               adj = c(0, 0.5),
                                                               cex = 1,
                                                               font = 2)
                                       } else if (sector %in% genes) {
                                         circlize::circos.text(x = mean(xlim),
                                                               y = ylim[1.5] + circlize::mm_y(1),
                                                               labels = sector,
                                                               facing = "clockwise",
                                                               niceFacing = TRUE,
                                                               adj = c(0, 0.5),
                                                               cex = 0.6,
                                                               font = 2)
                                       }
                                     })

    return(invisible(NULL))
  }
}
