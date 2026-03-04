#' Run GO enrichment preprocessing for MANGO
#'
#' Reads DEG lists, runs \code{clusterProfiler::enrichGO} (BP), writes per-case GO
#' tables, then calls the preprocessing shell script and writes the resulting
#' preprocessing list.
#'
#' @param DEG_list_path Path to a table listing DEG files (one file per row).
#' @param GO_list_path Output path to write the list of generated GO result files.
#' @param MANGO_PREPROCESSING_list_path Output path to write the list of generated
#'   MANGO preprocessing result files.
#' @param ref_genome Reference genome code. Use \code{"mm"} for mouse or \code{"hs"} for human.
#' @param core Number of cores to pass to the preprocessing shell script.
#' @param filepath Base directory where \code{CODE/init_MANGO_PREPROCESSING.sh} exists.
#'
#' @return Writes output files to disk. Returns \code{NULL} invisibly.
#'
#' @examples
#' stopifnot(is.function(MANGO_PREPROCESSING))
#' \dontrun{
#' MANGO_PREPROCESSING(
#'   DEG_list_path = "input_DEG_list.txt",
#'   GO_list_path = "input_DEG_GO_list.txt",
#'   MANGO_PREPROCESSING_list_path = "MANGO_PREPROCESSING_list.txt",
#'   ref_genome = "mm",
#'   core = 1,
#'   filepath = "."
#' )
#' }
#' @export



MANGO_PREPROCESSING = function(DEG_list_path,
                               GO_list_path,
                               MANGO_PREPROCESSING_list_path,
                               ref_genome,
                               core,
                               filepath){

    input_DEG_list = utils::read.table(DEG_list_path,sep = "\t", header=F,fill = TRUE)

    GO_list = "GO"
    GO_list = as.data.frame(GO_list)
    colnames(GO_list) = "list"
    GO_list = GO_list[-1,]

    for(i in 1:nrow(input_DEG_list)){
        input_DEG_list_sep = input_DEG_list[i,1]
        input_DEG_list_sep_DF = as.data.frame(input_DEG_list_sep)

        TXTbased = nrow(as.data.frame(input_DEG_list_sep_DF[grepl("txt", input_DEG_list_sep_DF[,1]),]))
        CSVbased = nrow(as.data.frame(input_DEG_list_sep_DF[grepl("csv", input_DEG_list_sep_DF[,1]),]))
        BEDbased = nrow(as.data.frame(input_DEG_list_sep_DF[grepl("bed", input_DEG_list_sep_DF[,1]),]))
        
        if(TXTbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(input_DEG_list_sep_DF[,1], split = ".txt"))[1,1],"_GO.txt")}
        if(CSVbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(input_DEG_list_sep_DF[,1], split = ".csv"))[1,1],"_GO.csv")}
        if(BEDbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(input_DEG_list_sep_DF[,1], split = ".bed"))[1,1],"_GO.bed")}
        
        OUTPUTpath_DF = as.data.frame(OUTPUTpath)
        colnames(OUTPUTpath_DF) = "list"
        
        GO_list = rbind(GO_list,OUTPUTpath_DF)
        
        if(ref_genome == "mm"){
            if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
              stop("Package 'org.Mm.eg.db' is required for ref_genome == 'mm'. Please install it.")
            }

            DF = utils::read.table(input_DEG_list_sep,sep = "\t", header=T,fill = TRUE)
    
            ego = clusterProfiler::enrichGO(gene          = DF[,1],
                                            OrgDb         = org.Mm.eg.db::org.Mm.eg.db,
                                            keyType       = 'SYMBOL',
                                            ont           = "BP")

            # 최소 안전화: enrichResult에 dplyr::mutate를 바로 적용하지 말고 df로 변환 후 적용
            ego_df = as.data.frame(ego)

            egoplotdata1 = dplyr::mutate(ego_df, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
            egoplotdata2 = dplyr::mutate(egoplotdata1, BgRatio_num = as.numeric(sub("/\\d+", "", BgRatio))/1)
            egoplotdata  = dplyr::mutate(egoplotdata2, GeneRatio_num = as.numeric(sub("/\\d+", "", GeneRatio))/1)
            egoplotdata  = as.data.frame(egoplotdata)

            utils::write.table(egoplotdata,
                        file = OUTPUTpath,
                        col.names=T, row.names=F, quote=F,sep="\t")
            }

        if(ref_genome == "hs"){
            if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
              stop("Package 'org.Hs.eg.db' is required for ref_genome == 'hs'. Please install it.")
            }

            DF = utils::read.table(input_DEG_list_sep,sep = "\t", header=T,fill = TRUE)
    
            ego = clusterProfiler::enrichGO(gene          = DF[,1],
                                            OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
                                            keyType       = 'SYMBOL',
                                            ont           = "BP")

            # 최소 안전화: df로 변환 후 mutate
            ego_df = as.data.frame(ego)

            egoplotdata1 = dplyr::mutate(ego_df, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
            egoplotdata2 = dplyr::mutate(egoplotdata1, BgRatio_num = as.numeric(sub("/\\d+", "", BgRatio))/1)
            egoplotdata  = dplyr::mutate(egoplotdata2, GeneRatio_num = as.numeric(sub("/\\d+", "", GeneRatio))/1)
            egoplotdata  = as.data.frame(egoplotdata)

            utils::write.table(egoplotdata,
                        file = OUTPUTpath,
                        col.names=T, row.names=F, quote=F,sep="\t")
            }
        }
    
    utils::write.table(GO_list,
                file = GO_list_path,
                col.names=T, row.names=F, quote=F,sep="\t")

    process_num = nrow(input_DEG_list)

    MANGO_PREPROCESSING_CMD = paste0("bash ",
                                     filepath,
                                     "/CODE/init_MANGO_PREPROCESSING.sh",
                                     " ",GO_list_path,
                                     " ",process_num,
                                     " ",core)

    system(MANGO_PREPROCESSING_CMD, intern = T)
    
    MANGO_PREPROCESSING_list = GO_list
    
    for(i in 1:nrow(GO_list)){
        sublist = GO_list[i,1]
        sublist = as.data.frame(sublist)

        TXTbased = nrow(as.data.frame(sublist[grepl("txt", sublist[,1]),]))
        CSVbased = nrow(as.data.frame(sublist[grepl("csv", sublist[,1]),]))
        BEDbased = nrow(as.data.frame(sublist[grepl("bed", sublist[,1]),]))

        if(TXTbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(sublist[,1], split = ".txt"))[1,1],"_MANGO_PREPROCESSING.txt")}
        if(CSVbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(sublist[,1], split = ".csv"))[1,1],"_MANGO_PREPROCESSING.csv")}
        if(BEDbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(sublist[,1], split = ".bed"))[1,1],"_MANGO_PREPROCESSING.bed")}

        MANGO_PREPROCESSING_list[i,1] = OUTPUTpath
        }

    utils::write.table(MANGO_PREPROCESSING_list,
                file = MANGO_PREPROCESSING_list_path,
                col.names=T, row.names=F, quote=F,sep="\t")
    }
