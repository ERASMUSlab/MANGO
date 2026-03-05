#' Run the MANGO analysis pipeline
#'
#' Orchestrates preprocessing (optional), term listing, comparison, and separation.
#' Supports both single- and multi-condition workflows.
#'
#' @param filepath Working directory path.
#' @param type DEG calling mode. One of \code{"standard"} or \code{"broad"}.
#'   \code{"standard"} uses \code{adjp=0.05} and \code{FC=2}; \code{"broad"} uses
#'   \code{adjp=0.05} and \code{FC=1.5}.
#' @param full_condition Character vector of condition names (must match sample condition labels).
#' @param number_of_rep Integer vector of replicate counts for each condition in \code{full_condition}.
#' @param mango_design Character vector specifying contrasts and direction to export.
#' @param DEG_list_name DEG list filename under \code{filepath}.
#' @param ref_genome Reference genome code (\code{"mm"} or \code{"hs"}).
#' @param core Number of cores.
#' @param PASSED_RATIO Passed dependency ratio threshold.
#' @param PASSED_NUM Passed dependency count threshold.
#' @param similarity Similarity cutoff (%).
#' @param FC Fold-change ratio cutoff.
#' @param condition condition index for dynamic analysis.
#' @param dynamic_analyisis \code{"T"} or \code{"F"} to enable multi-range filtering.
#' @param preprocessing \code{"T"} or \code{"F"} to run preprocessing inside this function.
#'
#' @return Writes output files to disk. Returns \code{NULL} invisibly.
#'
#' @examples
#' stopifnot(is.function(MANGO_ANALYSIS))
#' \dontrun{
#' MANGO_ANALYSIS(
#'   filepath=".",
#'   DEG_list_name="input_DEG_list.txt",
#'   ref_genome="mm",
#'   core=1,
#'   PASSED_RATIO=15,
#'   PASSED_NUM=4,
#'   similarity=70,
#'   type = "broad",
#'   full_condition = c("DAY0","DAY4","DAY7","DAY10","DAY14","DAY21"),
#'   number_of_rep = c(3,3,3,6,3,3),
#'   mango_design = c("DAY4_DAY0_UP","DAY7_DAY0_UP","DAY10_DAY0_UP","DAY14_DAY0_UP","DAY21_DAY0_UP")
#' )
#' }
#' @export


MANGO_ANALYSIS = function(filepath,
                          DEG_list_name,
                          ref_genome,
                          core = 1,
                          PASSED_RATIO,
                          PASSED_NUM,
                          similarity,
			  type = "standard",
                          full_condition,
                          number_of_rep,
                          mango_design,
                          FC = 2,
                          condition = 1,
			  dynamic_analyisis = "F",
                          preprocessing = "F"){

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

    if(TXTbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.txt"))[1,1],"MANGO_TERMLISTING.txt")}
    if(CSVbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.csv"))[1,1],"MANGO_TERMLISTING.csv")}
    if(BEDbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.bed"))[1,1],"MANGO_TERMLISTING.bed")}
    MANGO_TERMLISTING_path = OUTPUTpath

    if(TXTbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.txt"))[1,1],"MANGO_COMPARE.txt")}
    if(CSVbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.csv"))[1,1],"MANGO_COMPARE.csv")}
    if(BEDbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.bed"))[1,1],"MANGO_COMPARE.bed")}
    MANGO_COMPARE_path = OUTPUTpath

    if(TXTbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.txt"))[1,1],"MANGO_SEPERATE.txt")}
    if(CSVbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.csv"))[1,1],"MANGO_SEPERATE.csv")}
    if(BEDbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.bed"))[1,1],"MANGO_SEPERATE.bed")}
    MANGO_SEPERATE_path = OUTPUTpath

    if(TXTbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.txt"))[1,1],"MANGO_SEPERATE_forMULTI_range.txt")}
    if(CSVbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.csv"))[1,1],"MANGO_SEPERATE_forMULTI_range.csv")}
    if(BEDbased>0){OUTPUTpath = paste0(as.data.frame(strsplit(MANGO_PREPROCESSING_list_pathDF[,1], split = "MANGO_PREPROCESSING_list.bed"))[1,1],"MANGO_SEPERATE_forMULTI_range.bed")}
    MANGO_SEPERATE_forMULTI_range_path = OUTPUTpath

    if(preprocessing == "T"){

        MANGO_DEGcalling(filepath = filepath,
                         type = type,
                         full_condition = full_condition,
                         number_of_rep = number_of_rep,
                         mango_design = mango_design)

        MANGO_PREPROCESSING(DEG_list_path = DEG_list_path,
                            GO_list_path = GO_list_path,
                            MANGO_PREPROCESSING_list_path = MANGO_PREPROCESSING_list_path,
                            filepath = filepath,
                            ref_genome = ref_genome,
                            core = core)
        }

    fileTABLE = utils::read.table(MANGO_PREPROCESSING_list_path, sep="\t", header=T, fill=TRUE)

    if(nrow(fileTABLE) > 1){

    MANGO_TERMLISTING(fileTABLE_path = MANGO_PREPROCESSING_list_path,
                      outputPATH = MANGO_TERMLISTING_path,
                      PASSED_RATIO = PASSED_RATIO,
                      PASSED_NUM = PASSED_NUM,
                      similarity = similarity)

    MANGO_COMPARE(input_TERMLISTING = MANGO_TERMLISTING_path,
                  fileTABLE_path = MANGO_PREPROCESSING_list_path,
                  filepath = filepath,
                  DEG_list_name = DEG_list_name,
                  GOTABLE_path = GO_list_path,
                  outputPATH = MANGO_COMPARE_path)

    MANGO_SEPERATE(fileTABLE_path = MANGO_PREPROCESSING_list_path,
               GOTABLE_path  = GO_list_path,
               input_COMPARE = MANGO_COMPARE_path,
               outputPATH    = MANGO_SEPERATE_path,
               FC            = FC,
               PASSED_NUM    = PASSED_NUM)

    if(dynamic_analyisis == "T"){
	MANGO_SEPERATE_forDA(input_COMPARE = MANGO_COMPARE_path,
	                     fileTABLE_path = MANGO_PREPROCESSING_list_path,
                             FC = FC,
                             outputPATH = MANGO_SEPERATE_forMULTI_range_path,
                             condition = condition)
        }
        }


    if(nrow(fileTABLE) == 1){

        MANGO_TERMLISTING_single(PASSED_RATIO = PASSED_RATIO,
                                 PASSED_NUM = PASSED_NUM,
                                 similarity = similarity,
                                 fileTABLE_path = MANGO_PREPROCESSING_list_path,
                                 outputPATH = MANGO_TERMLISTING_path)
        
        MANGO_COMPARE_SINGLE(input_TERMLISTING = MANGO_TERMLISTING_path,
                             fileTABLE_path = MANGO_PREPROCESSING_list_path,
                             filepath = filepath,
                             DEG_list_name = DEG_list_name,
                             GOTABLE_path = GO_list_path,
			     PASSED_NUM = PASSED_NUM,
                             outputPATH = MANGO_SEPERATE_path)
        }

}

