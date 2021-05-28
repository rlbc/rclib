#' Import Brohl et al 2014 files
#'
#' This function import my files.
#' @keywords import
#' @export
#' @examples
#' import_files()
import_files <- function() {
  # Brohl (final _b)
  # Deixarei assim caso precise adicionar novos dados
  samples_b <- utils::read.csv(file.path("samples_brohl.csv"), header = TRUE, sep = ',')
  #samples_b <- samples_b %>% mutate(Study = "Brohl")
  files_b <- base::file.path("quantification/quants", samples_b$Directory, "quant_correct.sf")
  names(files_b) <- paste0(samples_b$Run)
}

#' RLE function
#'
#' This function calculates RLE and plots.
#' @param expTable Expression table.
#' @param logTable If expression table is already on log scale. Defaults to FALSE.
#' @param plotTable Plot results. Defaults to FALSE.
#' @keywords RLE
#' @export
#' @examples
#' calcrle()
calcrle <- function(expTable, logTable = FALSE, plotTable = FALSE) {
  if (logTable == FALSE) {
    expTable <- log2(expTable + 1)
  }
  
  #Calculate median value for rows (genes)
  med <- expTable %>% dplyr::rowwise() %>% dplyr::summarise(median = stats::median(c_across(V1:V12)))
  
  #Maps and calculate log2 - median to a new table
  rle <- purrr::map2_df(expTable, med, `-`)
  
  #Plots if TRUE
  if (plotTable == TRUE){
    rle <- rle %>% dplyr::mutate(gene = rownames(tbl_log))
    
    #Long format
    rle_long <- rle %>% dplyr::pivot_longer(!gene, names_to = "sample")
    
    ggplot2::ggplot(rle_long, aes(x=sample, y=value, fill=sample)) + 
      geom_boxplot(outlier.color = NA) + 
      coord_cartesian(ylim = c(-5, 2)) + 
      scale_x_discrete(limits = colnames(tbl_log)) + 
      geom_hline(yintercept = 0, linetype="dashed")
  }
}