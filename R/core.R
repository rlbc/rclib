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
    expTable <- as.data.frame(log2(expTable + 1))
  }
  
  #Calculate median value for rows (genes)
  med <- expTable %>% dplyr::rowwise() %>% dplyr::summarise(median = stats::median(c_across()))
  
  #Maps and calculate log2 - median to a new table
  rle <- purrr::map2_df(expTable, med, `-`)
  
  #Plots if TRUE
  if (plotTable == TRUE){
    rle <- rle %>% dplyr::mutate(gene = rownames(expTable))
    
    #Long format
    rle_long <- rle %>% tidyr::pivot_longer(!gene, names_to = "sample")
    
    ggplot2::ggplot(rle_long, aes(x=sample, y=value, fill=sample)) + 
      geom_boxplot(outlier.color = NA) + 
      coord_cartesian(ylim = c(-5, 2)) + 
      scale_x_discrete(limits = colnames(expTable)) + 
      geom_hline(yintercept = 0, linetype="dashed") +
      theme(axis.text.x = element_text(angle=45, hjust=1, size=rel(0.9)))
  }
}

#' norm_quantiles function
#'
#' This function normalizes by quantile.
#' @param expTable Expression table.
#' @param logTable If expression table is already on log scale. Defaults to FALSE.
#' @keywords normalization
#' @export
#' @examples
#' norm_quantiles(expression)
norm_quantiles <- function(expTable, logTable = FALSE) {
  if (logTable) {
    exp_norm <- preprocessCore::normalize.quantiles(log2(expTable))
  }
  else {
    exp_norm <- preprocessCore::normalize.quantiles(log2(expTable + 1))
  }
  return(exp_norm)
}

#' norm_deseq function
#'
#' This function normalizes by DESeq2.
#' @param txiObject txi object.
#' @keywords normalization
#' @export
#' @examples
#' norm_deseq(expression)
norm_deseq <- function(txiObject) {
  sampleTable <- data.frame(matrix("ewing", nrow = 26))
  rownames(sampleTable) <- colnames(txiObject$counts)
  
  ddsTxi <- DESeq2::DESeqDataSetFromTximport(txiObject, colData = sampleTable, design = ~ 1)
  ddsTxi <- DESeq2::estimateSizeFactors(ddsTxi)
  norm_dds <- counts(ddsTxi, normalized = TRUE)
  
  return(norm_dds)
}

#' calc_var function
#'
#' This function calculates the various variance and plots.
#' @param expTable Table with expression values.
#' @param logTable If expression table is already on log scale. Defaults to FALSE.
#' @keywords Variance
#' @export
#' @examples
#' calc_var(expression)
calc_var <- function(expTable, logTable = FALSE) {
  if (!logTable) {
    expTable <- log2(expTable + 1)
  }
  
  # Transpose table
  t_element <- as.data.frame(t(expTable))
    
  # Calculate variance
  v_element <- as.data.frame(sapply(t_element, stats::var))
  
  return(v_element)
}

#' calc_mean function
#'
#' This function calculates the various means.
#' @param expTable Table with expression values.
#' @param logTable If expression table is already on log scale. Defaults to FALSE.
#' @keywords Mean
#' @export
#' @examples
#' calc_mean(expression)
calc_mean <- function(expTable, logTable = FALSE) {
  if (!logTable) {
    expTable <- log2(expTable + 1)
  }
  
  # Transpose table
  t_element <- as.data.frame(t(expTable))
  
  # Calculate variance
  v_element <- as.data.frame(sapply(t_element, mean))
  
  colnames(v_element) <- "mean"
  
  return(v_element)
}

#' plotGenes function
#'
#' This function plots various genes.
#' @param expTable Expression table.
#' @param geneList List of genes to plot.
#' @keywords plot
#' @export
#' @examples
#' plotGenes(expression)
plotGenes <- function(expTable, geneList) {
  
}