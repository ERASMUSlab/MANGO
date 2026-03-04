#' MANGO default plotting theme
#'
#' A ggplot2 theme used across MANGO plotting functions (grid-on style).
#' @param base_size Base font size used by the theme.
#' @return A ggplot2 theme object.
#' @export
theme_dose <- function(base_size = 12){
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      plot.title = ggplot2::element_text(face = 'bold'),
      legend.title = ggplot2::element_text(face = 'bold')
    )
}

