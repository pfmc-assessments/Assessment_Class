#' Render Table of Contents
#' 
#' A simple function to extract headers from an xaringan RMarkdown
#' and build a table of contents. Returns a markdown list with links to the 
#' headers using the `name:` attribute of each slide.
#' 
#' @section Usage:
#' Just drop in a chunk where you want the toc to appear (set `echo=FALSE`):
#' 
#'     # Table of Contents
#' 
#'     ```{r echo=FALSE}
#'     render_toc("/path/to/the/file.Rmd")
#'     ```
#' 
#' @param filename Name of xaringan RMarkdown
render_toc <- function(filename) {
  x <- readLines(filename, warn = FALSE)
  x5 <- stringr::str_sub(x, 1, 5)
  xname <- stringr::str_which(x5, "name:")
  xtext <- stringr::str_which(x5, "text:")
  xtext <- stringr::str_trim(sapply(x[xtext], function(i){stringr::str_split(i, "text:")[[1]][2]}))
  has.xtext <- stringr::str_detect(x5[xname+1], "text:")
  header_slug <- stringr::str_trim(sapply(x[xname], function(i){stringr::str_split(i, "name:")[[1]][2]}))
  header_text <- header_slug
  header_text[has.xtext] <- xtext
  x <- paste0("* [", header_text, "](#", header_slug, ")")
  x <- c(".flexcolumn[",paste0("* [", header_text, "](#", header_slug, ")"), "]")
  knitr::asis_output(paste(x, collapse = "\n"))
}