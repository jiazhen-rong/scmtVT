#' Filter out variants that are purely noise
#'
#' @param x A character vector with one element.
#' @param split What to split on.
#'
#' @return A character vector.
#' @export
#'
#' @examples
#' x <- "alfa,bravo,charlie,delta"
#' strsplit1(x, split = ",")
FilterVariants <- function(x, split) {
  strsplit(x, split = split)[[1]]
}
