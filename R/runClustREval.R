# runClustREval
# Author: Xin Zhi Fang (xinzhi.fang@mail.utoronto.ca)
# Date: December 8, 2021

#' Launch Shiny App for clustREval
#'
#' A function that launches the Shiny app for clustREval
#' The code has been placed in \code{./inst/shiny-scripts}.
#'
#' @return 
#'
#' @export
#' @importFrom shiny runApp

runClustREval <- function() {
  appDir <- system.file("shiny-scripts",
                        package = "clustREval")
  shiny::runApp(appDir, display.mode = "normal")
  return()
}
# [END]