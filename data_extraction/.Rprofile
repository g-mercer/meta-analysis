source("renv/activate.R")

local({
  r <- getOption("repos")
  r["CRAN"] <- "https://cran.ma.imperial.ac.uk/"
  options(repos = r)
})

clr <- function() {
  ENV <- globalenv()
  ll <- ls(envir = ENV)
  ll <- ll[ll != "clr"]
  rm(list = ll, envir = ENV)
}
