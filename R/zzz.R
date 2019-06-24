.onAttach <- function(...) {
  #  pkgdesc <- utils::packageDescription("eggCounts", lib.loc = eggCountsLib)
  #  builddate <- gsub(';.*$', '', pkgdesc$Packaged)
  #  packageStartupMessage(paste("eggCounts (Version ", pkgdesc$Version, ")", sep = ""))
  packageStartupMessage(paste0("Loading spatialfusion (version ", utils::packageVersion("spatialfusion"),"):
- The compilation time for a Stan model can be up to 20s.
- We recommend using INLA method for larger datasets (several thousand observations).
- It is good practice to test your model on sub-sampled dataset first."))
}
