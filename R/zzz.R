.First.lib <- function(lib, pkg)
{
  library.dynam("MixSim", pkg, lib)

  require(MASS)

}

.Last.lib <- function(lib, pkg)
{
  library.dynam.unload("MixSim", pkg, lib)

  load.lib <- c("MASS")

  for(i in load.lib){
    pos <- match(paste("package:", i, sep = ""), search())
    if(! is.na(pos)){
      detach(pos = pos)
    }
  }
}

