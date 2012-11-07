.First.lib <- function(lib, pkg){
	library.dynam("MixSim", pkg, lib)
}

.Last.lib <- function(lib, pkg){
	library.dynam.unload("MixSim", pkg, lib)
}

