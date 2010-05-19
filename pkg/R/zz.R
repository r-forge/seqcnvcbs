.First.lib <- function(lib,pkg) {
   		library.dynam("SeqCNVCBS",pkg,lib)
   		cat("Dynamic library loaded\n")
   		require(clue)
	}
