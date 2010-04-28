readSeqJShen <-
function(filename) {
	## We do not know if the two reads have chromosome match; one read missing
	seqMatch = read.table(filename, as.is=T)
	seqChr = factor(seqMatch[,1])
	seqPos = seqMatch[,2]
	rm(seqMatch)
	return(list(seqF=seqPos, seqChr=seqChr))
}

