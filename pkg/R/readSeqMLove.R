readSeqMLove <-
function(filename) {
	## Assume all matched reads from output
	seqMatch = read.table(filename, as.is=T)
	seqChr = factor(seqMatch[,1])
	seqFront = seqMatch[,3]
	normalLFflag = (seqFront=="F")
	seqPosF = rep(0, dim(seqMatch)[1])
	seqPosF[normalLFflag] = as.numeric(seqMatch[normalLFflag,2])
	seqPosF[-normalLFflag] = as.numeric(seqMatch[-normalLFflag,4])
	seqPosR = rep(0, dim(seqMatch)[1])
	seqPosR[normalLFflag] = as.numeric(seqMatch[normalLFflag,4])
	seqPosR[-normalLFflag] = as.numeric(seqMatch[-normalLFflag,2])
	rm(seqMatch)
	return(list(seqF=seqPosF, seqR=seqPosR, seqChr=seqChr))
}

