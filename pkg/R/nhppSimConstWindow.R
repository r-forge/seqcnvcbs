nhppSimConstWindow <-
function(controlRates, nSpike=25, cptLen=seq(5, 50, by=5), nRepeat=5, minGain=1.5, maxGain=5, minLoss=0.01, maxLoss=0.5, pGain=0.6) {
	nCptLen = length(cptLen)
	controlSim = nhppSimulate(controlRates)
	counter = 1
	simCBSRes = vector("list", nCptLen*nRepeat)
	for(i in 1:nCptLen) {
		for(j in 1:nRepeat) {
			newSpikeRate = nhppSpikeConstWindow(controlRates, nSpike, cptLen[i], minGain, maxGain, minLoss, maxLoss, pGain)
			newSpikeSim = nhppSimulate(newSpikeRate$caseRates)
			newSpikeCBSN5 = ScanCBS(controlSim, newSpikeSim, statistic="normal", takeN=5, minStat=8, maxNCut=50, timing=T)
			newCBSRes = list(spikeCBSN5=newSpikeCBSN5, spikeMat=newSpikeRate$spikeMat)
			simCBSRes[[counter]] = newCBSRes
			counter = counter + 1
		}
	}
	return(simCBSRes)
}

