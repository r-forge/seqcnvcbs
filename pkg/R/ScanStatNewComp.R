ScanStatNewComp <-
function(combZCumSum, combXCumSum, combZPoint, combXPoint, p, nTotal, grid.cur, max.win, statistic) {
	if(statistic=="rabinowitz") {
		SijFactorR = p*(1-p)*(1-1/(nTotal-1))
		newRes = .Call("ScanStatNewCompRabinC", combZCumSum, combXCumSum, combZPoint, combXPoint, SijFactorR, p, nTotal, grid.cur, max.win)
	}
	else if(statistic=="normal") {
		SijFactorN = p*(1-p)
		newRes = .Call("ScanStatNewCompNormalC", combZCumSum, combXCumSum, combZPoint, combXPoint, SijFactorN, p, nTotal, grid.cur, max.win)
	}
	else if(statistic=="hybrid") {
		SijFactorH = p*(1-p)
		newRes = .Call("ScanStatNewCompHybridC", combZCumSum, combXCumSum, combZPoint, combXPoint, SijFactorH, p, nTotal, grid.cur, max.win)
	}
	else if(statistic=="binomial") {
		newRes = .Call("ScanStatNewCompBinomC", combZCumSum, combXCumSum, combZPoint, combXPoint, p, nTotal, grid.cur, max.win)
	}
	else {
		print("Name of Statistic is not Recognized, Rabinowitz used")
		SijFactorR = p*(1-p)*(1-1/(nTotal-1))
		newRes = .Call("ScanStatNewCompRabinC", combZCumSum, combXCumSum, combZPoint, combXPoint, SijFactorR, p, nTotal, grid.cur, max.win)
	}
	return(newRes)
}

