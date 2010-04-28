readSeq <-
function(filename, formatName) {
	if(formatName=="ELANDPaired") {
		return(readSeqELANDPaired(filename))
	}
	else if(formatName=="Chiang") {
		return(readSeqChiang(filename))
	}
	else if(formatName=="MLove") {
		return(readSeqMLove(filename))
	}
	else if(formatName=="JShen") {
		return(readSeqJShen(filename))
	}
	else {
		print("Input Format Not Supported")
	}
}

