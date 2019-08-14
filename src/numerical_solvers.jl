using LinearAlgebra

function SigmaPointsWeightsHelp(a,b,N)
	tmpPoints, tmpWeights = legendre(N)
	tmpPoints = ((b-a)/2)*tmpPoints .+ (a+b)/2
	tmpWeights = ((b-a)/2)*tmpWeights
	#tmpWeights = (1/(2*pi))*tmpWeights
	return [tmpPoints, tmpWeights]
end

function BKernelSymmHelp(SigmaPointsWeights,s) # SigmaPointsWeights[index] -> SigmaPointsWeights
	N = length(SigmaPointsWeights[1])
    tmp = [
          sqrt(Complex(SigmaPointsWeights[2][i]*tau(s, SigmaPointsWeights[1][i])))*
          BKernel(SigmaPointsWeights[1][i], s,SigmaPointsWeights[1][iPrime])*
          sqrt(Complex(SigmaPointsWeights[2][iPrime]*tau(s, SigmaPointsWeights[1][iPrime])))
    	  	for i in 1:N, iPrime in 1:N]
	return tmp
end

function BKernelPlotHelp(SigmaPointsWeights,s)
	N = length(SigmaPointsWeights[1])
    tmp = [
          BKernel(SigmaPointsWeights[1][i], s, SigmaPointsWeights[1][iPrime])
    	  	for i in 1:N, iPrime in 1:N]
	return tmp
end

KroneckerDelta(i,j) = (i == j) ? 1.0 : 0.0

function MHelp(BKernelTable)
	N = size(BKernelTable,1)
	tmp = Diagonal(fill(1.0,N)) - (BKernelTable ./ (2π))
	# tmp = [
    #       KroneckerDelta(i,iPrime) - (1/(2*pi))*BKernelTable[i,iPrime]
	# 	  for i in 1:N, iPrime in 1:N]
    return tmp
end

function LadderHelp(InputLadder, SigmaPointsWeights, s)
	N = length(SigmaPointsWeights[1])
    tmp = [
           (1/sqrt(Complex(SigmaPointsWeights[2][i]*tau(s, SigmaPointsWeights[1][i]))))*
           InputLadder[i,iPrime]*
           (1/sqrt(Complex(SigmaPointsWeights[2][iPrime]*tau(s, SigmaPointsWeights[1][iPrime]))))
    	   		for i in 1:N, iPrime in 1:N]
	return tmp
end
