using LinearAlgebra
using GaussQuadrature


const mPi = 140.

const g = 5000

const SigmaRes = 6*mPi^2

const epsilon = 5*10^(-5)

const NMandelstamS = 10

const NSigma = 1500

function KaellenLambda(x,y,z)
	x^2 + y^2 + z^2 - 2*x*y - 2*x*z - 2*y*z
end


function rho(sigma)
	sqrt(Complex(KaellenLambda(sigma,mPi^2,mPi^2)))/(8*pi*sigma)
end


function Lambda12(s, sigma)
	sqrt(Complex(KaellenLambda(s,sigma,mPi^2)))
end


function Rho(s, sigma)
	sqrt(Complex(KaellenLambda(s,sigma,mPi^2)))/(8*pi*s)
end


function t(sigma)
     (g^2)/(SigmaRes- sigma - im*rho(sigma)*g^2)
end


function tau(s, sigma)
	(1/3)*Rho(s,sigma)*t(sigma)
end

function BKernel(sigma, s, sigmaPrime)
	((2*s)/(Lambda12(s, sigma)*Lambda12(s, sigmaPrime)))*log((2*s*sigmaPrime + 
     2*im*s*epsilon - (s + mPi^2 - sigma)*(s + sigmaPrime - mPi^2) - 
     Lambda12(s, sigma)*Lambda12(s, sigmaPrime))/(2*s*sigmaPrime + 
     2*im*s*epsilon - (s + mPi^2 - sigma)*(s + sigmaPrime - mPi^2) + 
     Lambda12(s, sigma)*Lambda12(s, sigmaPrime)))
 end


function SigmaPointsWeightsHelp(a,b,N)
	tmpPoints, tmpWeights = legendre(N)
	tmpPoints = ((b-a)/2)*tmpPoints .+ (a+b)/2
	tmpWeights = ((b-a)/2)*tmpWeights
	#tmpWeights = (1/(2*pi))*tmpWeights
	return [tmpPoints, tmpWeights]
end


const MandelstamSPoints = [i for i in range((9 + 0.00005)*mPi^2, stop=16*mPi^2, step=(16*mPi^2 - (9 + 0.00005)*mPi^2)/NMandelstamS)]

const SigmaPointsWeights = [SigmaPointsWeightsHelp(4*mPi^2,(sqrt(MandelstamSPoints[i]) - mPi)^2,NSigma) for i in range(1, stop=NMandelstamS, step=1)]

function BKernelSymmHelp(index,N)	
    tmp = [
          sqrt(Complex(SigmaPointsWeights[index][2][i]*tau(MandelstamSPoints[index], SigmaPointsWeights[index][1][i])))*
          BKernel(SigmaPointsWeights[index][1][i], MandelstamSPoints[index],SigmaPointsWeights[index][1][iPrime])*
          sqrt(Complex(SigmaPointsWeights[index][2][iPrime]*tau(MandelstamSPoints[index], SigmaPointsWeights[index][1][iPrime])))
    	  for i in 1:N, iPrime in 1:N]
	return tmp
end


function BKernelPlotHelp(index,N)	
    tmp = [
          BKernel(SigmaPointsWeights[index][1][i], MandelstamSPoints[index],SigmaPointsWeights[index][1][iPrime])
    	  for i in 1:N, iPrime in 1:N]
	return tmp
end


function KroneckerDelta(i,j)
    if i == j
    	tmp = 1
    else
    	tmp = 0
    end
    return tmp
end

const BKernelTable = [ BKernelSymmHelp(i,NSigma) for i in range(1, stop=NMandelstamS) ]

const BKernelPlotTable = [ BKernelPlotHelp(i,NSigma) for i in range(1, stop=NMandelstamS) ]

function MHelp(index,N)
	tmp = [
          KroneckerDelta(i,iPrime) - (1/(2*pi))*BKernelTable[index][i,iPrime]
		  for i in 1:N, iPrime in 1:N]
    return tmp
end

const M = [ MHelp(i,NSigma) for i in range(1, stop=NMandelstamS) ]

const DetM = [ det(M[i]) for i in range(1, stop=NMandelstamS) ]

display("Beginning of Inversion ...")

const InvM = [ inv(M[i]) for i in range(1, stop=NMandelstamS) ]

display("... End of Inversion")

LadderBarTable = [ InvM[i]*BKernelTable[i] for i in range(1, stop=NMandelstamS) ]

function LadderHelp(index, InputLadder, N)
     tmp = [
           (1/sqrt(Complex(SigmaPointsWeights[index][2][i]*tau(MandelstamSPoints[index], SigmaPointsWeights[index][1][i]))))*
           InputLadder[i,iPrime]*
           (1/sqrt(Complex(SigmaPointsWeights[index][2][iPrime]*tau(MandelstamSPoints[index], SigmaPointsWeights[index][1][iPrime]))))
    	   for i in 1:N, iPrime in 1:N]
	return tmp
end

const LadderTable = [ LadderHelp(i, LadderBarTable[i], NSigma) for i in range(1, stop=NMandelstamS) ]

const LadderSymmetryCheckTable = [ LadderTable[i] - transpose(LadderTable[i]) for i in range(1, stop=NMandelstamS) ]

const RescatteringCorrectionsTable = [ LadderTable[i] - BKernelPlotTable[i] for i in range(1, stop=NMandelstamS) ]

display(RescatteringCorrectionsTable[5])

#display(BKernelTable[5])

#display(BKernelPlotTable[5])