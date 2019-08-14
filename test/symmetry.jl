using Test
using LinearAlgebra

push!(LOAD_PATH, "src/")
using LaddersAndResonances

const NSigma = 150
const s = 16*mPi^2 + 1e-5
const SigmaPointsWeights = SigmaPointsWeightsHelp(4*mPi^2,(sqrt(s) - mPi)^2,NSigma)

#
const B = BKernelSymmHelp(SigmaPointsWeights,s)
const BKernelPlotTable = BKernelPlotHelp(SigmaPointsWeights,s)
const M = MHelp(BKernelTable)
const DetM = det(M)

# @testset begin
@test (DetM ≈ 0) == false

InvM = inv(M)
LadderTable = LadderHelp(InvM*BKernelTable, SigmaPointsWeights, s)
LadderSymmetryCheckTable = LadderTable - transpose(LadderTable)
@test max(abs.(LadderSymmetryCheckTable)...) < 1e-8
# end
