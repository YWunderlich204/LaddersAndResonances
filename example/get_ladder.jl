push!(LOAD_PATH, "src/")

using LaddersAndResonances
using LinearAlgebra

const NSigma = 150

const s = 16*mPi^2 + 1e-5
const SigmaPointsWeights = SigmaPointsWeightsHelp(4*mPi^2,(sqrt(s) - mPi)^2,NSigma)
const BKernelTable = BKernelSymmHelp(SigmaPointsWeights,s)
const BKernelPlotTable = BKernelPlotHelp(SigmaPointsWeights,s)


#                              _|                _|
#  _|_|_|  _|_|      _|_|_|  _|_|_|_|  _|  _|_|      _|    _|
#  _|    _|    _|  _|    _|    _|      _|_|      _|    _|_|
#  _|    _|    _|  _|    _|    _|      _|        _|  _|    _|
#  _|    _|    _|    _|_|_|      _|_|  _|        _|  _|    _|

#
const M = MHelp(BKernelTable)
const DetM = det(M)
@show DetM
const InvM = inv(M)

const LadderBarTable = InvM*BKernelTable
const LadderTable = LadderHelp(LadderBarTable, SigmaPointsWeights, s)
#
const RescatteringCorrectionsTable = LadderTable - BKernelPlotTable
# display(RescatteringCorrectionsTable[5])
const LadderSymmetryCheckTable = LadderTable - transpose(LadderTable)
@show extrema(abs.(LadderSymmetryCheckTable))


#            _|              _|
#  _|_|_|    _|    _|_|    _|_|_|_|    _|_|_|
#  _|    _|  _|  _|    _|    _|      _|_|
#  _|    _|  _|  _|    _|    _|          _|_|
#  _|_|_|    _|    _|_|        _|_|  _|_|_|
#  _|
#  _|

using Plots
let σv = SigmaPointsWeights[1]
	plot(layout=grid(1,2), size=(1000,350))
	heatmap!(sp=1, σv, σv, real.(LadderTable), title="Re(Ladder)")
	heatmap!(sp=2, σv, σv, imag.(LadderTable), title="Im(Ladder)")
end
savefig(joinpath("plot","get_ladder.test.pdf"))
