push!(LOAD_PATH, "src/")

using LaddersAndResonances
using LinearAlgebra

const NSigma = 150

const s = 16*mPi^2 + 1e-5
const SigmaPointsWeights = SigmaPointsWeightsHelp(4*mPi^2,(sqrt(s) - mPi)^2,NSigma)
const BKernelTable = BKernelSymmHelp(SigmaPointsWeights,s)

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
# const BKernelPlotTable = BKernelPlotHelp(SigmaPointsWeights,s)
# const RescatteringCorrectionsTable = LadderTable - BKernelPlotTable
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
	plot(layout=grid(1,2), size=(1000,350),
		balanceplot(real.(LadderTable); xy=σv, title="Re L"),
		balanceplot(imag.(LadderTable); xy=σv, title="Im L"),
		)
end
savefig(joinpath("plot","get_ladder.test.pdf"))
