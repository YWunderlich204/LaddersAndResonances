using Test
using LinearAlgebra
using Plots

push!(LOAD_PATH, "src/")
using LaddersAndResonances

const NSigma = 150
const s = 16*mPi^2 + 1e-5
const σW = SigmaPointsWeightsHelp(4*mPi^2, (√s-mPi)^2, NSigma)
#
const Phi = Kibble.(σW[1],s,σW[1]')
balanceplot(Phi; xy=σW[1], title="Kibble function: Dalitz plot(>0) in legendre coordinates")
savefig(joinpath("plot","Kibble.pdf"))
#
const B = BKernelSymmHelp(σW,s)
plot(layout=grid(1,2), size=(1000,350),
    balanceplot(real.(B); xy=σW[1], title="Re B"),
    balanceplot(imag.(B); xy=σW[1], title="Im B", cscale=1))
savefig(joinpath("plot","B.pdf"))
