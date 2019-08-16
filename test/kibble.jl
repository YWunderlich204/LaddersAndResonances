using Test

push!(LOAD_PATH, "src/")
using LaddersAndResonances

let s = 16mPi^2,
    σ = (4mPi^2 + (√s-mPi)^2) / 2
    @test Kibble(σ,s,σ) > 0 && (Kibble(4mPi^2,s,4mPi^2) < 0)
end
