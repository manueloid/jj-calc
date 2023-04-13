include("./src/ConstantQuantities.jl")
cp = ControlParameterFull(0.1, 10) # I will use 10 particles to start with
final_times = range(0.2, 0.8, length=100) |> collect
esta_sq, sta_sq = squeezings(cp, final_times)

using PGFPlotsX

