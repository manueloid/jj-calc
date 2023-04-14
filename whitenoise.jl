include("./src/ConstantQuantities.jl")

# Testing the time_reverse function
tf = .2
cp = ControlParameterFull(tf,10)
h = 2.0/(cp.NParticles)
qts = ConstantQuantities(cp)
corrs = corrections(cp)
Λf(t) = Λ_esta(t, cp, corrs) # forward Λ
Λb(t) = time_reverse(t,Λf, cp) # backward Λ
ts = range(-.1tf, 1.1tf, length=100)
using Plots
plot(ts, Λf.(ts), label="Λ(t)")
plot!(ts, Λb.(ts), label="Λ_rev(t)")

# Testing the evolution function
ψ0f = qts.ψ0 # forward initial state
ψff = qts.ψf # forward final state
ψ0b = qts.ψf # backward initial state
ψfb = qts.ψ0 # backward final state
ts_sim = range(0.0, tf, length=100) # time steps for simulation
Hf = Hamiltonian(qts, Λf, h) # forward Hamiltonian
Hb = Hamiltonian(qts, Λb, h) # backward Hamiltonian

# Forward evolution
Ψ0 = evolution(ψ0f, Hf, ts_sim)[2] # wavefunction of the forward evolution at different time steps
ΨT = evolution(ψ0b, Hb, ts_sim)[2] # wavefunction of the backward evolution at different time steps

pr(ψa, ψb) = dagger.(ψa) .* ψb # projection of ψa onto ψb
pr(ψa) = dagger.(ψa) .* ψa # projection of ψa onto itself
pr(Ψ0) # checking normalization of Ψ0
pr(ΨT) # checking normalization of ΨT

correctness_f = [ dagger(Ψ) * ψ0f for Ψ in Ψ0 ] # checking that Ψ0 starts from ψ0f and then evolves
correctness_b = [ dagger(Ψ) * ψ0b for Ψ in ΨT ] # checking that ΨT starts from ψ0b and then evolves
plot(ts_sim, abs.(correctness_f), label="|<Ψ(t)|ψ0>| f")
plot!(ts_sim, abs.(correctness_b), label="|<Ψ(t)|ψ0>| b")

# Checking the correctness of the bf_evolution function
Ψ0_bf, ΨT_bf, H_bf = bf_evolution(cp, Λf; t = ts_sim) # state evolution using the bf_evolution function
Ψ0_bf - Ψ0 # checking that the two methods give the same result
ΨT_bf - ΨT # checking that the two methods give the same result

diff = H_bf - Hf.(ts_sim, 0) # checking that the Hamiltonian is correct
result = 2.0qts.Jx |> dense
[diff[i] - result for i in 1:length(diff)]


# plotting the sensitivities
final_times = range(0.1, 2.8, length=100)
esta_s, sta_s = sensitivities(10, final_times)
plot(final_times, abs.(esta_s), label="esta")
plot!(final_times, abs.(sta_s), label="sta")
