include("../src/ConstantQuantities.jl")
using Plots
cp = ControlParameterFull(0.9, 10)
ts(cp::ControlParameterFull, nsteps=1000) = range(0.0, stop=cp.final_time, length=nsteps)
qts = ConstantQuantities(cp)
corrs = corrections(cp)
esta(t) = Λ_esta(t, cp, corrs)
sta(t) = Λ_sta(t, cp)
ad(t) = Λ_ad(t, cp)
plot(ts(cp), esta.(ts(cp)), label="esta")
plot!(ts(cp), sta.(ts(cp)), label="sta")
plot!(ts(cp), ad.(ts(cp)), label="ad")

nsteps = 1000
t = ts(cp, nsteps)
α_esta = alpha(cp, qts, esta, nsteps)
α_sta = alpha(cp, qts, sta, nsteps)
α_ad = alpha(cp, qts, ad, nsteps)
plot(t, α_esta, label="esta")
plot!(t, α_sta, label="sta")
plot!(t, α_ad, label="ad")

ξN_esta = squeezing(cp, qts, esta, nsteps)
ξN_sta = squeezing(cp, qts, sta, nsteps)
ξN_ad = squeezing(cp, qts, ad, nsteps)
ξs_esta = ξN_esta .^ 2 ./ (α_esta .^ 2)
ξs_sta = ξN_sta .^ 2 ./ (α_sta .^ 2)
ξs_ad = ξN_ad .^ 2 ./ (α_ad .^ 2)
todecibel(x) = 10 * log10(x)
plot(t, todecibel.(ξs_esta), label="esta", legend=nothing)
plot!(t, todecibel.(ξs_sta), label="sta")
plot!(t, todecibel.(ξs_ad), label="ad")

fid_esta = fidelity(cp, qts, esta, nsteps)
fid_sta = fidelity(cp, qts, sta, nsteps)
fid_ad = fidelity(cp, qts, ad, nsteps)
plot(t, fid_esta)
plot!(t, fid_sta)
plot!(t, fid_ad)

# Now I need to plot for different final times
final_times = range(0.5, 4.0, length=100)
ξs_esta = zeros(length(final_times))
ξs_sta = zeros(length(final_times))
Threads.@threads for i in 1:length(final_times)
	cparam = cp_time(cp, final_times[i])
	corrs = corrections(cparam)
	esta(t) = Λ_esta(t, cparam, corrs)
	sta(t) = Λ_sta(t, cparam)
	α_esta = alpha(cparam, qts, esta, nsteps)[end]
	α_sta = alpha(cparam, qts, sta, nsteps)[end]
	ξs_esta[i] = squeezing(cparam, qts, esta, nsteps)[end] ^ 2 / α_esta ^ 2
	ξs_sta[i] = squeezing(cparam, qts, sta, nsteps)[end] ^ 2 / α_sta ^ 2
end

plot(final_times, todecibel.(ξs_esta), label="esta")
plot!(final_times, todecibel.(ξs_sta), label="sta")
