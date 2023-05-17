include("./src/ConstantQuantities.jl")
using ProgressMeter
# function to calculate the fidelity for different final times
# first I will calculate the fidelity for the esta protocol with the discrete Hamiltonian and 5 λ.

final_times = range(0.01pi, 0.2pi, length=100)
function fidelities(np::Int64, final_times::AbstractVector)
	cparam = ControlParameterFull(np, final_times[1])
	qts = ConstantQuantities(cparam)
	fid_esta = zeros(length(final_times))
	fid_sta = zeros(length(final_times))
	fid_esta1 = zeros(length(final_times))
	fid_esta_cont = zeros(length(final_times))
	fid_ad = zeros(length(final_times))
	senstn_esta = zeros(length(final_times))
	senstn_sta = zeros(length(final_times))
	sensmn_esta = zeros(length(final_times))
	sensmn_sta = zeros(length(final_times))
	Threads.@threads for index in 1:length(final_times)
		tf = final_times[index]
		cparam = cp_time(cparam, tf)
		cparam_cont = ControlParameterInt(np, tf)
		corrs = corrections(cparam)
		corrs1 = corrections(cparam; nlambda = 1 )
		corrs_cont = corrections(cparam_cont)
		esta(t) = Λ_esta(t, cparam, corrs)
		esta1(t) = Λ_esta(t, cparam, corrs1)
		esta_cont(t) = Λ_esta(t, cparam_cont, corrs_cont)
		sta(t) = Λ_sta(t, cparam)
		ad(t) = Λ_ad(t, cparam)
		fid_esta[index] = fidelity(cparam, qts, esta)
		fid_esta1[index] = fidelity(cparam, qts, esta1)
		fid_sta[index] = fidelity(cparam, qts, sta)
		fid_esta_cont[index] = fidelity(cparam_cont, qts, esta_cont)
		fid_ad[index] = fidelity(cparam, qts, ad)
		senstn_esta[index] = robustness_tn(cparam, qts, esta, 1e-7)
		senstn_sta[index] = robustness_tn(cparam, qts, sta, 1e-7)
		sensmn_esta[index] = robustness_mn(cparam, qts, esta, 1e-7)
		sensmn_sta[index] = robustness_mn(cparam, qts, sta, 1e-7)
	end
	return fid_esta, fid_esta1,fid_esta_cont, fid_sta, fid_ad, senstn_esta, senstn_sta, sensmn_esta, sensmn_sta
end
fidelities(10, final_times)

fidelities(30, final_times)
