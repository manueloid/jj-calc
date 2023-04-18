include("./src/ConstantQuantities.jl")
epsilons = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11,]
tf=range(0.1,1.0,length=100)
using Plots
plotlyjs()
cparam = ControlParameterFull(0.2,10)
pl=plot()
for eps in epsilons
	robustnesses_esta, robustnesses_sta = robustnesses_tn(cparam, tf, eps)
	plot!(tf, robustnesses_esta, label="δ=$eps")
end
pl
