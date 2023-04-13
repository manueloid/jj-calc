using DataFrames, CSV
data_dir(nlambda::Int64) = "./data/nlambda$(nlambda)/" # Saving directory for different lambda
filename(np::Int64, feature="fidelity", nlambda::Int64=5) = data_dir(nlambda) * "$(feature)_np$np.dat" # Saving filename for given number of lambda and given number of particles
data(np::Int64, feature="fidelity", nlambda::Int64=5) = CSV.read(filename(np, feature, nlambda), DataFrame; header=true) |> s -> sort(s, :tf)
using Colors, PGFPlotsX
# Colors for the plots
colors = (
red =     colorant"#FF0000", # eSTA Full Hamiltonian with Hessian 
yellow =     colorant"#FFa000", # eSTA Full Hamiltonian with original version
blue =     colorant"#0000FF", # eSTA Intermediate Hamiltonian with Hessian
green =     colorant"#00FF00",# eSTA Intermediate Hamiltonian original version
black =     colorant"#000000" # STA
)

# Array where all the line styles are stored
styles = (
    solid="solid",                         # eSTA Intermediate Hamiltonian with Hessian
    ldash="dash pattern={on 2pt off 2pt}",# eSTA Full Hamiltonian with Hessian 
    dash=" dashed",                        # eSTA Full Hamiltonian with original version
    dot="dash pattern={on 1pt off 1pt}",  # eSTA Intermediate Hamiltonian original version
    dot_dash="dash pattern={on 4pt off 1pt on 1pt off 1pt}" #STA
)

names = (
    "Full Hessian",
    "Full Original",
    "Intermediate Hessian",
    "Intermediate Original",
    "STA"
)
