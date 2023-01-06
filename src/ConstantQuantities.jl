using QuantumOptics
using JosephsonSTA

struct ConstantQuantities
    Jz::Operator
    Jx::Operator
    ψ0::Ket
    ψf::Ket
end

include("fidelities.jl")
