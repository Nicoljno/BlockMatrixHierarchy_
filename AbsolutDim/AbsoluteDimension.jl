using BlockMatrixHierarchy
using MosekTools
using JuMP
using Hypatia
using MathOptInterface
using LinearAlgebra
using SparseArrays
using Ket

function computational_basis_ρ(d::Integer)
    ρ = Matrix{ComplexF64}(I, d, d)               # identity
    return [ ρ[:,k] * ρ[k,:]'  for k in 1:d ]     # |k⟩⟨k|
end

function qft(d::Integer)
    ω = 2π * im / d
    F = [exp(ω*(j-1)*(k-1)) for j in 1:d, k in 1:d]
    return F / √d
end

function amplitude_dampening(ensemble::Vector{Any}, γ::AbstractVector)
    d = length(γ) + 1

    # Basic validation (optional but recommended)
    for (i, g) in enumerate(γ)
        (0.0 ≤ g ≤ 1.0) || throw(ArgumentError("γ[$i] = $g is not in [0,1]."))
    end

    Ks = Vector{AbstractMatrix{Float64}}(undef, d)

    diag_entries = Vector{Float64}(undef, d)
    diag_entries[1] = 1.0
    for i in 2:d
        diag_entries[i] = sqrt(1 - float(γ[i-1]))
    end
    Ks[1] = Diagonal(diag_entries)
    for k in 1:(d-1)
        Kk = zeros(Float64, d, d)
        Kk[k, k+1] = sqrt(float(γ[k]))
        Ks[k+1] = Kk
    end

    for i=1:length(ensemble)
        ensemble[i]=sum(Hermitian(Ks[j]*ensemble[i]*Ks[j]') for j=1:d)
    end
    return ensemble

end

d = 8
loc_d = d
num_parties::Int8  = 1
num_states::Int8  = d+1
num_states_parties::Int8 = -1
inputs::Int8  = 0
outputs::Int8  = 0
level::Int8 = 1
state_kind::Int8 = 1
element_kind::Int8 = 1

arr=[]

eq = ["ρ1-ρ1"]
ineq = ["ρ1-ρ1"]
#obj="1/3*(ρ1A11+ρ2A21+ρ3A31+ρ4B11+ρ5B21+ρ6B31)"

# d=3
#obj="1/3*(ρ1A11+ρ2A21+ρ3A31+ρ4B11+ρ5B21+ρ6B31)"
# d=4
#obj="1/4*(ρ1A11+ρ2A21+ρ3A31+ρ4A41+ρ5B11+ρ6B21+ρ7B31+ρ8B41)"
# d=5
obj="ρ1"

#ineq = ["ρ1A11", "ρ1B11", "ρ2A11", "ρ2B11", "ρ3A11", "ρ3B11", "ρ4A11", "ρ4B11"]
#ineq = [ineq; ["I/$d-ρ1A11", "I/$d-ρ1B11", "I/$d-ρ2A11", "I/$d-ρ2B11", "I/$d-ρ3A11", "I/$d-ρ3B11", "I/$d-ρ4A11", "I/$d-ρ4B11"] ]
#obj="1/8*(ρ1A11+ρ3A11+ρ1B11+ρ2B11-ρ4A11-ρ4B11-ρ2A11-ρ3B11)"

F = qft(d)

targets = []
for i=1:d
    push!(targets, ketbra(F[i,:]))
end
#targets=amplitude_dampening(targets, ones(d-1)*0.9096)
#targets = computational_basis_ρ(d)
#targets[d]=ketbra(ones(d))/d

"""
targets=[]
X = shift(d, 1)
vecs = eigvecs(X)
push!(targets, ketbra(vecs[:,1]))
Z = clock(d, 1)
vecs = eigvecs(Z)
push!(targets, ketbra(vecs[:,1]))
for i=3:d
    global vecs = eigvecs(X*clock(d,i-1))
    push!(targets, ketbra(vecs[:,1]))
end
#"""


arr=[]

for r = 7:7
    global model, var_dict, variables, G, gamma = BlockMatSDP(d, loc_d, num_parties, num_states, num_states_parties, inputs, outputs, level, eq, ineq, obj, state_kind, element_kind)
    global SymG = symbolic_matrix(gamma)

    global v = @variable(model, 0 <= v <= 1)
    for j=1:num_states-1

        # case 1: depolarizing noise model 
        @constraint(model, var_dict[gamma[findfirst(x->x=="ρ$(j)", SymG)]] == v*targets[j]+(1-v)*I(d)/d)
        @constraint(model, var_dict[gamma[findfirst(x->x=="ρ$(j)ρ$(d+1)", SymG)]] == var_dict[gamma[findfirst(x->x=="ρ$(j)", SymG)]])

        # case 2: phase dampening 
        #@constraint(model, var_dict[gamma[findfirst(x->x=="ρ$(j)", SymG)]] == targets[j])
        #@constraint(model, var_dict[gamma[findfirst(x->x=="ρ$(j)ρ$(d+1)", SymG)]] == targets[j])
    end
    @constraint(model, tr(var_dict[gamma[findfirst(x->x=="ρ$(d+1)", SymG)]]) == r)
    #@constraint(model, var_dict[gamma[findfirst(x->x=="ρ$(d+1)", SymG)]] == diagm(diag(var_dict[gamma[findfirst(x->x=="ρ$(d+1)", SymG)]])))
    #set_optimizer(model, Hypatia.Optimizer)
    @objective(model, Max, v)
    optimize!(model)

    push!(arr, [r, value(v)])
end
