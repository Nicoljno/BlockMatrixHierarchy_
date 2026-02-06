using BlockMatrixHierarchy
using MosekTools
using JuMP
import SCS            # or COSMO, or another solver
using MathOptInterface
using LinearAlgebra
using SparseArrays
using SplitApplyCombine
using Ket
using Plots
using FileIO, JLD2
using Hypatia
import ComplexOptInterface

function computational_basis_ρ(d::Integer)
    ρ = Matrix{ComplexF64}(I, d, d)               # identity
    return [ ρ[:,k] * ρ[k,:]'  for k in 1:d ]     # |k⟩⟨k|
end

function qft(d::Integer)
    ω = 2π * im / d
    F = [exp(ω*(j-1)*(k-1)) for j in 1:d, k in 1:d]
    return F / √d
end


d = 3
loc_d = d
num_parties::Int8  = 4*d
num_states::Int8  = 1
num_states_parties::Int8 = 4*d
inputs::Int8  = 0
outputs::Int8  = 0
level::Int8 = 1
state_kind::Int8 = 1
element_kind::Int8 = 1

arr=[]

eq = ["ρ11-ρ11"]
ineq = ["ρ11-ρ11"]

obj="ρ11"

F = qft(d)


#targets = computational_basis_ρ(d)
#targets[d]=ketbra(ones(d))/d

#"""
targets=[]
X = shift(d, 1)
vecs = eigvecs(X)
for i=1:d
    push!(targets, ketbra(vecs[:,i]))
end
Z = clock(d, 1)
vecs = eigvecs(Z)
count = 1
for i=d+1:2*d
    push!(targets, ketbra(vecs[:,count]))
    global count = count+1
end
for j=3:4
    global vecs = eigvecs(X*clock(d,j-1))
    global count = 1
    for i=j*d+1:(j+1)*d
        push!(targets, ketbra(vecs[:,count]))
        global count = count+1
    end
end
#"""


arr=[]

for r = d:d
    global model, var_dict, variables, G, gamma = BlockMatSDP(d, loc_d, num_parties, num_states, num_states_parties, inputs, outputs, level, eq, ineq, obj, state_kind, element_kind)
    global SymG = symbolic_matrix(gamma)

    global v = @variable(model, 0 <= v <= 1)
    for j=1:num_states_parties
        @constraint(model, var_dict[gamma[findfirst(x->x=="ρ1$(j)", SymG)]] == v*targets[j]+(1-v)*I(d)/d)
        #@constraint(model, var_dict[gamma[findfirst(x->x=="ρ1$(j)ρ1$(num_states_parties)", SymG)]] == var_dict[gamma[findfirst(x->x=="ρ1$(j)", SymG)]])
    end
    for j=1:num_states_parties
        for k=j+d:num_states_parties
            for l=1:num_states_parties
                for m=l+d:num_states_parties
                    @constraint(model, real(tr(var_dict[gamma[findfirst(x->x=="ρ1$(j)ρ1$(k)", SymG)]])) == real(tr(var_dict[gamma[findfirst(x->x=="ρ1$(l)ρ1$(m)", SymG)]])))
                end
            end
        end
    end
    #@constraint(model, tr(var_dict[gamma[findfirst(x->x=="ρ1$(num_states_parties)", SymG)]]) == r)
    #@constraint(model, var_dict[gamma[findfirst(x->x=="ρ$(d+1)", SymG)]] == diagm(diag(var_dict[gamma[findfirst(x->x=="ρ$(d+1)", SymG)]])))
    set_optimizer(model, Hypatia.Optimizer)
    @objective(model, Max, v)
    optimize!(model)

    push!(arr, [r, value(v)])
end
