using BlockMatrixHierarchy
using MosekTools
using JuMP
using MathOptInterface
using LinearAlgebra
using SparseArrays
using SplitApplyCombine
using Ket
using Hypatia
using FileIO, JLD2
using Plots
import ComplexOptInterface

d = 4
num_parties::Int8 = 2
num_states::Int8 = 1
num_states_parties::Int8 = 2
inputs::Int8 = 2
outputs::Int8 = 1
level::Int8 = -1
loc_d = 2
state_kind::Int8 = 1 
element_kind::Int8 = 0

function solve_SDP_step(n_loop, d, loc_d, num_parties, num_states, num_states_parties, inputs, outputs, level, eq, ineq, obj, state_kind, element_kind)
    ε = 0.01*(n_loop-1)
    εA = [0.01, 0.02]
    εB = [0.013, 0.023]
    
    """
    model, var_dict, variables, partial_variables, G, gamma, partial_G, partial_gamma = BlockMatSDP_loc(
        d, loc_d, num_parties, num_states, inputs, outputs, level, eq, ineq, obj
    )
    #"""
    
    #"""
    model, var_dict, variables, G, gamma = BlockMatSDP(
        d, loc_d, num_parties, num_states, num_states_parties, inputs, outputs, level, eq, ineq, obj, state_kind, element_kind
    )
    #"""
    SymG = BlockMatrixHierarchy.symbolic_matrix(gamma)

    # build constraints, etc.
    BlockMatrixHierarchy.build_conjugate_constraints(model, variables, var_dict)

    A_targ = [pauli(1), pauli(3)]
    B_targ = [pauli(1), pauli(3)]
    for i=1:2
        A_targ[i]=Hermitian(LinearAlgebra.kron(A_targ[i], I(2)))
        B_targ[i]=Hermitian(LinearAlgebra.kron(I(2), B_targ[i]))
    end
    #z=-pi/4+atan(1/(8*ε-8*ε^2-1))/4
    for i=1:loc_d
        #@constraint(model, var_dict[gamma[findfirst(x->x=="ρ11", SymG)]]==LinearAlgebra.kron(ketbra([cos(z), sin(z)]), I(2)))
        #@constraint(model, var_dict[gamma[findfirst(x->x=="ρ12", SymG)]]==LinearAlgebra.kron(I(2), ketbra([cos(z), sin(z)])))
        @constraint(model, real(tr(var_dict[gamma[findfirst(x->x=="A1$(i)", SymG)]]*A_targ[i])) >= loc_d^2*(1-2*εA[i]))
        @constraint(model, real(tr(var_dict[gamma[findfirst(x->x=="B1$(i)", SymG)]]*B_targ[i])) >= loc_d^2*(1-2*εB[i]))
        @constraint(model, real(tr(var_dict[gamma[findfirst(x->x=="A1$(i)", SymG)]])) == 0)
        @constraint(model, real(tr(var_dict[gamma[findfirst(x->x=="B1$(i)", SymG)]])) == 0)
        @constraint(model, real(tr(var_dict[gamma[findfirst(x->x=="ρ1$(i)", SymG)]])) == loc_d)
        for j=1:loc_d
            @constraint(model, real(tr(var_dict[gamma[findfirst(x->x=="A1$(i)B1$j", SymG)]]*A_targ[i]*B_targ[j])) >= loc_d^2*(1-2*εA[i])*(1-2*εB[j]))
        end
    end
    
    for v in variables
        tmp_A = [item for item in v if item.party==1]
        tmp_B = [item for item in v if item.party==2]
        if length(tmp_A)==1 && tmp_A isa Vector{BlockMatrixHierarchy.SymbolicElement}
            @constraint(model, partial_trace(var_dict[v], 1) == zeros(loc_d, loc_d))
        elseif length(tmp_A)==1 && tmp_A isa Vector{BlockMatrixHierarchy.SymbolicState}
            v_ = BlockMatrixHierarchy.SymbolicMonomial[item for item in v if item.party!=1]
            v_ = BlockMatrixHierarchy.partial_reduce_monomial(v_, 2)
            if v_ in variables
                @constraint(model, partial_trace(var_dict[v], 1) == partial_trace(var_dict[v_], 1)/loc_d )
            else
                @constraint(model, partial_trace(var_dict[v], 1) == partial_trace(var_dict[reverse(v_)], 1)'/loc_d )
            end
        elseif length(tmp_B)==1 && tmp_B isa Vector{BlockMatrixHierarchy.SymbolicElement}
            @constraint(model, partial_trace(var_dict[v], 2) == zeros(loc_d, loc_d))
        elseif length(tmp_B)==1 && tmp_B isa Vector{BlockMatrixHierarchy.SymbolicState}
            v_ = BlockMatrixHierarchy.SymbolicMonomial[item for item in v if item.party!=2]
            v_ = BlockMatrixHierarchy.partial_reduce_monomial(v_, 1)
            if v_ in variables
                @constraint(model, partial_trace(var_dict[v], 2) == partial_trace(var_dict[v_], 2)/loc_d )
            else
                @constraint(model, partial_trace(var_dict[v], 2) == partial_trace(var_dict[reverse(v_)], 2)'/loc_d )
            end  
        end
    end
    @constraint(model, LinearAlgebra.kron(I(loc_d), partial_trace(var_dict[gamma[findfirst(x->x=="ρ11ρ12", SymG)]], 1)) == var_dict[gamma[findfirst(x->x=="ρ12", SymG)]] )
    @constraint(model, LinearAlgebra.kron(partial_trace(var_dict[gamma[findfirst(x->x=="ρ11ρ12", SymG)]], 2), I(loc_d)) == var_dict[gamma[findfirst(x->x=="ρ11", SymG)]] )
    @constraint(model, var_dict[variables[1]] == I)
    set_optimizer(model, Mosek.Optimizer)
    optimize!(model)
    
    val = objective_value(model)
    
    return (εA, εB, val, model, var_dict, variables, G, gamma)
end

eq = ["ρ11ρ12-I/$d"]
ineq = ["ρ11A11B11ρ12+I/$d", "ρ11A12B11ρ12+I/$d", "ρ11A11B12ρ12+I/$d", "ρ11A12B12ρ12+I/$d"]
ineq = [ineq; ["I/$d-ρ11A11B11ρ12", "I/$d-ρ11A12B11ρ12", "I/$d-ρ11A11B12ρ12", "I/$d-ρ11A12B12ρ12"]]
#obj="ρ11A11B11ρ12+ρ11A12B11ρ12+ρ11A11B12ρ12-ρ11A12B12ρ12"      #CHSH
obj="ρ11A11B11ρ12+ρ11A12B12ρ12"                                 #ENTANGLEMENT WIT
#ineq = ["ρ1A11B11+I/$d", "ρ1A12B11+I/$d", "ρ1A11B12+I/$d", "ρ1A12B12+I/$d"]
#ineq = [ineq; ["I/$d-ρ1A11B11", "I/$d-ρ1A12B11", "I/$d-ρ1A11B12", "I/$d-ρ1A12B12"]]
#obj="ρ1A11B11+ρ1A12B11+ρ1A11B12-ρ1A12B12"

#"""
arr = []
for n_loop in 1:1
    global (εA, εB, val, model, var_dict, variables, G, gamma) = solve_SDP_step(n_loop, d, loc_d, num_parties, num_states, num_states_parties, inputs, outputs, level, eq, ineq, obj, state_kind, element_kind)
    push!(arr, [εA, εB, val])
    GC.gc()
end
#"""

"""
model, var_dict, variables, G, gamma = BlockMatSDP(
    d, loc_d, num_parties, num_states, num_states_parties, inputs, outputs, level, eq, ineq, obj, state_kind, element_kind
)

SymG = BlockMatrixHierarchy.symbolic_matrix(gamma)
#"""


#4*(1-2*ε)*sqrt(ε*(1-ε))+sqrt(2-16*ε*(1-ε)*(1-2*ε)^2)
1+4*(1-2*ε)*sqrt(ε*(1-ε))