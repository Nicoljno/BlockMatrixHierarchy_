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

d = 9
num_parties::Int8 = 2
num_states::Int8 = 1
num_states_parties::Int8 = 2
inputs::Int8 = 2
outputs::Int8 = 2
level::Int8 = -1
loc_d = 3
state_kind::Int8 = 1
element_kind::Int8 = 1

function solve_SDP_step(n_loop, d, loc_d, num_parties, num_states, num_states_parties, inputs, outputs, level, eq, ineq, obj, state_kind, element_kind)
    ε = 0.01*(n_loop-1)
    ε = 0.063
    
    
    model, var_dict, variables, G, gamma = BlockMatSDP(
        d, loc_d, num_parties, num_states, num_states_parties, inputs, outputs, level, eq, ineq, obj, state_kind, element_kind
    )
    
    SymG = symbolic_matrix(gamma)
    println(length(SymG))
    
    # build constraints, etc.
    build_conjugate_constraints(model, variables, var_dict)
    tmp1 = [ketbra(eigvecs(clock(loc_d))[i,:]) for i=1:loc_d]
    tmp = [ketbra(eigvecs(shift(loc_d))[i,:]) for i=1:loc_d]
    A_targ_ = [tmp, tmp1]
    B_targ_ = [tmp, tmp1]
    A_targ = [Vector{Matrix{ComplexF64}}(undef, outputs) for _ in 1:inputs]
    B_targ = [Vector{Matrix{ComplexF64}}(undef, outputs) for _ in 1:inputs]
    
    tmp = @variable(model, [1:d, 1:d] in HermitianPSDCone())
    @constraint(model, partial_transpose((var_dict[gamma[findfirst(x->x=="ρ11ρ12", SymG)]]), 1) == tmp)
    tmp = @variable(model, [1:d, 1:d] in HermitianPSDCone())
    @constraint(model, partial_transpose((var_dict[gamma[findfirst(x->x=="ρ11ρ12", SymG)]]), 2) == tmp)
    for i=1:inputs
        for j=1:outputs
            A_targ[i][j]=LinearAlgebra.kron(A_targ_[i][j], I(loc_d))
            B_targ[i][j]=LinearAlgebra.kron(I(loc_d), B_targ_[i][j])
        end
    end
    
    for j=1:inputs
        @constraint(model, real(tr(var_dict[gamma[findfirst(x->x=="ρ1$(j)", SymG)]])) == loc_d)
        for i=1:outputs
            @constraint(model, Hermitian(var_dict[gamma[findfirst(x->x=="ρ11", SymG)]] - var_dict[gamma[findfirst(x->x=="ρ11A$(i)$(j)ρ11", SymG)]]) in HermitianPSDCone() )
            @constraint(model, Hermitian(var_dict[gamma[findfirst(x->x=="ρ12", SymG)]] - var_dict[gamma[findfirst(x->x=="ρ12B$(i)$(j)ρ12", SymG)]]) in HermitianPSDCone() )
            @constraint(model, real(tr(var_dict[gamma[findfirst(x->x=="A$(i)$(j)ρ12", SymG)]]*A_targ[j][i])) >= (1-ε))
            @constraint(model, real(tr(var_dict[gamma[findfirst(x->x=="ρ11B$(i)$(j)", SymG)]]*B_targ[j][i])) >= (1-ε))
            @constraint(model, real(tr(var_dict[gamma[findfirst(x->x=="A$(i)$(j)", SymG)]]*A_targ[j][i])) >= (1-ε)/loc_d)
            @constraint(model, real(tr(var_dict[gamma[findfirst(x->x=="B$(i)$(j)", SymG)]]*B_targ[j][i])) >= (1-ε)/loc_d)
        end
    end

    for nA=1:inputs
        for nB=1:inputs
            @constraint(model, real(tr(sum(var_dict[gamma[findfirst(x->x=="ρ11A$(oA)$(nA)B$(oB)$(nB)ρ12", SymG)]] for oA=1:outputs for oB=1:outputs))) <= 1 )
            for o=1:outputs
                for o1=1:outputs
                    @constraint(model, real(tr(var_dict[gamma[findfirst(x->x=="ρ11A$(o1)$(nA)B$(o)$(nB)ρ12", SymG)]])) >=0 )
                    @constraint(model, real(tr(var_dict[gamma[findfirst(x->x=="A$(o1)$(nA)ρ11B$(o)$(nB)", SymG)]]*B_targ[nB][o])) >= real(tr(var_dict[gamma[findfirst(x->x=="ρ11A$(o1)$(nA)", SymG)]]))*(1-ε)/loc_d )
                    @constraint(model, real(tr(var_dict[gamma[findfirst(x->x=="A$(o1)$(nA)B$(o)$(nB)ρ12", SymG)]]*A_targ[nA][o1])) >= real(tr(var_dict[gamma[findfirst(x->x=="ρ12B$(o)$(nB)", SymG)]]))*(1-ε)/loc_d )
                    @constraint(model, real(tr(var_dict[gamma[findfirst(x->x=="A$(o1)$(nA)B$(o)$(nB)", SymG)]]*A_targ[nA][o1]*B_targ[nB][o])) >= (1-ε)^2 )
                end
                @constraint(model, real(tr(sum( var_dict[gamma[findfirst(x->x=="ρ11A$(i)$(nA)B$(o)$(nB)ρ12", SymG)]] for i=1:outputs ))) <= real(tr(var_dict[gamma[findfirst(x->x=="ρ11ρ12B$(o)$(nB)", SymG)]])) )
                @constraint(model, real(tr(sum( var_dict[gamma[findfirst(x->x=="ρ11A$(o)$(nA)B$(i)$(nB)ρ12", SymG)]] for i=1:outputs ))) <= real(tr(var_dict[gamma[findfirst(x->x=="ρ11A$(o)$(nA)ρ12", SymG)]])) )
            end
        end
    end

    #"""
    for v in variables
        tmp_A = [item for item in v if item.party==1]
        tmp_B = [item for item in v if item.party==2]
        if length(tmp_A)==1 && length(v)>1
            v_ = SymbolicMonomial[item for item in v if item.party!=1]
            v_ = partial_reduce_monomial(v_, 2)
            if v_ in variables
                @constraint(model, partial_trace(var_dict[v], 1) == partial_trace(var_dict[v_], 1)/loc_d )
            elseif reverse(v_) in variables
                @constraint(model, partial_trace(var_dict[v], 1) == partial_trace(var_dict[reverse(v_)], 1)'/loc_d )
            end
        end
        if length(tmp_B)==1 && length(v)>1
            v_ = SymbolicMonomial[item for item in v if item.party!=2]
            v_ = partial_reduce_monomial(v_, 1)
            if v_ in variables
                @constraint(model, partial_trace(var_dict[v], 2) == partial_trace(var_dict[v_], 2)/loc_d )
            elseif reverse(v_) in variables
                @constraint(model, partial_trace(var_dict[v], 2) == partial_trace(var_dict[reverse(v_)], 2)'/loc_d )
            end
        end
    end
    model = BlockMatrixHierarchy.trace_constraints(model, variables, var_dict)
    
    set_optimizer(model, Hypatia.Optimizer)
    #set_optimizer_attribute(model, "MSK_IPAR_NUM_THREADS", 4)
    println("SDP start")
    optimize!(model)
    
    val = objective_value(model)
    
    return (ε, val, model, var_dict, variables, G, gamma)
end

#eq = ["ρ11ρ12-I/$d"]
#ineq = ["ρ11A11ρ12B11+I/$d"]

#obj="ρ11A11ρ12B11+ρ11A12ρ12B12+ρ11A21ρ12B21+ρ11A22ρ12B22+ρ11A31ρ12B31+ρ11A32ρ12B32"

eq = ["ρ11-ρ11"]
ineq = ["ρ11"]
#obj="2*ρ11A11B11ρ12+2*ρ11A12B12ρ12+2*ρ11ρ12-ρ11A12ρ12-ρ11A11ρ12-ρ11ρ12B12-ρ11ρ12B11"
obj="2*ρ11A11B11ρ12+2*ρ11A12B12ρ12+2*ρ11A21B21ρ12+2*ρ11A22B22ρ12+2*ρ11ρ12-ρ11A11ρ12-ρ11A12ρ12-ρ11A21ρ12-ρ11A22ρ12-ρ11ρ12B11-ρ11ρ12B12-ρ11ρ12B21-ρ11ρ12B22+ρ11A11B21ρ12+ρ11A21B11ρ12+ρ11A22B12ρ12+ρ11A12B22ρ12"
#obj = "ρ11A11B11ρ12+ρ11A12B12ρ12+ρ11A21B21ρ12+ρ11A22B22ρ12+ρ11A31B31ρ12+ρ11A32B32ρ12"
arr = []
for n_loop in 1:1
    global (ε, val, model, var_dict, variables, G, gamma) = solve_SDP_step(n_loop, d, loc_d, num_parties, num_states, num_states_parties, inputs, outputs, level, eq, ineq, obj, state_kind, element_kind)
    push!(arr, [ε, val])
    GC.gc()
end

