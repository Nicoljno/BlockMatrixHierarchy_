using BlockMatrixHierarchy
using MosekTools
using JuMP
using MathOptInterface
using LinearAlgebra
using SparseArrays
using SplitApplyCombine
using Ket
using Plots
using FileIO, JLD2
using Hypatia
import ComplexOptInterface

zeta(i, j, n) = (abs(i - j) == 1 || abs(i - j) == n - 1) ? -1 : 0
function upper_triangular(A::AbstractMatrix{String}; fill = "")
    mask = triu(trues(size(A)...))      # Bool mask: true on upper triangle
    return ifelse.(mask, A, fill)       # keep A where mask=true, else `fill`
end


d = 2
loc_d = 2
num_parties::Int8  = 1
num_states::Int8  = 1
num_states_parties::Int8 = -1
inputs::Int8  = 21
outputs::Int8  = 3
level::Int8 = 1
state_kind::Int8 = 1
element_kind::Int8 = -1

arr=[]

eq = ["ρ1-ρ1"]
ineq = ["ρ1-ρ1"]

obj="ρ1A11+ρ1A12+ρ1A13+ρ1A14+ρ1A15+ρ1A16+ρ1A17+ρ1A18+ρ1A19+ρ1A110+ρ1A111+ρ1A112+ρ1A113+ρ1A114+ρ1A115+ρ1A116+ρ1A117+ρ1A118+ρ1A119+ρ1A120+ρ1A121"#

val=[]

for eps_=0:0
    global model, var_dict, variables, G, gamma = BlockMatSDP_real(d, loc_d, num_parties, num_states, num_states_parties, inputs, outputs, level, eq, ineq, obj, state_kind, element_kind)
    global SymG = symbolic_matrix(gamma)

    @constraint(model, real(tr(var_dict[gamma[findfirst(x->x=="ρ1", SymG)]])) == 1)

    for i = 1:inputs
        @constraint(model, var_dict[gamma[findfirst(x->x=="A1$(i)A1$(i)", SymG)]] == tr(var_dict[gamma[findfirst(x->x=="ρ1A1$(i)", SymG)]])*I(d))
        @constraint(model, tr(var_dict[gamma[findfirst(x->x=="A1$i", SymG)]])==0)
        @constraint(model, tr(var_dict[gamma[findfirst(x->x=="A2$i", SymG)]])==0)
        @constraint(model, tr(var_dict[gamma[findfirst(x->x=="A3$(i)", SymG)]])==0)
    end

    for i=1:inputs-1
        @constraint(model, var_dict[gamma[findfirst(x->x=="A2$(i)A2$(i)", SymG)]]==var_dict[gamma[findfirst(x->x=="A1$(i+1)A1$(i+1)", SymG)]])
    end
    @constraint(model, var_dict[gamma[findfirst(x->x=="A2$(inputs)A2$(inputs)", SymG)]]==var_dict[gamma[findfirst(x->x=="A11A11", SymG)]])

    # O^2=I && O^2<O>=<O>*I
    for i=1:inputs
        @constraint(model, var_dict[gamma[findfirst(x->x=="A3$(i)A3$(i)", SymG)]] == I(d))
        @constraint(model, var_dict[gamma[findfirst(x->x=="A1$(i)A3$(i)", SymG)]] == tr(var_dict[gamma[findfirst(x->x=="ρ1A3$(i)", SymG)]]) * I(d) )    
    end

    for i=1:inputs-1
        var_dict[gamma[findfirst(x->x=="A1$(i)A2$(i)", SymG)]] = @variable(model, [1:d,1:d])
        local lim=Hermitian(eps_*0.01*I(d)*tr(var_dict[gamma[findfirst(x->x=="A1$(i)A2$(i)", SymG)]]/d))
        @constraint(model, tr(var_dict[gamma[findfirst(x->x=="ρ1A2$(i)", SymG)]])*d == tr(var_dict[gamma[findfirst(x->x=="A1$(i)A2$(i)", SymG)]]) )
        @constraint(model, I-Hermitian(var_dict[gamma[findfirst(x->x=="A1$(i)A2$(i)", SymG)]]) in HermitianPSDCone())
        @constraint(model, lim-Hermitian(var_dict[gamma[findfirst(x->x=="A1$(i)A1$(i+1)", SymG)]]+var_dict[gamma[findfirst(x->x=="A1$(i)A1$(i+1)", SymG)]]') in HermitianPSDCone())
        @constraint(model, lim+Hermitian(var_dict[gamma[findfirst(x->x=="A1$(i)A1$(i+1)", SymG)]]+var_dict[gamma[findfirst(x->x=="A1$(i)A1$(i+1)", SymG)]]') in HermitianPSDCone())
        @constraint(model, eps_*0.01*I(d)-Hermitian(var_dict[gamma[findfirst(x->x=="A3$(i)A3$(i+1)", SymG)]]+var_dict[gamma[findfirst(x->x=="A3$(i)A3$(i+1)", SymG)]]') in HermitianPSDCone())
        @constraint(model, eps_*0.01*I(d)+Hermitian(var_dict[gamma[findfirst(x->x=="A3$(i)A3$(i+1)", SymG)]]+var_dict[gamma[findfirst(x->x=="A3$(i)A3$(i+1)", SymG)]]') in HermitianPSDCone())
        
        #@constraint(model, tr(var_dict[gamma[findfirst(x->x=="ρ1A2$(i)", SymG)]]) <= tr(var_dict[gamma[findfirst(x->x=="ρ1A3$(i)", SymG)]]) )    
        #@constraint(model, tr(var_dict[gamma[findfirst(x->x=="ρ1A2$(i)", SymG)]]) <= tr(var_dict[gamma[findfirst(x->x=="ρ1A3$(i+1)", SymG)]]) )    

        @constraint(model, Hermitian(var_dict[gamma[findfirst(x->x=="ρ1A1$(i)", SymG)]]+var_dict[gamma[findfirst(x->x=="ρ1A1$(i+1)", SymG)]]- d*var_dict[gamma[findfirst(x->x=="ρ1A2$(i)", SymG)]]) in HermitianPSDCone()  )
    end
    @constraint(model, eps_*0.01*I(d)-Hermitian(var_dict[gamma[findfirst(x->x=="A31A3$(inputs)", SymG)]]+var_dict[gamma[findfirst(x->x=="A31A3$(inputs)", SymG)]]') in HermitianPSDCone())
    @constraint(model, eps_*0.01*I(d)+Hermitian(var_dict[gamma[findfirst(x->x=="A31A3$(inputs)", SymG)]]+var_dict[gamma[findfirst(x->x=="A31A3$(inputs)", SymG)]]') in HermitianPSDCone())
    @constraint(model, Hermitian(var_dict[gamma[findfirst(x->x=="ρ1A11", SymG)]] + var_dict[gamma[findfirst(x->x=="ρ1A1$(inputs)", SymG)]] - d*var_dict[gamma[findfirst(x->x=="ρ1A2$(inputs)", SymG)]]) in HermitianPSDCone()  )    

    #@constraint(model, tr(var_dict[gamma[findfirst(x->x=="ρ1A2$(inputs)", SymG)]]) <= tr(var_dict[gamma[findfirst(x->x=="ρ1A3$(inputs)", SymG)]])  )    
    #@constraint(model, tr(var_dict[gamma[findfirst(x->x=="ρ1A2$(inputs)", SymG)]]) <= tr(var_dict[gamma[findfirst(x->x=="ρ1A31", SymG)]])  )    
    var_dict[gamma[findfirst(x->x=="A1$(inputs)A2$(inputs)", SymG)]] = @variable(model, [1:d,1:d])
    @constraint(model, tr(var_dict[gamma[findfirst(x->x=="ρ1A2$(inputs)", SymG)]])*d == tr(var_dict[gamma[findfirst(x->x=="A1$(inputs)A2$(inputs)", SymG)]]) )
    local lim=Hermitian(eps_*0.01*I(d)*tr(var_dict[gamma[findfirst(x->x=="A1$(inputs)A2$(inputs)", SymG)]]/d))
    @constraint(model, I - Hermitian(var_dict[gamma[findfirst(x->x=="A1$(inputs)A2$(inputs)", SymG)]]) in HermitianPSDCone())
    @constraint(model, lim-Hermitian(var_dict[gamma[findfirst(x->x=="A11A1$(inputs)", SymG)]]+var_dict[gamma[findfirst(x->x=="A11A1$(inputs)", SymG)]]') in HermitianPSDCone())
    @constraint(model, lim+Hermitian(var_dict[gamma[findfirst(x->x=="A11A1$(inputs)", SymG)]]+var_dict[gamma[findfirst(x->x=="A11A1$(inputs)", SymG)]]') in HermitianPSDCone())

    model = trace_constraints(model, variables, var_dict)
    set_optimizer(model, Mosek.Optimizer)
    optimize!(model)
    push!(val, [eps_*0.01, objective_value(model)])
end

#writedlm("uncertainty_N$(inputs)_d$(d).txt", val)