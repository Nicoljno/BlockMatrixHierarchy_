function build_gamma(monom_list::Vector{Vector{SymbolicMonomial}}) ::Matrix{Vector{SymbolicMonomial}}
    n = length(monom_list)
    gamma = fill(SymbolicMonomial[], n, n)
    for i in eachindex(monom_list), j in eachindex(monom_list)
        m_j_dagger = reverse(monom_list[j])
        gamma[i,j] = vcat(monom_list[i], m_j_dagger)
        gamma[i,j] = reduce_monomial(gamma[i,j])
    end
    return gamma
end

function build_gamma(monom_list::Vector{Vector{SymbolicMonomial}}, remove::Int64) ::Matrix{Vector{SymbolicMonomial}}
    n = length(monom_list)
    gamma = fill(SymbolicMonomial[], n, n)
    for i=1:n
        for j=1:n
            m_j_dagger = reverse(monom_list[j])
            gamma[i,j] = vcat(monom_list[i], m_j_dagger)
            gamma[i,j] = partial_trace_reduce_monomial(gamma[i,j], remove)
        end
    end
    return gamma
end
export build_gamma

function build_localizing_gamma(monom_list::Vector{Vector{SymbolicMonomial}}, loc_element::SymbolicState) ::Matrix{Vector{SymbolicMonomial}}
    n = length(monom_list)
    gamma = fill(SymbolicMonomial[], n, n)
    for i=1:n
        m_i_dagger = reverse(monom_list[i])
        for j=1:n
            gamma[i,j] = vcat(m_i_dagger, loc_element, monom_list[j])
            gamma[i,j] = reduce_monomial(gamma[i,j])
        end
    end
    return gamma
end


function build_boolean_indicator_matrices(gamma::Matrix{Vector{SymbolicMonomial}}, vars::Vector{Vector{SymbolicMonomial}})
    BlockDict = Dict{Vector{SymbolicMonomial}, SparseMatrixCSC{Bool}}()
    for v in vars
        BlockDict[v] = sparse(map(x -> x == v, gamma))
    end
    return BlockDict
end
export build_boolean_indicator_matrices

function build_real_matrices(d::Int64, vars::Vector{Vector{SymbolicMonomial}}, monomial_list_level1::Vector{Vector{SymbolicMonomial}})
    model = Model(Mosek.Optimizer)
    var_dict = Dict{Vector{SymbolicMonomial}, Matrix}()

    for (_, v) in enumerate(vars)
        revred = reduce_monomial(reverse(v))
        if revred == v
            # Create a PSD or symmetric matrix variable.
            if any(v == reduce_monomial(vcat(A, B, reverse(A))) for A in vars for B in monomial_list_level1)
                tmp=@variable(model, [i=1:d, j=1:d] in PSDCone())
            else
                tmp = @variable(model, [i=1:d, j=1:d], Symmetric)
            end
            var_dict[v] = tmp
        else
            # Create a matrix variable.
            tmp = @variable(model, [i=1:d, j=1:d])
            var_dict[v] = tmp
        end
    end

    return model, var_dict
end

function build_real_matrices(d::Int64, loc_d::Int64, vars::Vector{Vector{SymbolicMonomial}}, monomial_list_level1::Vector{Vector{SymbolicMonomial}}; solver=Mosek.Optimizer)
    model = Model(solver)
    var_dict = Dict{Vector{SymbolicMonomial}, Matrix}()
    
    sdp_vars = [item for item in monomial_list_level1 if !any(x->x.kind==-1, item)]

    for (idx, v) in enumerate(vars)

        revred = reduce_monomial(reverse(v))
        party_number = numerate_parties(v)
        if party_number == -1
            if revred == v
                # Create a Hermitian complex matrix variable.
                if any(v == reduce_monomial(vcat(A, B, reverse(A))) for A in vars for B in sdp_vars)
                    tmp=@variable(model, [i=1:d, j=1:d] in PSDCone())
                else
                    tmp = @variable(model, [i=1:d, j=1:d], Symmetric)
                end
                var_dict[v] = tmp
            else
                # Create a complex matrix variable.
                if any(v == reduce_monomial(vcat(A, B, reverse(A))) for A in vars for B in sdp_vars)
                    tmp = @variable(model, [1:d, 1:d] in PSDCone())
                else
                    tmp = @variable(model, [i=1:d, j=1:d])
                end
                var_dict[v] = tmp
            end
        elseif party_number == 1
            if revred == v
                # Create a Hermitian complex matrix variable.
                if any(v == reduce_monomial(vcat(A, B, reverse(A))) for A in vars for B in sdp_vars)
                    tmp=@variable(model, [i=1:loc_d, j=1:loc_d] in PSDCone())
                else
                    tmp = @variable(model, [i=1:loc_d, j=1:loc_d], Symmetric)
                end
                var_dict[v] = LinearAlgebra.kron(tmp, I(loc_d))
            else
                # Create a complex matrix variable.
                if any(v == reduce_monomial(vcat(A, B, reverse(A))) for A in vars for B in sdp_vars)
                    tmp = @variable(model, [1:loc_d, 1:loc_d] in PSDCone())
                else
                    tmp = @variable(model, [i=1:loc_d, j=1:loc_d])
                end
                var_dict[v] = LinearAlgebra.kron(tmp, I(loc_d))
            end
        elseif party_number == 2
            if revred == v
                # Create a Hermitian complex matrix variable.
                if any(v == reduce_monomial(vcat(A, B, reverse(A))) for A in vars for B in sdp_vars)
                    tmp=@variable(model, [i=1:loc_d, j=1:loc_d] in PSDCone())
                else
                    tmp = @variable(model, [i=1:loc_d, j=1:loc_d], Symmetric)
                end
                var_dict[v] = LinearAlgebra.kron(I(loc_d), tmp)
            else
                # Create a complex matrix variable.
                if any(v == reduce_monomial(vcat(A, B, reverse(A))) for A in vars for B in monomial_list_level1)
                    tmp = @variable(model, [1:loc_d, 1:loc_d] in PSDCone())
                else
                    tmp = @variable(model, [i=1:loc_d, j=1:loc_d])
                end
                var_dict[v] = LinearAlgebra.kron(I(loc_d), tmp)
            end
        end
    end

    return model, var_dict
end
export build_real_matrices

function build_real_variables(model::Model, ft_vars::Vector{Vector{SymbolicMonomial}})
    
    ft_var_dict = Dict{Vector{SymbolicMonomial}, JuMP.VariableRef}()
    
    for (idx, v) in enumerate(ft_vars)

        tmp = @variable(model)
        ft_var_dict[v] = tmp
    end
    return model, ft_var_dict
end
export build_real_variables

function build_complex_matrices(
    d::Int64,
    vars::Vector{Vector{SymbolicMonomial}},
    monomial_list_level1::Vector{Vector{SymbolicMonomial}};
    Θ::Union{Symbol,AbstractString} = :Full,
)
    build_complex_matrices(Val(Symbol(Θ)), d, vars, monomial_list_level1)
end

function build_complex_matrices(
    d::Int64, loc_d::Int64,
    vars::Vector{Vector{SymbolicMonomial}},
    monomial_list_level1::Vector{Vector{SymbolicMonomial}};
    Θ::Union{Symbol,AbstractString} = :Full,
)
    build_complex_matrices(Val(Symbol(Θ)), d, loc_d, vars, monomial_list_level1)
end

function build_complex_matrices(::Val{:Full}, d::Int64, vars::Vector{Vector{SymbolicMonomial}}, monomial_list_level1::Vector{Vector{SymbolicMonomial}})
    model = Model(Mosek.Optimizer)
    var_dict = Dict{Vector{SymbolicMonomial}, Matrix}()
    
    for (idx, v) in enumerate(vars)
        revred = reduce_monomial(reverse(v))
        if revred == v
            # Create a Hermitian complex matrix variable.
            if any(v == reduce_monomial(vcat(A, B, reverse(A))) for A in vars for B in monomial_list_level1)
                tmp=@variable(model, [i=1:d, j=1:d] in HermitianPSDCone())
            else
                tmp = @variable(model, [i=1:d, j=1:d], Hermitian)
            end
            var_dict[v] = tmp
        else
            # Create a complex matrix variable.
            tmp = @variable(model, [i=1:d, j=1:d] in ComplexPlane())
            var_dict[v] = tmp
        end
    end

    return model, var_dict
end

function build_complex_matrices(::Val{:DirectSum}, d::Int64, vars::Vector{Vector{SymbolicMonomial}}, monomial_list_level1::Vector{Vector{SymbolicMonomial}})
    model = Model(Mosek.Optimizer)
    T = GenericAffExpr{ComplexF64, VariableRef}
    var_dict = Dict{Vector{SymbolicMonomial}, Matrix}()
    for (_, v) in enumerate(vars)
        revred = reduce_monomial(reverse(v))
        tmp=zeros(T, d, d)
        if revred == v
            # Create a Hermitian complex matrix variable.
            if any(v == reduce_monomial(vcat(A, B, reverse(A))) for A in vars for B in monomial_list_level1)
                tmp[1:d-1,1:d-1] = @variable(model, [i=1:d-1, j=1:d-1] in HermitianPSDCone())
                tmp[d,d] = @variable(model)
            else
                tmp[1:d-1,1:d-1] = @variable(model, [i=1:d-1, j=1:d-1], Hermitian)
                tmp[d,d] = @variable(model)
            end
            var_dict[v] = tmp
        else
            # Create a complex matrix variable.
            tmp[1:d-1,1:d-1] = @variable(model, [i=1:d-1, j=1:d-1] in ComplexPlane())
            tmp[d,d] = @variable(model, set = ComplexPlane())
            var_dict[v] = tmp
        end
    end

    return model, var_dict
end

function build_complex_matrices(::Val{:Full}, d::Int64, loc_d::Int64, vars::Vector{Vector{SymbolicMonomial}}, monomial_list_level1::Vector{Vector{SymbolicMonomial}})
    model = Model(Mosek.Optimizer)
    var_dict = Dict{Vector{SymbolicMonomial}, Matrix}()
    
    sdp_vars = [item for item in monomial_list_level1 if !any(x->x.kind==0, item)]
    bar_d = Int64(d/loc_d)
    for (_, v) in enumerate(vars)

        revred = reduce_monomial(reverse(v))
        party_number = numerate_parties(v)
        if party_number == -1
            if revred == v
                # Create a Hermitian complex matrix variable.
                if any(v == reduce_monomial(vcat(A, B, reverse(A))) for A in vars for B in sdp_vars)
                    tmp = @variable(model, [i=1:d, j=1:d] in HermitianPSDCone())
                else
                    tmp = @variable(model, [i=1:d, j=1:d], Hermitian)
                end
                var_dict[v] = tmp
            else
                # Create a complex matrix variable.
                tmp = @variable(model, [i=1:d, j=1:d] in ComplexPlane())
                var_dict[v] = tmp
            end
        elseif party_number == 1
            if revred == v
                # Create a Hermitian complex matrix variable.
                if any(v == reduce_monomial(vcat(A, B, reverse(A))) for A in vars for B in sdp_vars)
                    tmp = @variable(model, [i=1:loc_d, j=1:loc_d] in HermitianPSDCone())
                else
                    tmp = @variable(model, [i=1:loc_d, j=1:loc_d], Hermitian)
                end
                var_dict[v] = LinearAlgebra.kron(tmp, I(bar_d))
            else
                # Create a complex matrix variable.
                tmp = @variable(model, [i=1:loc_d, j=1:loc_d] in ComplexPlane())
                var_dict[v] = LinearAlgebra.kron(tmp, I(bar_d))
            end
        elseif party_number == 2
            if revred == v
                # Create a Hermitian complex matrix variable.
                if any(v == reduce_monomial(vcat(A, B, reverse(A))) for A in vars for B in sdp_vars)
                    tmp=@variable(model, [i=1:loc_d, j=1:loc_d] in HermitianPSDCone())
                else
                    tmp = @variable(model, [i=1:loc_d, j=1:loc_d], Hermitian)
                end
                var_dict[v] = LinearAlgebra.kron(I(bar_d), tmp)
            else
                # Create a complex matrix variable.
                tmp = @variable(model, [i=1:loc_d, j=1:loc_d] in ComplexPlane())
                var_dict[v] = LinearAlgebra.kron(I(bar_d), tmp)
            end
        end
    end

    return model, var_dict
end

function build_complex_matrices(model::Model, d::Int64, loc_d::Int64, vars::Vector{Vector{SymbolicMonomial}}, monomial_list_level1::Vector{Vector{SymbolicMonomial}}, remove::Int64)
    
    pt_var_dict = Dict{Vector{SymbolicMonomial}, Matrix}()
    
    sdp_vars = [item for item in monomial_list_level1 if !any(x->x.kind==0, item)]
    
    for (idx, v) in enumerate(vars)

        revred = partial_trace_reduce_monomial(reverse(v), remove)
        party_number = numerate_parties(v)
        
        if party_number != remove
            if revred == v
                # Create a Hermitian complex matrix variable.
                if any(v == partial_trace_reduce_monomial(vcat(A, B, reverse(A)), remove) for A in vars for B in sdp_vars)
                    tmp=@variable(model, [i=1:loc_d, j=1:loc_d] in HermitianPSDCone())
                else
                    tmp = @variable(model, [i=1:loc_d, j=1:loc_d], Hermitian)
                end
                pt_var_dict[v] = tmp
            else
                # Create a complex matrix variable.
                tmp = @variable(model, [i=1:loc_d, j=1:loc_d] in ComplexPlane())
                pt_var_dict[v] = tmp
            end
        elseif party_number == remove
            if revred == v
                tmp=@variable(model, [i=1:loc_d, j=1:loc_d], Hermitian)
            else
                tmp=@variable(model, set = ComplexPlane())
            end
            pt_var_dict[v] = tmp.*I(loc_d)
        end
    end

    return model, pt_var_dict
end
export build_complex_matrices


function add_complex_matrices(model::Model, d::Int64, vars::Vector{Vector{SymbolicMonomial}}, monomial_list_level1::Vector{Vector{SymbolicMonomial}}, var_dict::Dict{Vector{SymbolicMonomial}, Matrix})
    for (idx, v) in enumerate(vars)
        revred = reduce_monomial(reverse(v))
        if revred == v
            # Create a Hermitian complex matrix variable.
            if any(v == reduce_monomial(vcat(A, B, reverse(A))) for A in vars for B in monomial_list_level1)
                tmp=@variable(model, [i=1:d, j=1:d] in HermitianPSDCone())
            else
                tmp = @variable(model, [i=1:d, j=1:d], Hermitian)
            end
            var_dict[v] = tmp
        else
            # Create a complex matrix variable.
            if any(v == reduce_monomial(vcat(A, B, reverse(A))) for A in vars for B in monomial_list_level1)
                tmp = @variable(model, [1:d, 1:d] in HermitianPSDCone())
            else
                tmp = @variable(model, [i=1:d, j=1:d] in ComplexPlane())
            end
            var_dict[v] = tmp
        end
    end

    return model, var_dict
end


function build_objective(model::Model, var_dict::Dict{Vector{SymbolicMonomial}, Matrix}, gamma::Matrix{Vector{SymbolicMonomial}}, expr::String)
    expr_clean = replace(expr, r"\s+" => "")
    variables = [el for el in split(expr_clean, r"[[:punct:]]+") if length(el)>1]
    variables = sort(variables, by=length, rev=true)
    vars = [var_dict[gamma[findfirst(x->x==variables[i], symbolic_matrix(gamma))]] for i in eachindex(variables)]
    variables = [Symbol(el) for el in variables]
    return eval_expr(expr_clean, variables, vars)
end

function build_objective(model::Model, ft_var_dict::Dict{Vector{SymbolicMonomial}, VariableRef}, full_trace_gamma::Matrix{Vector{SymbolicMonomial}}, expr::String)
    expr_clean = replace(expr, r"\s+" => "")
    variables = [el for el in split(expr_clean, r"[[:punct:]]+") if length(el)>1]
    variables = sort(variables, by=length, rev=true)
    vars = [ft_var_dict[full_trace_gamma[findfirst(x->x==variables[i], symbolic_matrix(full_trace_gamma))]] for i in eachindex(variables)]
    variables = [Symbol(el) for el in variables]
    return eval_expr(expr_clean, variables, vars)
end


function build_constraint(model::Model, var_dict::Dict{Vector{SymbolicMonomial}, Matrix}, gamma::Matrix{Vector{SymbolicMonomial}}, eq::Vector{String}, ineq::Vector{String})::Model
    
    symbolic_gamma = symbolic_matrix(gamma)
    
    eq = [replace(eq[i], r"\s+" => "") for i in eachindex(eq)]
    eq_variables = [[el for el in split(eq[i], r"[[:punct:]]+") if el in symbolic_gamma] for i in eachindex(eq)]
    eq_variables = [sort(eq_variables[i], by=length, rev=true) for i in eachindex(eq_variables)]
    eq_vars = [[var_dict[gamma[findfirst(x->x==eq_variables[j][i], symbolic_matrix(gamma))]] for i in eachindex(eq_variables[j])] for j in eachindex(eq_variables)]
    eq_variables = [[Symbol(el) for el in eq_variables[i]] for i in eachindex(eq_variables)]

    ineq = [replace(ineq[i], r"\s+" => "") for i in eachindex(ineq)]
    ineq_variables = [[el for el in split(ineq[i], r"[[:punct:]]+") if el in symbolic_gamma] for i in eachindex(ineq)]
    ineq_variables = [sort(ineq_variables[i], by=length, rev=true) for i in eachindex(ineq_variables)]
    ineq_vars = [[var_dict[gamma[findfirst(x->x==ineq_variables[j][i], symbolic_matrix(gamma))]] for i in eachindex(ineq_variables[j])] for j in eachindex(ineq_variables)]
    ineq_variables = [[Symbol(el) for el in ineq_variables[i]] for i in eachindex(ineq_variables)]
    
    [@constraint(model, real(tr(eval_expr(eq[i], eq_variables[i], eq_vars[i]))) == 0 ) for i in eachindex(eq_variables)]
    [@constraint(model, real(tr(eval_expr(ineq[i], ineq_variables[i], ineq_vars[i]))) >= 0 ) for i in eachindex(ineq_variables)]
    return model

end

function build_conjugate_constraints(model::Model, variables::Vector{Vector{SymbolicMonomial}}, var_dict::Dict{Vector{SymbolicMonomial}, Matrix})::Model
    processed = Set()
    for v in variables
        rev_v = reverse(v)
        if rev_v != v && rev_v in variables && !(rev_v in processed)
            @constraint(model, var_dict[v]==var_dict[rev_v]')
            push!(processed, v)
            push!(processed, rev_v)
        end
    end
    return model
end
export build_conjugate_constraints

function trace_constraints(model::Model, variables::Vector{Vector{SymbolicMonomial}}, var_dict::Dict{Vector{SymbolicMonomial}, Matrix})::Model
    trace_variables = [x for x in variables if length(x)>2]
    for (idx, monom) in enumerate(trace_variables)
        reduced_variable = trace_reduce_monomial(monom)
        reduced_variable, exists = find_in_dict(reduced_variable, var_dict)
        if exists
            @constraint(model, tr(var_dict[trace_variables[idx]])==tr(var_dict[reduced_variable]))
        end
    end
    return model
end

function trace_constraints(model::Model, variables::Vector{Vector{SymbolicMonomial}}, var_dict::Dict{Vector{SymbolicMonomial}, VariableRef})::Model
    trace_variables = [x for x in variables if length(x)>2]
    for (idx, monom) in enumerate(trace_variables)
        reduced_variable = trace_reduce_monomial(monom)
        reduced_variable, exists = find_in_dict(reduced_variable, var_dict)
        if exists
            @constraint(model, var_dict[trace_variables[idx]]==var_dict[reduced_variable])
        end
    end
    return model
end
export trace_constraints

function partial_trace_constraints(model::Model, variables::Vector{Vector{SymbolicMonomial}}, var_dict::Dict{Vector{SymbolicMonomial}, Matrix})::Model
    partial_trace_variables = [x for x in variables if length(x)>2]
    for (idx, monom) in enumerate(partial_trace_variables)
        reduced_variable_B = partial_trace_reduce_monomial(monom, 1)
        reduced_variable_A = partial_trace_reduce_monomial(monom, 2)
        reduced_variable_A, exists_A = find_in_dict(reduced_variable_A, var_dict, 2)
        reduced_variable_B, exists_B = find_in_dict(reduced_variable_B, var_dict, 1)        
        if exists_A && monom != reduced_variable_A
            #println([show_symbolic(monom[i]) for i in eachindex(monom)], [show_symbolic(reduced_variable_A[i]) for i in eachindex(reduced_variable_A)], "B")
            @constraint(model, partial_trace(var_dict[monom], 2) == partial_trace(var_dict[reduced_variable_A], 2) )
        elseif exists_B && monom != reduced_variable_B
            #println([show_symbolic(monom[i]) for i in eachindex(monom)], [show_symbolic(reduced_variable_B[i]) for i in eachindex(reduced_variable_B)], "A")
            @constraint(model, partial_trace(var_dict[monom], 1) == partial_trace(var_dict[reduced_variable_B], 1) )
        end
    end
    return model
end
export partial_trace_constraints