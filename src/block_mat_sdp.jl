function BlockMatSDP(d::Int64, 
    loc_d::Int64, 
    num_parties::Int8, 
    num_states::Int8, 
    num_states_parties::Int8,
    inputs::Int8, 
    outputs::Int8, 
    level::Int8, 
    eq::Vector{String}, 
    ineq::Vector{String}, 
    obj::String, 
    state_kind::Int8, 
    element_kind::Int8;
    T = GenericAffExpr{ComplexF64, VariableRef},
    partial_trace_level = nothing,
    Θ::Union{Symbol,AbstractString}  = :Full,
    central_block = [NONE_ELEMENT],
    custom_monomials = Vector{Vector{SymbolicMonomial}}(),
    )
    
    if partial_trace_level===nothing
        if num_states_parties > 0
            all_states = build_all_states(state_kind, num_states, num_states_parties)
        else
            all_states = build_all_states(state_kind, num_states)
        end
        
        if inputs>0
            all_measurements = build_all_measurements(num_parties, outputs, inputs, element_kind)
            monomials = vcat([NONE_ELEMENT], all_states, all_measurements)
        else
            monomials = vcat([NONE_ELEMENT], all_states)
        end

        @assert level≥-1
        monomial_list_level1 = build_monomial_list(Int8(1), monomials; custom_monomials)
        if level == -1
            monomial_list = build_monomial_list_1_plus(monomials; custom_monomials)
        else
            monomial_list = build_monomial_list(level, monomials; custom_monomials)
        end
        gamma = build_gamma(monomial_list) ::Matrix{Vector{SymbolicMonomial}}
        size_of_gamma = size(gamma, 1) ::Int
        variables = unique([gamma[i, j] for i in 1:size_of_gamma for j in i:size_of_gamma]) ::Vector{Vector{SymbolicMonomial}}
        blocks = build_boolean_indicator_matrices(gamma, variables)
        if d == loc_d
            model, var_dict = build_complex_matrices(d, variables, monomial_list_level1; Θ)
        else 
            model, var_dict = build_complex_matrices(d, loc_d, variables, monomial_list_level1; Θ)
        end

        if element_kind!=0
            deleteat!(variables, 1)
        end

        
        G = zeros(T, length(monomial_list)*d, length(monomial_list)*d)
        
        @inbounds for v in variables
            # --- strict upper triangle: triu(blocks[v], 1) ---
            for j in 1:size_of_gamma
                base_j = (j - 1) * d + 1
                for i in 1:j
                    if blocks[v][i, j]
                        base_i = (i - 1) * d + 1
                        G[base_i:base_i+d-1, base_j:base_j+d-1] .= var_dict[v]
                        G[base_j:base_j+d-1, base_i:base_i+d-1] .= var_dict[v]'
                    end
                end
            end

            for i in 1:size_of_gamma
                if blocks[v][i, i]
                    base = (i - 1) * d + 1
                    G[base:base+d-1, base:base+d-1] .= var_dict[v]
                end
            end
        end
        G[1:d, 1:d] .= theta_identity(Θ, d, model)

        G = Hermitian(G, :U)
        println("Constraint on G done")
        @constraint(model, G in HermitianPSDCone())

        model = build_constraint(model, var_dict, gamma, eq, ineq)
        model = trace_constraints(model, variables, var_dict)
        if num_parties == 2
            model = partial_trace_constraints(model, variables, var_dict)
        end

        opt_obj = real(tr(build_objective(model, var_dict, gamma, obj)))
        @objective(model, Max, opt_obj)

        return model, var_dict, variables, G, gamma
    else
        if num_states_parties > 0
            all_states = build_all_states(state_kind, num_states, num_states_parties)
        else
            all_states = build_all_states(state_kind, num_states)
        end
        
        all_measurements = build_all_measurements(num_parties, outputs, inputs, element_kind)

        monomials = vcat([NONE_ELEMENT], all_states, all_measurements)
        @assert level≥-1
        @assert partial_trace_level≥-1
        monomial_list_level1 = build_monomial_list(Int8(1), monomials)
        if level == -1
            monomial_list = build_monomial_list_1_plus(monomials)
        else
            monomial_list = build_monomial_list(level, monomials)
        end
        gamma = build_gamma(monomial_list)

        if partial_trace_level == -1
            partial_trace_monomial_list = build_monomial_list_1_plus(monomials)
        else
            partial_trace_monomial_list = build_monomial_list(partial_trace_level, monomials)
        end

        partial_trace_gamma_1 = build_gamma(partial_trace_monomial_list, 1)
        partial_trace_gamma_2 = build_gamma(partial_trace_monomial_list, 2)
        
        variables = unique([gamma[i, j] for i in 1:size(gamma)[1] for j in i:size(gamma)[1]])
        partial_trace_variables_1 = unique([partial_trace_gamma_1[i, j] for i in 1:size(partial_trace_gamma_1)[1] for j in i:size(partial_trace_gamma_1)[1]])
        partial_trace_variables_2 = unique([partial_trace_gamma_2[i, j] for i in 1:size(partial_trace_gamma_2)[1] for j in i:size(partial_trace_gamma_2)[1]])

        blocks = build_boolean_indicator_matrices(gamma, variables)
        partial_trace_blocks_1 = build_boolean_indicator_matrices(partial_trace_gamma_1, partial_trace_variables_1)
        partial_trace_blocks_2 = build_boolean_indicator_matrices(partial_trace_gamma_2, partial_trace_variables_2)

        if d == loc_d
            model, var_dict = build_complex_matrices(d, variables, monomial_list_level1)
        else 
            model, var_dict = build_complex_matrices(d, loc_d, variables, monomial_list_level1)
        end

        model, pt_var_dict_1 = build_complex_matrices(model, d, loc_d, partial_trace_variables_1, monomial_list_level1, 1)
        model, pt_var_dict_2 = build_complex_matrices(model, d, loc_d, partial_trace_variables_2, monomial_list_level1, 2)


        G = sum(LinearAlgebra.kron(triu(blocks[v], 1), var_dict[v]) for v in variables) + sum(LinearAlgebra.kron(triu(blocks[v], 1)', var_dict[v]') for v in variables) + sum(LinearAlgebra.kron(Diagonal(blocks[v]), var_dict[v]) for v in variables)
        G[1:d, 1:d] = I(d)

        pt_G_1 = sum(LinearAlgebra.kron(triu(partial_trace_blocks_1[v], 1), pt_var_dict_1[v]) for v in partial_trace_variables_1) + sum(LinearAlgebra.kron(triu(partial_trace_blocks_1[v], 1)', pt_var_dict_1[v]') for v in partial_trace_variables_1) + sum(LinearAlgebra.kron(Diagonal(partial_trace_blocks_1[v]), pt_var_dict_1[v]) for v in partial_trace_variables_1)
        pt_G_1[1:loc_d, 1:loc_d] = loc_d.*I(loc_d)

        pt_G_2 = sum(LinearAlgebra.kron(triu(partial_trace_blocks_2[v], 1), pt_var_dict_2[v]) for v in partial_trace_variables_2) + sum(LinearAlgebra.kron(triu(partial_trace_blocks_2[v], 1)', pt_var_dict_2[v]') for v in partial_trace_variables_2) + sum(LinearAlgebra.kron(Diagonal(partial_trace_blocks_2[v]), pt_var_dict_2[v]) for v in partial_trace_variables_2)
        pt_G_2[1:loc_d, 1:loc_d] = loc_d.*I(loc_d)

        for v in variables
            if v in keys(pt_var_dict_1)
                @constraint(model, partial_trace(var_dict[v],1) == pt_var_dict_1[v])
            elseif reverse(v) in keys(pt_var_dict_1)
                @constraint(model, partial_trace(var_dict[v]',1) == pt_var_dict_1[reverse(v)])
            end
            if v in keys(pt_var_dict_2)
                @constraint(model, partial_trace(var_dict[v],2) == pt_var_dict_2[v])
            elseif reverse(v) in keys(pt_var_dict_2)
                @constraint(model, partial_trace(var_dict[v]',2) == pt_var_dict_2[reverse(v)])
            end
        end

        
        @constraint(model, Hermitian(G) in HermitianPSDCone())
        @constraint(model, Hermitian(pt_G_1) in HermitianPSDCone())
        @constraint(model, Hermitian(pt_G_2) in HermitianPSDCone())

        model = build_constraint(model, var_dict, gamma, eq, ineq)
        model = trace_constraints(model, partial_trace_variables_1, pt_var_dict_1)

        if num_parties > 1
            model = partial_trace_constraints(model, variables, var_dict)
        end

        opt_obj = real(tr(build_objective(model, pt_var_dict_1, partial_trace_gamma_1, obj)))

        @objective(model, Max, opt_obj)

        return model, var_dict, pt_var_dict_1, partial_trace_variables_1, variables, G, gamma, pt_G_1, partial_trace_gamma_1
    end
end

function BlockMatSDP_real(d::Int64, 
    loc_d::Int64, 
    num_parties::Int8, 
    num_states::Int8, 
    num_states_parties::Int8,
    inputs::Int8, 
    outputs::Int8, 
    level::Int8, 
    eq::Vector{String}, 
    ineq::Vector{String}, 
    obj::String, 
    state_kind::Int8, 
    element_kind::Int8;
    T = GenericAffExpr{ComplexF64, VariableRef},
    partial_trace_level = nothing,
    Θ::Union{Symbol,AbstractString}  = :Full,
    central_block = [NONE_ELEMENT],
    custom_monomials = Vector{Vector{SymbolicMonomial}}(),
    )
    
    if num_states_parties > 0
        all_states = build_all_states(state_kind, num_states, num_states_parties)
    else
        all_states = build_all_states(state_kind, num_states)
    end
    
    all_measurements = build_all_measurements(num_parties, outputs, inputs, element_kind)

    monomials = vcat([NONE_ELEMENT], all_states, all_measurements)
    @assert level≥-1
    monomial_list_level1 = build_monomial_list(Int8(1), monomials; custom_monomials)
    if level == -1
        monomial_list = build_monomial_list_1_plus(monomials; custom_monomials)
    else
        monomial_list = build_monomial_list(level, monomials; custom_monomials)
    end
    gamma = build_gamma(monomial_list)
    size_of_gamma = size(gamma, 1) ::Int
    
    variables = unique([gamma[i, j] for i in 1:size(gamma)[1] for j in i:size(gamma)[1]])
    blocks = build_boolean_indicator_matrices(gamma, variables)
    if d == loc_d
        model, var_dict = build_real_matrices(d, variables, monomial_list_level1)
    else 
        model, var_dict = build_real_matrices(d, loc_d, variables, monomial_list_level1)
    end

    T = GenericAffExpr{Float64, VariableRef}
    G = zeros(T, length(monomial_list)*d, length(monomial_list)*d)
    @inbounds for v in variables
        # --- strict upper triangle: triu(blocks[v], 1) ---
        for j in 1:size_of_gamma
            base_j = (j - 1) * d + 1
            for i in 1:j
                if blocks[v][i, j]
                    base_i = (i - 1) * d + 1
                    G[base_i:base_i+d-1, base_j:base_j+d-1] .= var_dict[v]
                    G[base_j:base_j+d-1, base_i:base_i+d-1] .= var_dict[v]'
                end
            end
        end

        for i in 1:size_of_gamma
            if blocks[v][i, i]
                base = (i - 1) * d + 1
                G[base:base+d-1, base:base+d-1] .= var_dict[v]
            end
        end
    end
    G[1:d, 1:d] = I(d)
    G = Hermitian(G, :U)
    println("Constraint on G done")
    @constraint(model, G in HermitianPSDCone())

    model = build_constraint(model, var_dict, gamma, eq, ineq)
    model = trace_constraints(model, variables, var_dict)
    #if num_parties > 1
    #    model = partial_trace_constraints(model, variables, var_dict)
    #end

    opt_obj = real(tr(build_objective(model, var_dict, gamma, obj)))
    @objective(model, Max, opt_obj)

    return model, var_dict, variables, G, gamma
end

# Outdated
function BlockMatSDP_loc(d::Int8, loc_d::Int64, num_parties::Int8, num_states::Int8, inputs::Int8, outputs::Int8, level::Int8, eq::Vector{String}, ineq::Vector{String}, obj::String)   
    all_states = build_all_states(Int8(1), num_states)
    tau = SymbolicState(Int8(-1), Int8(loc_d+1), Int8(-1))
    all_states=vcat(all_states, tau)
    all_partial_states = build_all_states(Int8(-1), num_states)
    all_measurements = build_all_measurements(num_parties, outputs, inputs, Int8(1))
    all_partial_measurements = build_all_measurements(num_parties, outputs, inputs, Int8(1))

    
    monomials = vcat([NONE_ELEMENT], all_states, all_measurements)
    partial_monomials = vcat([NONE_ELEMENT], all_partial_states, all_partial_measurements)

    @assert level≥-1
    monomial_list_level1 = build_monomial_list(Int8(1), monomials)
    if level == -1
        monomial_list = build_monomial_list_1_plus(monomials)
        partial_monomial_list = build_monomial_list_1_plus(partial_monomials)
    else
        monomial_list = build_monomial_list(level, monomials)
        partial_monomial_list = build_monomial_list(level, partial_monomials)
    end
    gamma = build_gamma(monomial_list)
    partial_gamma = build_gamma(partial_monomial_list)
    for i in eachindex(partial_gamma[:,1])
        partial_gamma[1, i] = gamma[1, i]
        partial_gamma[i, 1] = gamma[i, 1]
    end
    
    variables = unique([gamma[i, j] for i in 1:size(gamma)[1] for j in i:size(gamma)[1]])
    partial_variables = unique([partial_gamma[i, j] for i in 1:size(partial_gamma)[1] for j in i:size(partial_gamma)[1]])
    add_partial_variables = [v for v in partial_variables if v ∉ variables]
    
    blocks = build_boolean_indicator_matrices(gamma, variables)
    partial_blocks = build_boolean_indicator_matrices(partial_gamma, partial_variables)
    
    model, var_dict = build_complex_matrices(d, variables, monomial_list_level1)
    model, var_dict = add_complex_matrices(model, d, add_partial_variables, monomial_list_level1, var_dict)
    
    G = Matrix{GenericAffExpr{ComplexF64, VariableRef}}(undef, length(monomial_list)*d, length(monomial_list)*d) 

    @inbounds for v in variables
        # --- strict upper triangle: triu(blocks[v], 1) ---
        for j in 1:size_of_gamma
            base_j = (j - 1) * d + 1
            for i in 1:j
                if blocks[v][i, j]
                    base_i = (i - 1) * d + 1
                    G[base_i:base_i+d-1, base_j:base_j+d-1] .= var_dict[v]
                    G[base_j:base_j+d-1, base_i:base_i+d-1] .= var_dict[v]'
                end
            end
        end

        for i in 1:size_of_gamma
            if blocks[v][i, i]
                base = (i - 1) * d + 1
                G[base:base+d-1, base:base+d-1] .= var_dict[v]
            end
        end
    end
    G[1:d, 1:d] = I(d)

    partial_G = []
    push!(partial_G, sum(LinearAlgebra.kron(triu(partial_blocks[v], 1), partial_trace(var_dict[v], 1)) for v in partial_variables) +
        sum(LinearAlgebra.kron(triu(partial_blocks[v], 1)', partial_trace(var_dict[v], 1)') for v in partial_variables) + 
        sum(LinearAlgebra.kron(Diagonal(partial_blocks[v]), partial_trace(var_dict[v], 1)) for v in partial_variables)
    )
    push!(partial_G, sum(LinearAlgebra.kron(triu(partial_blocks[v], 1), partial_trace(var_dict[v], 2)) for v in partial_variables) +
        sum(LinearAlgebra.kron(triu(partial_blocks[v], 1)', partial_trace(var_dict[v], 2)') for v in partial_variables) + 
        sum(LinearAlgebra.kron(Diagonal(partial_blocks[v]), partial_trace(var_dict[v], 2)) for v in partial_variables)
    )
    partial_G[1][1:loc_d,1:loc_d] = I(loc_d)
    partial_G[2][1:loc_d,1:loc_d] = I(loc_d)

    
    model = build_constraint(model, var_dict, gamma, eq, ineq)
    model = trace_constraints(model, variables, var_dict)
    opt_obj = real(tr(build_objective(model, var_dict, gamma, obj)))
    @objective(model, Max, opt_obj)
    return model, var_dict, variables, partial_variables, G, gamma, partial_G, partial_gamma
end



export BlockMatSDP, BlockMatSDP_loc, BlockMatSDP_real