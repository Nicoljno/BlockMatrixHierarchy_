function theta_identity(Θ::Union{Symbol,AbstractString}, d::Int64, model::JuMP.Model)
    s = Symbol(Θ)
    if s === :Full
        return Matrix{Float64}(I, d, d)               # I(d)
    elseif s === :DirectSum || s === :Direct_sum
        T = JuMP.GenericAffExpr{Float64,JuMP.VariableRef}
        M = zeros(T, d, d)                             # holds numbers & variables
        for i in 1:d-1
            M[i,i] = 1.0                               # I(d-1)
        end
        M[d,d] = @variable(model)                      # the d_ scalar (real)
        return M
    elseif s === :Trace
        return @variable(model)
    else
        throw(ArgumentError("Unknown Θ=$(Θ). Use :Full or :DirectSum."))
    end
end

function embed(matrix::Any, D::Int64)
    new_M = zeros((D, D))
    new_M[1:size(matrix, 1), 1:size(matrix, 1)].=matrix
    return new_M
end

function embed(matrix::Matrix{Float64}, D::Int64)
    new_M = zeros((D, D))
    new_M[1:size(matrix, 1), 1:size(matrix, 1)].=matrix
    return new_M
end

function embed(matrix::Matrix{ComplexF64}, D::Int64)
    new_M = Complex.(zeros((D, D)))
    new_M[1:size(matrix, 1), 1:size(matrix, 1)].=matrix
    return new_M
end

function embed(H::Hermitian{ComplexF64,<:AbstractMatrix{ComplexF64}}, D::Int64)
    n = size(H, 1)
    A = zeros(ComplexF64, D, D)
    A[1:n, 1:n] .= Matrix(H)
    return Hermitian(A)
end

export embed

function find_in_dict(monom::Vector{SymbolicMonomial}, var_dict::Dict{Vector{SymbolicMonomial}, Matrix})::Tuple{Vector{SymbolicMonomial}, Bool}
    n=length(monom)-1
    for _ in 1:n
        monom = reduce_monomial(vcat(monom[end], monom[1:n]))
        if monom in keys(var_dict)
            reduced = monom
            return reduced, true
        end
    end   
    return monom, false
end

function find_in_dict(monom::Vector{SymbolicMonomial}, var_dict::Dict{Vector{SymbolicMonomial}, Matrix}, remove::Int64)::Tuple{Vector{SymbolicMonomial}, Bool}
    n=length(monom)-1
    if monom in keys(var_dict)
        reduced = monom
        return reduced, true
    end

    for _ in 1:n
        if monom[end].party == remove
            monom = reduce_monomial(vcat(monom[end], monom[1:n]))
        end

        if monom in keys(var_dict)
            reduced = monom
            return reduced, true
        end
    end   
    return monom, false
end

function find_in_dict(monom::Vector{SymbolicMonomial}, var_dict::Dict{Vector{SymbolicMonomial}, VariableRef})::Tuple{Vector{SymbolicMonomial}, Bool}
    n=length(monom)-1
    for _ in 1:n
        monom = reduce_monomial(vcat(monom[end], monom[1:n]))
        if monom in keys(var_dict)
            reduced = monom
            return reduced, true
        end
    end   
    return monom, false
end


function find_in_dict(monom::Vector{SymbolicMonomial}, var_dict::Dict{Vector{SymbolicMonomial}, VariableRef}, remove::Int64)::Tuple{Vector{SymbolicMonomial}, Bool}
    n=length(monom)-1
    if monom in keys(var_dict)
        reduced = monom
        return reduced, true
    end

    for _ in 1:n
        if monom[end].party == remove
            monom = reduce_monomial(vcat(monom[end], monom[1:n]))
        end

        if monom in keys(var_dict)
            reduced = monom
            return reduced, true
        end
    end   
    return monom, false
end
export find_in_dict

function numerate_parties(monom::Vector{SymbolicMonomial})
    #parties = [m.party for m in monom]
    parties = [getfield(m, :party) for m in Base.vec(monom) if hasfield(typeof(m), :party)]
    parties = unique(parties)
    if -1 in parties || length(parties) == 2
        return -1
    elseif 1 in parties
        return 1
    elseif 2 in parties
        return 2
    else
        return -1
    end
end
export numerate_parties