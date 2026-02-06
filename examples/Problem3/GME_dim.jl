using LinearAlgebra
using SparseArrays
using JuMP
using MosekTools
using Combinatorics
using Ket 
using MathOptInterface
using Hypatia
using SCS
MOI = MathOptInterface

# ----------------------------------------------------------
# Helper: bipartitions(v) — returns minimal non-repeating bipartitions of v
function bipartitions(v::Vector{Int})
    n = length(v)
    parts = Vector{Tuple{Vector{Int},Vector{Int}}}()
    for mask = 1:(2^(n-1)-1)  # skip empty and full
        part1 = Int[]
        part2 = Int[]
        for i = 1:n
            if ((mask >> (i-1)) & 1) == 1
                push!(part1, v[i])
            else
                push!(part2, v[i])
            end
        end
        if length(part1) <= length(part2)
            push!(parts, (part1, part2))
        else
            push!(parts, (part2, part1))
        end
    end
    return parts
end

# ----------------------------------------------------------
# Dicke state generator (kept structurally identical to MATLAB version)
function dicke_state_d_dim(k::Vector{Int})
    n = sum(k)
    d = length(k)

    arr = Int[]
    for val = 1:d
        append!(arr, fill(val-1, k[val]))
    end

    P = unique(collect(permutations(arr)))

    dim = d^n
    DickeState = zeros(Float64, dim)

    for p in P
        idx = 1
        for position = 1:n
            idx = (idx - 1) * d + p[position] + 1
        end
        DickeState[idx] += 1.0
    end

    DickeState ./= sqrt(length(P))
    return DickeState
end

v = [1, 2, 3, 4, 5, 6]          # Partitions
d = 4                           # Local dimension
r = d - 1                       # SR
n = length(v)                   # Number of qudits
dim = d^n                       # Dimension of the global space

id = I(d)
comp = Vector{Vector{Float64}}(undef, d)
for l = 0:(d-1)
    comp[l+1] = id[:, l+1]
end

# Use Ket.gellmann instead of GenPauli
X = gellmann(2, 1, d)
Z = gellmann(1, 2, d)

k = [2, 2, 1, 1]

test_results = []

for n_test = 1:1
    ket = sparse(dicke_state_d_dim(k))
    rho_pure = ket * ket'
    # Visibility (JuMP variable defined below)
    Id = I(dim)
    biparts = bipartitions(v)

    model = Model(SCS.Optimizer)

    @variable(model, 0 <= vis <= 1)


    Xvars = Vector{Any}(undef, length(biparts))

    for i in eachindex(biparts)
        part1 = biparts[i][1]
        part2 = biparts[i][2]
        kl = length(part1)
        dim_k = d^kl

        Xmat = @variable(model, [1:dim_k, 1:dim_k] in HermitianPSDCone())
        Xvars[i] = Xmat

        @constraint(model, real(tr(Xmat)) <= r / dim_k)
        tmp = @variable(model, [1:dim_k, 1:dim_k] in HermitianPSDCone())
        @constraint(model, Hermitian(I(dim_k) - Xmat) .== tmp)
    end

    # sum_i trace(X_i) <= r
    @constraint(model, real(sum(tr(Xvars[i]) for i in eachindex(Xvars))) <= r)

    # Build M_vars (extensions and permutations)
    M_vars = Vector{Any}(undef, length(Xvars))
    dims_vec = fill(d, n)  # [d, d, ..., d]

    for i in eachindex(Xvars)
        part1 = biparts[i][1]
        part2 = biparts[i][2]
        local kk = length(part1)
        dim_k = d^kk
        dim_comp = d^(length(part2))
        I_comp = Matrix(I, dim_comp, dim_comp)

        Xi_ext = LinearAlgebra.kron(Xvars[i], I_comp)

        ordered_sites = vcat(part1, part2)
        # Use permute_systems instead of PermuteSystems
        M_vars[i] = permute_systems(Xi_ext, ordered_sites, dims_vec)
    end

    # M_sum = sum_i M_vars{i}
    M_sum = sum(M_vars)
    rho_mixed = rho_pure.*vis + Id.*((1 - vis)/dim)

    # Constraints mirroring MATLAB
    @constraint(model, Hermitian(M_sum - rho_mixed) in HermitianPSDCone())
    @constraint(model, Hermitian(Matrix(I, dim, dim) - M_sum) in HermitianPSDCone())

    @objective(model, Min, -vis)

    optimize!(model)
    push!(test_results, [d, n, objective_value(model)])

end

