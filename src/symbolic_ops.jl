letter_from_party(p::Int8) = Char('A' + (p - 1))
letter_from_state(p::Int8) = Char('a' + (p - 1))
party_from_letter(l::Char) = Int(l)-Int('A')+1
state_from_letter(l::Char) = Int(l)-Int('a')+1

function show_symbolic(x::SymbolicElement)
    # If party == -1 or something else, placeholder
    if x.party > 0
        party_letter = letter_from_party(x.party)
    else
        party_letter = ""  # placeholder
    end

    # Build the output string
    # e.g. "A21" means "A" + "2" + "1"
    return string(party_letter, x.output, x.input)
end

# will have to be updated if the state is not shared between all parties
function show_symbolic(x::SymbolicState)
    if x.party > 0
        party_letter = x.party
    else
        party_letter = ""  # placeholder
    end
    return string("ρ", x.number, party_letter)
end

function show_symbolic(x::SymbolicScalar)
    return string("z", x.number)
end

function show_symbolic(x::SymbolicNone)
    return " "
end

function show_symbolic(m::SymbolicMonomial)
    if m isa SymbolicElement
        return show_symbolic(m::SymbolicElement)
    elseif m isa SymbolicState
        return show_symbolic(m::SymbolicState)
    elseif m isa SymbolicScalar
        return show_symbolic(m::SymbolicScalar)
    else
        error("Unknown subtype of SymbolicMonomial: $(typeof(m))")
    end
end

function show_symbolic(vec::Vector{SymbolicMonomial})
    strs = map(show_symbolic, vec)
    # Concatenate end-to-end.
    return join(strs)
end
export show_symbolic

function symbolic_matrix(gamma::Matrix{Vector{SymbolicMonomial}})
    return [show_symbolic(gamma[i,j]) for i in 1:size(gamma,1), j in 1:size(gamma,2)]
end
export symbolic_matrix

#Complex variables
function eval_expr(expr_string::String, chars::Vector{Symbol}, vals::Vector{Matrix{GenericAffExpr{ComplexF64, VariableRef}}})
    # Parse the string into a Julia expression
    ex = Meta.parse(expr_string)
    assignments = [Expr(:(=), Symbol(chars[i]), vals[i]) for i in eachindex(chars)]
    # Evaluate the expression in a local scope
    return eval(:(let $(assignments...)
            $ex
    end))
end


#Real variables
function eval_expr(expr_string::String, chars::Vector{Symbol}, vals::Vector{Matrix{VariableRef}})
    # Parse the string into a Julia expression
    ex = Meta.parse(expr_string)
    assignments = [Expr(:(=), Symbol(chars[i]), vals[i]) for i in eachindex(chars)]
    # Evaluate the expression in a local scope
    return eval(:(let $(assignments...)
            $ex
    end))
end

#Real scalar variables
function eval_expr(expr_string::String, chars::Vector{Symbol}, vals::Vector{VariableRef})
    # Parse the string into a Julia expression
    ex = Meta.parse(expr_string)
    assignments = [Expr(:(=), Symbol(chars[i]), vals[i]) for i in eachindex(chars)]
    # Evaluate the expression in a local scope
    return eval(:(let $(assignments...)
            $ex
    end))
end

function parse_monomials(s::AbstractString, element_kind::Integer, state_kind::Integer)::Vector{SymbolicMonomial}

    ek = Int8(element_kind)
    sk = Int8(state_kind)

    out = SymbolicMonomial[]

    # Step 1: split into substrings that don't start with a number
    # Equivalent to: each token starts with a non-digit, then has zero or more digits.
    for m in eachmatch(r"[^\d]\d*", s)
        token = m.match
        chars = collect(token)  # safe indexing even with unicode like 'ρ'

        # Step 2: build SymbolicState / SymbolicElement
        if chars[1] == 'ρ'
            length(chars) >= 2 || throw(ArgumentError("State token '$token' missing number"))
            number = Int8(chars[2] - '0')
            push!(out, SymbolicState(sk, number, Int8(-1)))
        elseif chars[1] == 'z'
            length(chars) >= 2 || throw(ArgumentError("State token '$token' missing number"))
            number = Int8(chars[2] - '0')
            push!(out, SymbolicScalar(sk, number))
        else
            length(chars) >= 3 || throw(ArgumentError("Element token '$token' must look like Axy (e.g. A11)"))
            party  = party_from_letter(chars[1])
            output = Int8(chars[2] - '0')
            input  = Int8(chars[3] - '0')
            push!(out, SymbolicElement(party, output, input, ek))
        end
    end

    return out
end

function parse_objective_terms(obj::AbstractString, element_kind::Integer, state_kind::Integer
    )::Vector{Vector{SymbolicMonomial}}

    # remove whitespace and strip a single outer (...) if present
    s = replace(obj, r"\s+" => "")
    if startswith(s, "(")
        # strip only the first matching outer parentheses
        close = findfirst(')', s)
        if close !== nothing
            inside = s[2:close-1]
            s = inside
        end
    end

    # split by '+', ignore empty pieces
    pieces = filter(!isempty, split(s, '+'))

    # parse + unique (preserve first-seen order)
    seen = Set{Vector{SymbolicMonomial}}()
    out  = Vector{Vector{SymbolicMonomial}}()

    for p in pieces
        v = parse_monomials(p, element_kind, state_kind)
        if !(v in seen)
            push!(seen, v)
            push!(out, v)
        end
    end

    return out
end