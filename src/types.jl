abstract type SymbolicMonomial end
export SymbolicMonomial

struct SymbolicElement <: SymbolicMonomial
    party::Int8    # e.g. 1, 2,... or -1 if not applicable
    output::Int8   # e.g. 1, 2,... or -1 if not applicable
    input::Int8    # e.g. 1, 2,... or -1 if not applicable
    kind::Int8     #  1 => projector, 0 => observable, -1 => none
end
export SymbolicElement

struct SymbolicState <: SymbolicMonomial
    kind::Int8     #  1 => pure, 0 => mixed, -1 => none
    number::Int8   #  e.g. 1, 2,... or -1 if not applicable
    party::Int8    #  -1 => shared state, 1, 2, ... if applicable
end
export SymbolicState

struct SymbolicScalar <: SymbolicMonomial
    number::Int8   #  e.g. 1, 2,... or -1 if not applicable
end
export SymbolicScalar

struct SymbolicNone <: SymbolicMonomial
end
export SymbolicNone

const NONE_ELEMENT = SymbolicNone()

function isless_monom(monom1::Vector{SymbolicMonomial}, monom2::Vector{SymbolicMonomial})::Bool
    n1 = length(monom1)
    n2 = length(monom2)
    # First check if a monomial is longer than the other
    if n1<n2
        return true
    elseif n2>n1
        return false
    end
    # For equal lengths, compare element by element.
    for (x, y) in zip(monom1, monom2)
        # If the elements are identical (by reference), continue.
        if x === y
            continue
        end

        # Define a numeric code for each type: SymbolicNone=0, SymbolicState=1, SymbolicElement=2.
        t_x = x isa SymbolicNone ? 0 : (x isa SymbolicState ? 1 : 2)
        t_y = y isa SymbolicNone ? 0 : (y isa SymbolicState ? 1 : 2)

        # If the types differ, decide by the numeric code.
        if t_x != t_y
            return t_x < t_y
        end

        # If both are SymbolicState, compare by kind then by number.
        if x isa SymbolicState && y isa SymbolicState
            if x.kind != y.kind
                return x.kind < y.kind
            elseif x.number != y.number
                return x.number < y.number
            elseif x.party != y.party
                return x.party < y.party
            else
                continue
            end
        end

        # If both are SymbolicElement, compare by party, input, output, and then kind.
        if x isa SymbolicElement && y isa SymbolicElement
            if x.party != y.party
                return x.party < y.party
            elseif x.input != y.input
                return x.input < y.input
            elseif x.output != y.output
                return x.output < y.output
            elseif x.kind != y.kind
                return x.kind < y.kind
            else
                continue
            end
        end

        # For SymbolicNone (or any types we haven't otherwise differentiated), assume equality.
    end
    
    # If all elements compare equal, then monom1 is not less than monom2.
    return false
end

function build_all_states(kind::Int8, num_states::Int8, num_parties::Int8)
    all_states = SymbolicMonomial[]
    for n in 1:num_states
        for n_p in 1:num_parties
            push!(all_states, SymbolicState(kind, Int8(n), Int8(n_p)))
        end
    end
    return all_states
end

function build_all_states(kind::Int8, num_states::Int8)
    all_states = SymbolicMonomial[]
    for n in 1:num_states
        push!(all_states, SymbolicState(kind, Int8(n), Int8(-1)))
    end
    return all_states
end
export build_all_states

function build_all_measurements(parties::Vector{Int8}, inputs::Vector{Int8}, outputs::Vector{Int8}, kind::Vector{Int8})
    all_measurements = SymbolicMonomial[]
    for p in parties
        for i::Int8=1:inputs[p]
            for o::Int8=1:outputs[p]
                push!(all_measurements, SymbolicElement(p, i, o, kind[p]))
            end
        end
    end
    return all_measurements
end

function build_all_measurements(parties::Int8, inputs::Int8, outputs::Int8, kind::Int8)
    all_measurements = SymbolicMonomial[]
    for p::Int8 = 1:parties
        for i::Int8=1:inputs
            for o::Int8=1:outputs
                push!(all_measurements, SymbolicElement(p, i, o, kind))
            end
        end
    end
    return all_measurements
end
export build_all_measurements

function build_monomial_list(level::Int8, monomials::Vector{SymbolicMonomial}; custom_monomials::Vector{Vector{SymbolicMonomial}} = Vector{Vector{SymbolicMonomial}}())
    @assert level ≥ 1 
    cartesian = IterTools.product((monomials for _ in 1:level) ...)
    monomial_list = Vector{Vector{SymbolicMonomial}}()
    for tup in cartesian 
        combo = collect(SymbolicMonomial, tup)  # combo has to be compatible with reduce_monomial
        filter!(!=(NONE_ELEMENT), combo)        # filtering NONE_ELEMENT
        combo = reduce_monomial(combo)          # apply reduction rules
        push!(monomial_list, combo)
    end

    if custom_monomials != Vector{Vector{SymbolicMonomial}}()
        for cm in custom_monomials
            push!(monomial_list, cm)
        end
    end
    
    monomial_list = unique(monomial_list)       # remove repetitions
    return sort(monomial_list, by = length)
end

function build_monomial_list(level::Int8, monomials::Vector{SymbolicMonomial}, remove::Int64)
    @assert level ≥ 1 
    cartesian = IterTools.product((monomials for _ in 1:level) ...)
    monomial_list = Vector{Vector{SymbolicMonomial}}()
    for tup in cartesian 
        combo = collect(SymbolicMonomial, tup)  # combo has to be compatible with reduce_monomial
        filter!(!=(NONE_ELEMENT), combo)        # filtering NONE_ELEMENT
        combo = partial_trace_reduce_monomial(combo, remove)          # apply reduction rules
        push!(monomial_list, combo)
    end
    
    monomial_list = unique(monomial_list)       # remove repetitions
    return sort(monomial_list, by = length)
end
export build_monomial_list

function build_monomial_list_1_plus(monomials::Vector{SymbolicMonomial}; custom_monomials::Vector{Vector{SymbolicMonomial}})
    cartesian = IterTools.product((monomials for _ in 1:1) ...)
    monomial_list = Vector{Vector{SymbolicMonomial}}()
    #"""
    for tup in cartesian 
        combo = collect(SymbolicMonomial, tup)  # combo has to be compatible with reduce_monomial
        filter!(!=(NONE_ELEMENT), combo)        # filtering NONE_ELEMENT
        combo = reduce_monomial(combo)          # apply reduction rules
        push!(monomial_list, combo)
    end

    for cm in custom_monomials
        push!(monomial_list, cm)
    end
    #"""
    # Additionally add monomials of the form SymbolicState + SymbolicElement.
    #"""
    for s in monomials
        if s isa SymbolicState
            for e in monomials
                if e isa SymbolicElement
                    push!(monomial_list, reduce_monomial([s, e]))
                end
            end
        end
    end
    #"""
    # Additionally add monomials of the form SymbolicElement + SymbolicState.
    """
    for s in monomials
        if s isa SymbolicState
            for e in monomials
                if e isa SymbolicElement
                    push!(monomial_list, reduce_monomial([e, s]))
                end
            end
        end
    end
    #"""


    # Additionally add monomials of the form SymbolicState + SymbolicState.
    #"""
    for s in monomials
        if s isa SymbolicState
            for e in monomials
                if e isa SymbolicState
                    push!(monomial_list, reduce_monomial(SymbolicMonomial[e, s]))
                end
            end
        end
    end
    #"""
    
    # Add monomials of the form SymbolicElement + SymbolicElement
    """
    for s in monomials
        if s isa SymbolicElement
            for e in monomials
                if e isa SymbolicElement
                    push!(monomial_list, reduce_monomial(SymbolicMonomial[e, s]))
                end
            end
        end
    end
    #"""
    monomial_list = unique(monomial_list)       # remove repetitions
    println(length(monomial_list))
    return sort(monomial_list, by = length)
end

function build_monomial_list_1_plus(monomials::Vector{SymbolicMonomial}, pers::Int8)
    cartesian = IterTools.product((monomials for _ in 1:1) ...)
    monomial_list = Vector{Vector{SymbolicMonomial}}()
    for tup in cartesian 
        combo = collect(SymbolicMonomial, tup)  # combo has to be compatible with reduce_monomial
        filter!(!=(NONE_ELEMENT), combo)        # filtering NONE_ELEMENT
        combo = reduce_monomial(combo)          # apply reduction rules
        push!(monomial_list, combo)
        #println(combo)
    end
    
    #push!(monomial_list, SymbolicMonomial[])
    #push!(monomial_list, SymbolicMonomial[SymbolicState(1, 1, 1)])
    #push!(monomial_list, SymbolicMonomial[SymbolicState(1, 1, 2)])
    #push!(monomial_list, SymbolicMonomial[SymbolicState(1, 1, 1), SymbolicState(1, 1, 2)])
    #push!(monomial_list, SymbolicMonomial[SymbolicElement(1, 1, 1, 1), SymbolicElement(2, 1, 1, 1)])
    #push!(monomial_list, SymbolicMonomial[SymbolicElement(1, 2, 2, 1), SymbolicElement(2, 2, 2, 1)])
    #push!(monomial_list, SymbolicMonomial[SymbolicElement(1, 0, 0, 1)])
    #push!(monomial_list, SymbolicMonomial[SymbolicElement(2, 0, 0, 1)])

    # Additionally add monomials of the form SymbolicState + SymbolicElement.
    #"""
    for s in monomials
        if s isa SymbolicState
            for e in monomials
                if e isa SymbolicElement #&& e.party != s.party
                    push!(monomial_list, [s,e])
                end
            end
        end
    end
    #"""

    # Additionally add monomials of the form SymbolicElement + SymbolicState.
    """
    for s in monomials
        if s isa SymbolicState
            for e in monomials
                if e isa SymbolicElement && e.party != s.party 
                    if e.input == 1
                        push!(monomial_list, reduce_monomial([e, s]))
                    end
                end
            end
        end
    end
    #"""


    # Additionally add monomials of the form SymbolicState + SymbolicState.
    #"""
    for s in monomials
        if s isa SymbolicState
            for e in monomials
                if e isa SymbolicState
                    push!(monomial_list, reduce_monomial(SymbolicMonomial[e, s]))
                end
            end
        end
    end
    #"""
    
    # Add monomials of the form SymbolicElement + SymbolicElement
    #"""
    for s in monomials
        if s isa SymbolicElement
            for e in monomials
                if e isa SymbolicElement && e.output==s.output && e.input==s.input #&& e.party == s.party
                    push!(monomial_list, reduce_monomial(SymbolicMonomial[e, s]))
                end
            end
        end
    end
    #"""

    monomial_list = unique(monomial_list)       # remove repetitions
    println([show_symbolic(monomial_list[i]) for i in eachindex(monomial_list)])
    return monomial_list
end
export build_monomial_list_1_plus
