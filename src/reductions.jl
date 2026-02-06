function reorder_monomial_by_party(monomial::Vector{SymbolicMonomial}) :: Vector{SymbolicMonomial}
    result = SymbolicMonomial[]    # will hold the reordered monomial
    block  = SymbolicMonomial[]     # accumulate consecutive operators


    for item in monomial
        if item.party == -1
            # 1) Sort the block of operators in ascending .party
            sort!(block, by=x -> x.party)
            # 2) Append them to result
            append!(result, block)
            empty!(block)
            # 3) Now append the state
            push!(result, item)
        elseif item.party > 0
            push!(block, item)
        else 
            error("Unknown item type: $item")
        end
    end

    # Flush any leftover operators at the end
    sort!(block, by=x -> x.party)
    append!(result, block)


    return result
end

function reduce_monomial(monom::Vector{SymbolicMonomial})::Vector{SymbolicMonomial}
    reduced = SymbolicMonomial[]
    
    monom = reorder_monomial_by_party(monom)

    for item in monom
        if isempty(reduced)
            # There's no previous item => just push
            push!(reduced, item)
            continue
        end

        # Compare 'item' with the last item in 'reduced'
        last_item = reduced[end]

        if last_item isa SymbolicElement && item isa SymbolicElement
            # We have consecutive operators => apply rules (1) and (2)
            a = last_item::SymbolicElement
            b = item::SymbolicElement

            # If projectors => apply othogonality and projectivity 
            if a.kind == 1 && b.kind == 1 && a.party == b.party
                # Rule 1: contradictory outputs => entire monomial => [NONE_ELEMENT]
                if a.input == b.input && a.output != b.output 
                    return SymbolicMonomial[]
                end
                
                # Rule 2: projective => if same (party, input, output) and both kind=1 => skip 'b'
                if a.input == b.input && a.output == b.output 
                    # P^2 => P => do not push 'b'
                    continue
                end
                #rule to add if A is invariant under projector
                """
                if a.output == 0 
                    deleteat!(reduced, length(reduced))
                elseif b.output == 0
                    continue
                end
                """

            end

            # If operators FOR NOW ASSUME QUBIT
            if a.kind == 0 && b.kind == 0 && a.party == b.party
                if a.input == b.input 
                    deleteat!(reduced, length(reduced))
                    continue
                end
            end

            # Otherwise, no conflict => push the new operator
            push!(reduced, b)

        elseif last_item isa SymbolicState && item isa SymbolicState && item.party == last_item.party
            # We have consecutive states => apply rule (3)
            a = last_item::SymbolicState
            b = item::SymbolicState

            # If both kind=1 => skip the second
            if a.kind == 1 && b.kind == 1 && a.number == b.number
                # S^2 => S => do not push 'b'
                continue
            end

            # Otherwise, push the new state
            push!(reduced, b)

        elseif last_item isa SymbolicState && item isa SymbolicElement && item.party == last_item.party
            a = last_item::SymbolicState
            b = item::SymbolicElement
            if b.output == 0
                continue
            end

            # Otherwise, push the new state
            push!(reduced, b)
        elseif last_item isa SymbolicElement && item isa SymbolicState && item.party == last_item.party
            a = last_item::SymbolicElement
            b = item::SymbolicState
            
            if a.output == 0
                deleteat!(reduced, length(reduced))
            end
            push!(reduced, item)

        else
            # no adjacency rule => just push
            push!(reduced, item)
        end
    end

    return reduced
end
export reduce_monomial

function partial_reduce_monomial(monom::Vector{SymbolicMonomial}, remove::Int64)::Vector{SymbolicMonomial}

    reduced = SymbolicMonomial[]
    
    monom = reorder_monomial_by_party(monom)

    for item in monom
        if isempty(reduced)
            # There's no previous item => just push
            push!(reduced, item)
            continue
        end

        # Compare 'item' with the last item in 'reduced'
        last_item = reduced[end]

        if last_item isa SymbolicElement && item isa SymbolicElement && last_item.party == remove && item.party == remove
            # We have consecutive operators => apply rules (1) and (2)
            a = last_item::SymbolicElement
            b = item::SymbolicElement

            # If projectors => apply othogonality and projectivity 
            if a.kind == 1 && b.kind == 1
                # Rule 1: contradictory outputs => entire monomial => [NONE_ELEMENT]
                if a.input == b.input && a.output != b.output 
                    return SymbolicMonomial[]
                end

                # Rule 2: projective => if same (party, input, output) and both kind=1 => skip 'b'
                if a.input == b.input && a.output == b.output 
                    # P^2 => P => do not push 'b'
                    continue
                end
            end

            # If operators FOR NOW ASSUME QUBIT
            if a.kind == 0 && b.kind == 0
                if a.input == b.input 
                    deleteat!(reduced, length(reduced))
                    continue
                end
            end

            # Otherwise, no conflict => push the new operator
            push!(reduced, b)

        elseif last_item isa SymbolicState && item isa SymbolicState && item.party == remove && last_item.party == remove
            # We have consecutive states => apply rule (3)
            a = last_item::SymbolicState
            b = item::SymbolicState

            # If both kind=1 => skip the second
            if a.kind == 1 && b.kind == 1 && a.number == b.number
                # S^2 => S => do not push 'b'
                continue
            end

            # Otherwise, push the new state
            push!(reduced, b)

        else
            # If last_item is a state and item is an element (or vice versa),
            # no adjacency rule => just push
            push!(reduced, item)
        end
    end

    return reduced
end
export partial_reduce_monomial


function trace_reduce_monomial(monom::Vector{SymbolicMonomial})::Vector{SymbolicMonomial}

    reduced = monom
    n=length(monom)-1

    for _ in 1:n
        monom = vcat(monom[end], monom[1:n])
        monom_ = reduce_monomial(monom)
        if isless_monom(monom_, reduced)
            reduced = monom_
        end
        if monom_ == SymbolicMonomial[]
            break
        end
    end   
    return reduced 

end
export trace_reduce_monomial


function partial_trace_reduce_monomial(monom::Vector{SymbolicMonomial}, remove::Int64)::Vector{SymbolicMonomial}

    reduced = partial_reduce_monomial(monom, remove)
    n=length(monom)-1

    for _ in 1:n

        if monom[end].party == remove
            monom = vcat(monom[end], monom[1:n])
            monom_ = partial_reduce_monomial(monom, remove)
        else
            break
        end

        if isless_monom(monom_, reduced)
            reduced = monom_
        end

        if monom == SymbolicMonomial[]
            break
        end

    end   
    return reduced 

end
export partial_trace_reduce_monomial
