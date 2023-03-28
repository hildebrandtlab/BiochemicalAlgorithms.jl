using DataStructures: Stack, OrderedSet

export reconstruct_fragments!

# Identify two reference atoms.
# Performs a breadth-first search for two additional heavy atoms
# starting from the center atom. These atoms are used as 
# anchor points for attaching the next atom.
function _get_two_reference_atoms(ref_center_atom::DBAtom{T}, allowed::Set{DBAtom{T}}, ref::DBVariant{T}) where {T<:Real}
    # a hash set to remember all those atoms we have already visited
    atom_list = OrderedSet{DBAtom{T}}()
    push!(atom_list, ref_center_atom)

    # abort if we found the three first atoms (beyond the center atom)
    # or we are running out of fresh atoms

    (current, atom_list_rest) = Iterators.peel(atom_list)
    while length(atom_list) < 3 && !isnothing(current)
        for bond in bonds(current, ref)
            next_atom = get_partner(bond, current, ref)

            if (next_atom ∈ allowed) && (next_atom ∉ atom_list)
                push!(atom_list, next_atom)
                if length(atom_list) > 2
                    break
                end
            end
        end
        
        # try the bonds of the next atom in the list
        (current, atom_list_rest) = Iterators.peel(atom_list_rest)
    end

    atom_list = collect(atom_list)

    (
        length(atom_list) == 3, 
        length(atom_list) > 1 ? atom_list[2] : nothing, 
        length(atom_list) > 2 ? atom_list[3] : nothing
    )
end

function reconstruct_fragment_!(f::Fragment{T}, template::DBVariant) where {T<:Real}
    num_inserted_atoms = 0

    # Get a copy of the atom names occurring in the current fragment....
    name_to_atom = Dict(
        a.name => a for a in atoms(f)
    )

    # And add the atoms from the template missing in the reference
    tpl_to_frag = Dict{DBAtom{T}, Atom{T}}()
    transformed = Set{DBAtom{T}}()

    for tpl_atom in template.atoms
        if haskey(name_to_atom, tpl_atom.name)
            # remember that the coordinates of this one are correct
			res_atom = name_to_atom[tpl_atom.name]
            push!(transformed, tpl_atom)
            tpl_to_frag[tpl_atom] = res_atom
        else
            # We create a copy of the existing atom and insert it into
            # the residue. Coordinates are bogus, but we'll correct that
            # later on.
			new_atom = Atom(
                f,
                maximum(atoms_df(parent_system(f)).number)+1, # does this make sense?
                tpl_atom.element,
                tpl_atom.name,
                "",
                tpl_atom.r,
                zero(Vector3{T}),
                zero(Vector3{T}),
                0,
                zero(T),
                zero(T),
                false,
                false
            )

            tpl_to_frag[tpl_atom] = new_atom

            # update the atom count
            num_inserted_atoms += 1
        end
    end

    # We've now made sure that all atoms of the tplate exist in the 
    # reconstructed residue as well (careful, not the other way round!)
    # we can now start to adjust the atom coordinates.

    # If no atoms were in common, there's not much we can do...
    # Trivial solution: no atoms are actually matched to each 
    # other, so we just leave the coordinates the way they
    # are (copy of the tpl coordinates) and return.
    if (!isempty(transformed))
        # Otherwise, we start adjusting coordinates
        # We use a set for BFS
        visited = Set{DBAtom}()
        stack = Stack{DBAtom}()

        push!(stack, first(transformed))

        while (!isempty(stack))
            current = pop!(stack)
            push!(visited, current)

            @debug "Center is $(current.idx) ($(get_full_name(current.name))) visited = " *
                "$(current ∈ visited) transformed = $(current ∈ transformed)" *
                " @ $(current.r)"
                
            @debug "Residue atom is @ $(tpl_to_frag[current].r) (dist = " *
                "$(norm(tpl_to_frag[current].r - current.r)))"
        

            for bond in bonds(current, template)
                next = get_partner(bond, current, template)

                @debug "Examining $(next.idx) ($(get_full_name(next))) visited = " *
                    "$(next ∈ visited) transformed = $(next ∈ transformed)"

                if next ∉ visited
                    push!(stack, next)
                    push!(visited, next)

                    if next ∉ transformed
                        @debug "Searching reference atoms for $(get_full_name(next))"

                        hit, a1, a2 = _get_two_reference_atoms(current, transformed, template)

                        @debug """Reference atoms:  $(isnothing(a1) ? "-" : _get_full_name(a1))""" *
							   """ / $(isnothing(a2) ? "-" : _get_full_name(a2))"""

                        translation, rotation = if hit
                            # we can map all three atoms, great!
                            translation, rotation = match_points(
                                current.r, a1.r, a2.r,
                                tpl_to_frag[current].r, tpl_to_frag[a1].r, tpl_to_frag[a2].r
                            )
                        
                            translation, rotation
                        else
                            # We could map the two center atoms only, which corresponds to 
                            # a simple translation by the difference of the two atom positions.
                            tpl_to_frag[current].r - current.r, T(1)I(3)
                        end

						# Transform the coordinates of the atom we're interest in
                        tpl_to_frag[next].r = rotation * tpl_to_frag[next].r + translation

                        # Remember that we already took care of that guy.
                        push!(transformed, next)

                        @debug "$(get_full_name(next)) is transformed: $(tpl_to_res[next].r) / $(next.r)"
                        @debug "Distance = $(distance(tpl_to_res[next].r, next.r))"
					
                    end                  
                end
            end
        end
    end

    num_inserted_atoms
end

function reconstruct_fragments!(ac::AbstractAtomContainer{T}, fdb::FragmentDB) where {T<:Real}
    num_inserted_atoms = 0

    # iterate over all fragments
    for f in fragments(ac)
            
        # check whether our DB knows the fragment and, if so, retrieve the template
        template = get_reference_fragment(f, fdb)

        if isnothing(template)
            @warn "reconstruct_fragments!(): could not find reference fragment for $(f.name):$(f.number)"

            continue
        end

        num_inserted_atoms += reconstruct_fragment_!(f, template)
    end

    @info "reconstruct_fragments!(): added $(num_inserted_atoms) atoms."
    
    num_inserted_atoms
end