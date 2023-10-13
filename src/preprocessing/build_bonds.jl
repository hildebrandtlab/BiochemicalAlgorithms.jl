using LinearAlgebra

export build_bonds!

function build_fragment_bonds!(
        f::Fragment{T}, 
        connections::Dict{Int, DBConnection},
        fdb::FragmentDB) where {T<:Real}

    # check whether our DB knows the fragment and, if so, retrieve the template
    template = get_reference_fragment(f, fdb)

    if isnothing(template)
        @warn "build_bonds!(): could not find reference fragment for $(f.name)."
        return 0
    end

    @debug "build_bonds!(): building bonds for $(f.name) from template $(template.name)"

    template_names = Dict{AbstractString, DBAtom{T}}(strip(a.name) => a for a in template.atoms)

    num_bonds_built = 0

    # iterate over all fragment atoms
    for atom in atoms(f)
        atom_name = strip(atom.name)

        if !haskey(template_names, atom_name)
            continue
        end

        tpl_atom = template_names[atom_name]

        # we've matched the fragment atom to an atom in the template
        # now, let's check their bonds
        tpl_bonds = filter(b -> b.a1 == tpl_atom.name || b.a2 == tpl_atom.name, template.bonds)

        for tpl_bond in tpl_bonds
            tpl_partner = 
                (tpl_bond.a1 == tpl_atom.name) ? tpl_bond.a2 : ((tpl_bond.a2 == tpl_atom.name) ? tpl_bond.a1 : nothing)

            if isnothing(tpl_partner)
                continue
            end

            # does the fragment contain the partner atom?
            partner_atoms = filter(a -> strip(a.name) == strip(tpl_partner), atoms(f))

            if length(partner_atoms) > 1
                throw("build_bonds!(): corrupt FragmentDB - data encountered. Fragment atoms not unique!")
            end

            if length(partner_atoms) == 1
                partner_atom = partner_atoms[1]

                # ok, find out if the fragment already contains the bond
                bond_index = findfirst(
                    b ->    (b.a1 == atom.idx && b.a2 == partner_atom.idx)
                         || (b.a1 == partner_atom.idx && b.a2 == atom.idx),
                    bonds(f))

                if isnothing(bond_index)
                    Bond(parent_system(f), atom.idx, partner_atom.idx, tpl_bond.order, Properties())
                
                    num_bonds_built += 1
                else
                    bond = bonds(f)[bond_index]

                    bond.order = tpl_bond.order
                end

                break
            end
        end
    end

    # now, store the connections for this fragment so we can connect them together later
    cons = fdb.fragments[f.name].connections

    for con in cons
        # connections have the form:
        # <name> <atom_name> <match_name> <distance> <tolerance>
        # bonds are only built if the atoms are within the tolerance of the expected distance

        # try to find an atom with the name of this connection in the fragment
        atom_idx = findfirst(a -> a.name == con.atom_name, atoms(f))

        if !isnothing(atom_idx)
            atom = atoms(f)[atom_idx].idx

            if haskey(connections, atom)
                @warn "build_bonds!(): duplicate connection encountered!"
            end

            connections[atom] = con
        end
    end

    num_bonds_built
end

function try_build_connection!(a1::Atom, con_1::DBConnection, a2::Atom, con_2::DBConnection)
    # are the connection types compatible?
    if con_1.name != con_2.match_name || con_2.name != con_1.match_name
        return false
    end

    # are the distances compatible?
    distance = norm(a1.r - a2.r)

    if abs(con_1.distance - distance) > con_1.tolerance ||
       abs(con_2.distance - distance) > con_2.tolerance
        return false
    end

    if con_1.order != con_2.order
        @warn "build_bonds!(): inconsistent bond orders"
    end

    # mark disulfide bridges
    props = Properties()
    if a1.name == "SG" && a2.name == "SG"
        props[:DISULPHIDE_BOND] = true
    end

    Bond(parent_system(a1), a1.idx, a2.idx, con_1.order, props)

    return true
end

function build_bonds!(m::AbstractAtomContainer{T}, fdb::FragmentDB) where {T<:Real}
    # while building up individual fragments, we remember inter-fragment connections
    connections = Dict{Int, DBConnection}()

    # is_bound_to is currently quite slow and relies on dynamic dispatch; to speed things
    # up, we precompute the bonds into a dictionary
    bond_cache = Set{Tuple{Int, Int}}()

    for bond in eachbond(m)
        push!(bond_cache, (bond.a1, bond.a2))
        push!(bond_cache, (bond.a2, bond.a1))
    end

    num_bonds_built = 0

    # first, build the bonds inside each fragment
    for frag in fragments(m)
        num_bonds_built += build_fragment_bonds!(frag, connections, fdb)
    end

    # now, connect the fragments that have suitable connections
    num_connections = length(connections)

    if num_connections > 2
        # we test all possible pairings of connection entries
        con_keys = collect(keys(connections))

        for i in 1:num_connections
            for j in i+1:num_connections
                con_i = connections[con_keys[i]]
                con_j = connections[con_keys[j]]

                # only create the bond if it does not exist
                if (con_keys[i], con_keys[j]) ∉ bond_cache
                    if try_build_connection!(atom_by_idx(m, con_keys[i]), con_i, atom_by_idx(m, con_keys[j]), con_j)
                        push!(bond_cache, (con_keys[i], con_keys[j]))
                        push!(bond_cache, (con_keys[j], con_keys[i]))
                        num_bonds_built += 1
                        break
                    end
                end
            end
        end
    end

    @info "build_bonds!(): built $(num_bonds_built) bonds"
end