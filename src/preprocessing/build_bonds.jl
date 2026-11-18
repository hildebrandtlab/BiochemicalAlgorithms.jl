export
    build_bonds!

# Lookup tables for build_fragment_bonds!, built once per build_bonds! call.
# `atoms(f)` and `bonds(f)` each filter the whole system, so calling them for
# every fragment made build_bonds! O(N^2).
#  - `frag_atoms`: fragment idx => idx of its atoms in table order (same rows,
#    same order as `atoms(f)`)
#  - `bond_of`: (atom1 idx, atom2 idx), both orientations => idx of the first
#    bond in table order between the two atoms (what `findfirst` on `bonds(f)`
#    returned)
function _bond_lookup(sys::System{T}) where T
    frag_atoms = Dict{Int, Vector{Int}}()
    at = atoms(sys)
    for (aidx, fidx) in zip(at.idx, at.fragment_idx)
        isnothing(fidx) && continue
        push!(get!(Vector{Int}, frag_atoms, fidx), aidx)
    end
    bond_of = Dict{Tuple{Int, Int}, Int}()
    bt = bonds(sys)
    for (bidx, a1, a2) in zip(bt.idx, bt.atom1_idx, bt.atom2_idx)
        get!(bond_of, (a1, a2), bidx)
        get!(bond_of, (a2, a1), bidx)
    end
    frag_atoms, bond_of
end

function build_fragment_bonds!(
    f::Fragment{T},
    connections::OrderedDict{Int, DBConnection},
    fdb::FragmentDB{T},
    lookup = nothing
) where T
    sys = parent_system(f)
    frag_atoms, bond_of = isnothing(lookup) ? _bond_lookup(sys) : lookup

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
    fatoms = AtomTable{T}(sys, get(frag_atoms, f.idx, Int[]))
    for atom in fatoms
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
            partner_atoms = filter(a -> strip(a.name) == strip(tpl_partner), fatoms)

            if length(partner_atoms) > 1
                @show _show(stdout, fatoms)
                @show _show(stdout, partner_atoms)
                throw("build_bonds!(): corrupt FragmentDB - data encountered. Fragment atoms not unique!")
            end

            if length(partner_atoms) == 1
                partner_atom = partner_atoms[1]

                # ok, find out if the fragment already contains the bond
                bond_idx = get(bond_of, (atom.idx, partner_atom.idx), nothing)

                if isnothing(bond_idx)
                    new_idx = Bond(atom, partner_atom, tpl_bond.order).idx
                    bond_of[(atom.idx, partner_atom.idx)] = new_idx
                    bond_of[(partner_atom.idx, atom.idx)] = new_idx

                    num_bonds_built += 1
                else
                    bond = bond_by_idx(sys, bond_idx)

                    bond.order = tpl_bond.order
                end

                continue
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
        atom_idx = findfirst(a -> a.name == con.atom_name, fatoms)

        if !isnothing(atom_idx)
            atom = fatoms[atom_idx].idx

            if haskey(connections, atom)
                @warn "build_bonds!(): duplicate connection encountered!"
            end

            connections[atom] = con
        end
    end

    num_bonds_built
end

function try_build_connection!(a1::Atom{T}, con1::DBConnection{T}, a2::Atom{T}, con2::DBConnection{T}) where T
    # are the connection types compatible?
    if con1.name != con2.match_name || con2.name != con1.match_name
        return false
    end

    # are the distances compatible?
    distance = norm(a1.r .- a2.r)

    if abs(con1.distance - distance) > con1.tolerance ||
       abs(con2.distance - distance) > con2.tolerance
        return false
    end

    if con1.order != con2.order
        @warn "build_bonds!(): inconsistent bond orders"
    end

    # mark disulphide bridges
    flags = Flags()
    if a1.name == "SG" && a2.name == "SG"
        push!(flags, :TYPE__DISULPHIDE_BOND)
    end

    Bond(a1, a2, con1.order; flags = flags)

    return true
end

"""
    build_bonds!(::AbstractAtomContainer{Float32})
    build_bonds!(::AbstractAtomContainer{T}, ::FragmentDB{T})

Attempts to construct missing bonds in the given container, according to the
default/given fragment database, and returns the number of built bonds.

!!! note
    This preprocessor expects fragment and atom names to be normalized and
    fragment terminals to be labeled (cf. [`normalize_names!`](@ref) and
    [`label_terminal_fragments!`](@ref)).
"""
function build_bonds!(m::AbstractAtomContainer{T}, fdb::FragmentDB{T}) where T
    # while building up individual fragments, we remember inter-fragment connections
    connections = OrderedDict{Int, DBConnection}()

    # is_bound_to is currently quite slow and relies on dynamic dispatch; to speed things
    # up, we precompute the bonds into a dictionary
    bond_cache = Set{Tuple{Int, Int}}()

    for bond in bonds(m)
        push!(bond_cache, (bond.atom1_idx, bond.atom2_idx))
        push!(bond_cache, (bond.atom2_idx, bond.atom1_idx))
    end

    num_bonds_built = 0

    # first, build the bonds inside each fragment
    lookup = _bond_lookup(parent_system(m))
    for frag in fragments(m)
        num_bonds_built += build_fragment_bonds!(frag, connections, fdb, lookup)
    end

    # now, connect the fragments that have suitable connections
    num_connections = length(connections)

    if num_connections > 2
        # we test all possible pairings of connection entries
        con_keys = collect(keys(connections))
        con_vals = collect(values(connections))
        # look the atoms up once, not once per pair
        con_atoms = [atom_by_idx(m, k) for k in con_keys]

        for i in 1:num_connections
            for j in i+1:num_connections
                con_i = con_vals[i]
                con_j = con_vals[j]

                # only create the bond if it does not exist
                if (con_keys[i], con_keys[j]) ∉ bond_cache
                    if try_build_connection!(con_atoms[i], con_i, con_atoms[j], con_j)
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

    num_bonds_built
end

@inline function build_bonds!(m::AbstractAtomContainer{Float32})
    build_bonds!(m, default_fragmentdb())
end
