using LinearAlgebra
using Mendeleev
using Rotations

export add_hydrogens!

function place_peptide_bond_h_!(frag::Fragment{T}) where T
    BOND_LENGTH_N_H::T = 1.020

    if is_amino_acid(frag) && !is_n_terminal(frag)
        prev_residue = get_previous(frag)

        if isnothing(prev_residue) # should not happen, as we are not the n-terminal
            @warn "Inconsistent chain encountered!"

            return false
        end

        n_atom = atom_by_name(frag, "N")
        o_atom = atom_by_name(frag, "O")
        c_atom = atom_by_name(frag, "C")

        if isnothing(n_atom) || isnothing(o_atom) || isnothing(c_atom) || nbonds(n_atom) >= 3
            return false
        end

        # place the hydrogen according to the planarity of the peptide bond
        OC = o_atom.r - c_atom.r
        length = norm(OC)

        if length ≤ eps(T)
            return false
        end

        # add the new hydrogen
        h_atom = Atom(
            frag,
            maximum(atoms(parent_system(frag)).number)+1, # does this make sense?
            Elements.H,
            name = "H",
            r = n_atom.r .- (OC .* BOND_LENGTH_N_H) / length,
        )

        Bond(parent_system(frag), n_atom.idx, h_atom.idx, BondOrder.Single)

        return true
    end

    return false
end

function _get_connectivity(atom)
    try
       element = atom.element

        group = chem_elements[Symbol(element)].group.no

        if group < 1 || (group > 2 && group < 13)
            return 0
        end

        electrons = 0
        if group < 3
            electrons = group - atom.formal_charge
        else
            electrons = 18 - group + atom.formal_charge
        end

        if electrons < 0
            @error "Could not calculate number of electrons for $(parent_fragment(atom).name):$(atom.name)!"
        end

        return electrons
    catch
        return 0
    end
end

function _count_bond_orders(a::Atom{T}) where {T<:Real}
    nr = zero(T)

    for b in bonds(a)
        if b.order == BondOrder.Aromatic
            nr += T(1.5)
        elseif b.order >= BondOrder.Single && b.order <= BondOrder.Quadruple
            nr += Int(b.order)
        end
    end

    Int(floor(nr))
end

# Calculate the reference bond length value using a modified Shoemaker-Stevenson rule
# (taken from MMFF94 force field)
function _get_bond_length(element::ElementType, mmff_params::MMFF94Parameters{T}) where {T<:Real}
    # currently only supports atoms up to Xenon
    if (element >= Elements.Xe)
        return one(T)
    end

    radius = mmff_params.radii[Int(element)]

    # if no radius available for the element:
    if (radius == 0)
        return one(T);
    end

    radius_h = mmff_params.radii[Int(Elements.H)]

    #  c and n are constants defined in R.Blom and A. Haaland,
    #  J. Molec. Struc, 1985, 128, 21-27.
    # calculate proportionality constant c
    c = T(0.05)

    # POWER
	n = T(1.4)

    diff_e = abs(
        mmff_params.electronegativities[Int(Elements.H)] -
        mmff_params.electronegativities[Int(element)])

    # FORMULA
	radius + radius_h - c * diff_e^n
end


function _add_hydrogen!(atom::Atom{T}, r::Vector3{T}, atom_nr::Int) where {T<:Real}
    # add the new hydrogen
    h_atom = Atom(
        parent_fragment(atom),
        maximum(atoms(parent_system(atom)).number)+1, # does this make sense?
        Elements.H,
        name = "H" * ((atom_nr > 1) ? (string(atom_nr) * atom.name) : atom.name),
        r = r
    )

    Bond(parent_system(atom), atom.idx, h_atom.idx, BondOrder.Single)
end

function _get_normal(v::Vector3{T}) where {T<:Real}
    n = normalize(cross(v, Vector3{T}(1, 0, 0)))

    if isnan(n[1])
        n = normalize(cross(v, Vector3{T}(0, 1, 0)))

        if isnan(n[1])
            n = normalize(cross(v, Vector3{T}(0, 0, 1)))
        end
    end

    n
end

function _handle_atom!(atom::Atom{T}, mmff_params::MMFF94Parameters{T}, is_ring_atom::Bool, atom_nr::Int) where {T<:Real}
    # prevent adding hydrogens to, e.g., aromatic Carboxy group
    if nbonds(atom) == 1 && first(bonds(atom)).order == BondOrder.Aromatic
        return 0, atom_nr
    end

    # determine the number of electrons that have to be delivered through bonds
    con = _get_connectivity(atom)

    sum_bond_orders = _count_bond_orders(atom)

    h_to_add = con - sum_bond_orders

    if h_to_add <= 0
        return 0, atom_nr
    end

    bond_length = _get_bond_length(atom.element, mmff_params)

    nr_bonds = nbonds(atom)

    # one bond and one Hydrogen missing: (e.g. H-F)
    if con == 1
        p = atom.r - Vector3{T}(bond_length, 0, 0)

        _add_hydrogen!(atom, p, atom_nr)

        return 1, atom_nr+1
    end

    # linear compounds
    if (h_to_add == 1 && nr_bonds == 1 && sum_bond_orders > 2)
        partner_atom = get_partner(first(bonds(atom)), atom)

        diff = normalize(partner_atom.r - atom.r)

        if isnan(diff[1])
            diff = Vector3{T}(0, 1, 0)
        end

        diff *= bond_length;

        _add_hydrogen!(atom, atom_position - diff, atom_nr)

        return 1, atom_nr+1
    end

    # two partner atoms and a planar 106 degree angle: (e.g. H-O-H)
    if (con == 2)
        if (h_to_add == 2)
            # add first bond
            p = atom.r - Vector3{T}(bond_length, 0, 0)

            _add_hydrogen!(atom, p, atom_nr)

            # add second bond
            added_h, atom_nr = _handle_atom!(atom, mmff_params, is_ring_atom, atom_nr+1)

            return 1+added_h, atom_nr
        else
            # h_to_add == 1
            bv = atom.r - get_partner(first(bonds(atom)), atom).r
            axis = _get_normal(bv)

            rotation = AngleAxis{T}(deg2rad(106), axis...)
            bv = normalize(rotation * bv)

            if isnan(bv[1])
                bv = Vector3{T}(0, 0, 1)
            end

            bv *= bond_length

            _add_hydrogen!(atom, atom.r - bv, atom_nr)

            return 1, atom_nr+1
        end
    end

    # Ring atoms
    if is_ring_atom
        bs = bonds(atom)

        v1 = normalize(get_partner(bs[1], atom).r - atom.r)
        v2 = normalize(get_partner(bs[2], atom).r - atom.r)

        if isnan(v1[1])
            v1 = Vector3{T}(1, 0, 0)
        end

        if isnan(v2[1])
            v2 = Vector3{T}(0, 1, 0)
        end

        v3 = normalize(-(v1 + v2))

        if isnan(v3[1])
            v3 = Vector3{T}(0, 0, 1)
        end

        v3 *= bond_length

        # e.g. Nitrogen in Ring
        if (con == 3 || (con == 4 && atom.formal_charge == 1))
            if (h_to_add == 1)
                _add_hydrogen!(atom, atom.r + v3, atom_nr);

                return 1, atom_nr+1
            end
        end

        # Carbon in Ring
        if (con == 4)
            vx = v2 - v1

            if (h_to_add == 2)
                rotation = AngleAxis{T}(deg2rad(60), vx...)
                _add_hydrogen!(atom, atom.r + rotation * v3, atom_nr)

                atom_nr += 1

                rotation = AngleAxis{T}(deg2rad(-60), vx...)
                _add_hydrogen!(atom, atom.r + rotation * v3, atom_nr)

                return 2, atom_nr+1
            end

            if (h_to_add == 1)
                if nbonds(atom) == 3
                    # maybe another Hydrogen was already added?
                    for bond in bonds(atom)
                        partner = get_partner(bond, atom)

                        if partner.element != Elements.H
                            continue
                        end

                        rotation = AngleAxis{T}(deg2rad(60), vx...)
                        partner.r = atom_position + rotation * v3

                        rotation = AngleAxis{T}(deg2rad(-60), vx...)
                        _add_hydrogen!(atom, atom.r + rotation * v3, atom_nr)

                        return 1, atom_nr+1
                    end

					# planar and 1 atom to add:
					if nbonds(atom) == 2
                        _add_hydrogen!(atom, atom.r + v3, atom_nr)

                        return 1, atom_nr+1
                    else
                        v3 = normalize(-cross(v1, v2))

                        if isnan(v3[1])
                            v3 = Vector3{T}(0, 0, 1)
                        end

						v3 *= bond_length

                        _add_hydrogen!(atom, atom.position + v3, atom_nr)
                        return 1, atom_nr+1
                    end
                end
            end
        end

        # does the atom have at least one bond that is not a single bond?
        if any(map(b -> b.order > BondOrder.Single, bonds(a)))
            bv = normalize(get_partner(first(bonds(atom)), atom).r - atom.r)

			# e.g. (C[-H][-H]=O) or (H-N=O)
			if ((con == 4 && h_to_add == 2) ||
					(con == 3 && h_to_add == 1))

                if isnan(bv[1])
                    bv = Vector3{T}(-1, 0, 0)
                end

                axis = _get_normal(bv)

                rotation = AngleAxis{T}(deg2rad(120), axis...)

				bv = rotation * b
				bv *= bond_length

                _add_hydrogen!(atom, atom.r + bv, atom_nr)

                # add second bond ?
                if (h_to_add == 2)
                    added_h, atom_nr = _handle_atom!(atom, mmff_params, is_ring_atom, atom_nr+1)

                    return 1+added_h, atom_nr
                else
                    return 1, atom_nr+1
                end
            end

            # e.g. (C[-H][-H]=O)
			if (con == 4 && h_to_add == 1)
                bs = bonds(atom)

                p1 = normalize(get_partner(bs[1], atom).r - atom.r)
                p2 = normalize(get_partner(bs[2], atom).r - atom.r)

                if isnan(p1[1])
                    p1 = Vector3{T}(0, 1, 0)
                end

                if isnan(p2[1])
                    p2 = Vector3{T}(0, 0, 1)
                end

                v = normalize(p1 + p2)

                if isnan(v[1])
                    v = Vector3{T}(1, 0, 0)
                end

				v *= bond_length

                _add_hydrogen!(atom, atom.r - v, atom_nr)

                return 1, atom_nr+1
            end
        end

        # three partner atoms and a 106 degree angle: (NH3)
		if (con == 3)
			if (h_to_add == 3)
				# add first bond
                p = atom.r - Vector3{T}(bond_length, 0, 0)

                _add_hydrogen!(atom, p, atom_nr)


                # and continue recursively
                added_h, atom_nr = _handle_atom!(atom, mmff_params, is_ring_atom, atom_nr+1)

                return 1+added_h, atom_nr
            end

			if (h_to_add == 2)
				# add second bond
                bv = normalize(get_partner(first(bonds(atom)), atom).r - atom.r)

                if isnan(bv[1])
                    bv = Vector3{T}(0, 1, 0)
                end

                axis = _get_normal(bv)

                rotation = AngleAxis{T}(deg2rad(112.754181), axis...)
				axis = rotation * bv

                rotation = AngleAxis{T}(deg2rad(120), axis...)
                new_pos = rotation * bv

                _add_hydrogen!(atom, atom.r + new_pos * bond_length, atom_nr)
                _add_hydrogen!(atom, atom.r + (rotation * new_pos) * bond_length, atom_nr+1)

                return 2, atom_nr+2
            end

			if (h_to_add == 1)
				#TODO: This can be improved further. However the approximation
				#      should provide a good placement.
                bs = bonds(atom)

                p1 = normalize(get_partner(bs[1], atom).r)
                p2 = normalize(get_partner(bs[2], atom).r)

                # connection line between the two partner atoms:
				d = p2 - p1

                if norm(d) ≤ eps(T)
                    _add_hydrogen!(atom, atom.r - Vector3{T}(0, 1, 0), atom_nr)

                    return 1, atom_nr+1
                end

				# Point between two partner aoms:
				p = p1 + d / T(2.0)
                d2 = p - atom.r

                rotation = AngleAxis{T}(deg2rad(117.3), d...)
                v = normalize(rotation * d2)

                if isnan(v[1])
                    v = Vector3{T}(0, 0, 1)
                end

                v *= bond_length

                _add_hydrogen!(atom, atom.r + v, atom_nr)

                return 1, atom_nr+1
            end
        end

        # Carbon without double bonds and not in ring:
		# tetrahedral: e.g. CH4
		if (con == 4)
            bs = bonds(atom)

			if (h_to_add == 4)
				# add first hydrogen randomly
                _add_hydrogen!(atom, atom.r + Vector3{T}(bond_length, 0, 0), atom_nr)

                # continue with the next case:
                added_h, atom_nr = _handle_atom!(atom, mmff_params, is_ring_atom, atom_nr+1)

                return 1+added_h, atom_nr
            end

            v = normalize(get_partner(bs[1], atom).r - atom.r)

            if isnan(v[1])
                v = Vector3{T}(0, 1, 0)
            end

			if (h_to_add == 3)
				# Rotate the partner atom around the target atom
				# in order to obtain the first hydrogen
				axis = _get_normal(v)

                rotation = AngleAxis{T}(deg2rad(109.471221), axis...)
                new_pos = (rotation * v) * bond_length

                _add_hydrogen!(atom, atom.r + new_pos, atom_nr)

				# Create two copies of the first hydrogen by rotating
				# for 120 degrees.
                rotation = AngleAxis(deg2rad(120), v...)
				new_pos = rotation * new_pos

                _add_hydrogen!(atom, atom.r + new_pos, atom_nr+1)

                new_pos = rotation * new_pos

                _add_hydrogen!(atom, atom.r + new_pos, atom_nr+2)

                return 3, atom_nr+3
            end

            v2 = normalize(get_partner(bs[2], atom).r - atom.r)

            if isnan(v2[1])
                v2 = Vector3{T}(0, 0, 1)
            end

			if (h_to_add == 2)
				# Create a normal to the plane defined by the atom and its two partners
                v12 = normalize(get_partner(bs[2], atom).r - get_partner(bs[1], atom).r)

                if isnan(v12[1])
                    v12 = Vector3{T}(0, 1, 0)
                end

                normal = cross(v, v2)

				# The new hydrogen atoms are obtained by rotating the normal around the
				# connection between the two partner atoms
                rotation = AngleAxis{T}(deg2rad(-(180 - 109.471221)/2.0), v12...)
                new_pos = (rotation * normal) * bond_length

                _add_hydrogen(atom, atom.r + new_pos, atom_nr)

                rotation = AngleAxis{T}(deg2rad(-109.471221), v12...)
                new_pos = (rotation * new_pos) * bond_length

                _add_hydrogen!(atom, atom.r + new_pos, atom_nr+1)

                return 2, atom_nr+2
            end

			if (h_to_add == 1)
                v3 = normalize(get_partner(bs[3], atom).r - atom.r)

                if isnan(v3[1])
                    v3 = Vector3{T}(1, 0, 0)
                end

                v4 = normalize(v + v2 + v3)

                if isnan(v4[1])
                    v4 = Vector3{T}(1, 0, 0)
                end

				v4 *= bond_length

                _add_hydrogen!(atom, atom.r - v4, atom_nr)

                return 1, atom_nr+1
			end
		end
    end

    return 0, atom_nr
end

function add_hydrogens!(ac::AbstractAtomContainer{T}) where {T<:Real}
    # we will need those parameters later, and it makes sense to only parse the files once
    mmff_params = MMFF94Parameters{T}()

    # the same holds for the set of ring atoms
    ring_atoms = is_ring_atom(ac)

    # we need to collect all atoms here, before any hydrogen has been added
    all_atoms = atoms(ac)

    # first, place all peptide bond hydrogens
    num_added_hydrogens = 0
    for frag in fragments(ac)
        if is_amino_acid(frag)
            if place_peptide_bond_h_!(frag)
                num_added_hydrogens += 1
            end
        end
    end

    # then, iterate over all atoms
    for (i,atom) in enumerate(all_atoms)
        atom_nr = 1

        added_h, atom_nr = _handle_atom!(atom, mmff_params, ring_atoms[i], atom_nr)

        num_added_hydrogens += added_h
    end

    @info "Added $(num_added_hydrogens) atoms."

    num_added_hydrogens
end