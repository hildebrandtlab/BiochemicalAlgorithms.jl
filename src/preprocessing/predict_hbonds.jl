using CellListMap

export predict_hbonds!, backbone_hydrogen_bonds

"""
    $(TYPEDSIGNATURES)

    Predict hydrogen bonds for a given AtomContainer.

    The `method` parameter selects one of the implemented strategies for H-bond prediction. Bonds can be created between
    the donor and acceptor atoms (e.g., `N` and `O`), which is the default, or between the hydrogen and the acceptor (e.g., `H` and `O`). 
    This behaviour is controlled by the `h_bond_from_donor`-switch.

    Available methods:

        - `:KABSCH_SANDER`: only predicts *protein backbone* hydrogen bonds as required for secondary
          structure prediction with the Kabsch-Sander algorithm `DSSP` (Kabsch W & Sander C (1983). Dictionary of protein secondary 
          structure: pattern recognition of hydrogen-bonded and geometrical features. Biopolymers, 22, 2577-2637".)

        - `:WISHART_ET_AL`: predicts hydrogen bonds between amid and α-hydrogens (`H`/`HA`) and carbonyl oxygens in the backbone (`O`)
          or sidechain oxygens (`OD`, `OE`, `OG`, `OH`). This method follows the criterion given in (Neal, S., Nip, A. M., Zhang, H., 
          and Wishart, D. S. (2003). Rapid and accurate calculation of protein 1H, 13C and 15N chemical shifts. J Biomol NMR, 26(3):215-240.


    Example:

    ```julia
    predict_hbonds!(sys, :WISHART_ET_AL)
    ```
"""
function predict_hbonds!(ac::AbstractAtomContainer{T}, method::Symbol, h_bond_from_donor=true) where {T}
    if (method == :KABSCH_SANDER)
        return predict_hbonds_kabsch_sander!(ac, h_bond_from_donor)
    elseif (method == :WISHART_ET_AL)
        return predict_hbonds_wishart_et_al!(ac, h_bond_from_donor)
    else
        error("predict_hbonds!: Unknown prediction method $(method)!")
    end
end

function predict_hbonds_kabsch_sander!(ac::AbstractAtomContainer{T}, h_bond_from_donor=true) where {T}
    MAX_HBOND_LENGTH = 5.2
    ϵ = -0.5

    # we only iterate over the amino acids in the atom container that are complete
    # in the sense that they contain an "O", an "N", and a "C" atom
    amino_acids = filter(f -> (is_amino_acid(f)  && (["O", "N", "C"] ⊆ atoms_df(f).name)), fragments(ac))

    N_rs = getproperty.(atom_by_name.(amino_acids, "N"), :r)
    O_rs = getproperty.(atom_by_name.(amino_acids, "O"), :r)

    # now, find all pairs of residues with an N-O distance less than 5.2Å
    candidates = neighborlist(N_rs, O_rs, MAX_HBOND_LENGTH)

    candidate_fragments = map(c -> (amino_acids[c[1]], amino_acids[c[2]], c[3]), candidates)
    
    # filter out residues that are adjacent in the chain, and interactions inside the same residue
    # TODO: we should really check for a peptide bond here, as the numbering might be misleading
    filtered_candidates = filter(c -> (
            (c[1] ≠ c[2]) &&
            !is_previous(c[1], c[2]) &&
            !is_next(c[1], c[2])
        ),
        candidate_fragments)

    n_added_hbonds = 0

    # now, for each filtered candidate pair, approximate the hydrogen bond energy
    # TODO: do we need to treat N-terminals differently? right now, we skip them
    for (c1, c2) in filtered_candidates
        if (is_n_terminal(c1))
            continue
        end
        
        N = atom_by_name(c1, "N")
        H = atom_by_name(c1, "H")
        O = atom_by_name(c2, "O")
        C = atom_by_name(c2, "C")

        if isnothing(H)
            continue
        end
       
        energy = 0.42 * 0.20 * 332.0 * (
            1/distance(O, N) + 1/distance(C, H) - 1/distance(O, H) - 1/distance(C, N)
        )

        if energy >= ϵ
            continue
        end

        # does this h-bond already exist?
        existing_hbond_partners = get_partner.(hydrogen_bonds(O), Ref(O))

        if (N ∈ existing_hbond_partners) || (H ∈ existing_hbond_partners)
            continue
        end

        n_added_hbonds += 1

        donor_idx    = h_bond_from_donor ? N.idx : H.idx
        acceptor_idx = O.idx

        Bond(parent_system(N), donor_idx, acceptor_idx, BondOrder.Single, Properties(), Flags([:TYPE__HYDROGEN]))

        @debug "H-Bond found between $(get_full_name(c1)) and $(get_full_name(c2))"
    end

    @info "Added $(n_added_hbonds) hydrogen bonds."
end

"""
    $(TYPEDSIGNATURES)

    Return all backbone hydrogen bonds in the given atom container.

    The criterion for a hydrogen bond to be a backbone bond is as follows:
        - the bond connects a donor nitrogen or hydrogen with the name "N" or "H" with an oxygen named "O"
        - the bond connects atoms in two different fragments

    Example:

    ```julia
    backbone_hydrogen_bonds(sys)
    ```
"""
function backbone_hydrogen_bonds(ac::AbstractAtomContainer{T}) where {T}
    filter(
        b -> begin
            a1, a2 = get_partners(b)
            (
                (parent_fragment(a1) != parent_fragment(a2)) &&
                    (
                        (a1.name ∈ ["N", "H"]) && (a2.name == "O") ||
                        (a2.name ∈ ["N", "H"]) && (a1.name == "O")
                    )
            )
        end,
        hydrogen_bonds(ac)
    )
end