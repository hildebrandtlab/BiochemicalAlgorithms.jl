export
    load_pubchem_json

@enumx _PCBondAnnotation begin
    BA_CROSSED = 1     # Double Bond that can be both Cis/Trans
    BA_DASHED          # Hydrogen-Bond (3D Only?)
    BA_WAVY            # Unknown Stereochemistry
    BA_DOTTED          # Complex/Fractional
    BA_WEDGE_UP        # Above-Plane
    BA_WEDGE_DOWN      # Below-Plane
    BA_ARROW           # Dative
    BA_AROMATIC        # Aromatic
    BA_RESONANCE       # Resonance
    BA_BOLD            # Fat Bond (Non-Specific User Interpreted Information)
    BA_FISCHER         # Interpret Bond Stereo using Fischer Conventions
    BA_CLOSECONTACT    # Identification of Atom-Atom Close Contacts (3D Only)
    BA_UNKNOWN = 255   # Unspecified or Unknown Atom-Atom Annotation
end

const _CT_3D = 2
const _CT_UNITS_NANOMETERS = 11

function _convert_coordinates(pb_coords_vec::JSON3.Array, ::Type{T}) where {T <: Real}
    if isempty(pb_coords_vec)
        return []
    end

    result = []

    for pb_coords in pb_coords_vec
        # first, figure out what unit to use for the Coordinates
        # we interpret pixels, points, and stdbonds as Angstroms
        scale = _CT_UNITS_NANOMETERS in pb_coords.type ? 10 : 1

        # do we have 3d information, or only 2d?
        is_3d = _CT_3D in pb_coords.type

        for c in pb_coords.conformers
            converted = Array{Vector3{T}}(undef, length(c.x))

            for i in 1:length(c.x)
                converted[i] = Vector3{T}(c.x[i], c.y[i], (is_3d && haskey(c, :z)) ? c.z[i] : 0) * scale
            end

            push!(result, converted)
        end
    end

    result
end

function _parse_atoms!(mol::Molecule, compound::JSON3.Object, ::Type{T} = Float32) where {T <: Real}
    if haskey(compound, :atoms) && haskey(compound, :coords)
        conformers = _convert_coordinates(compound.coords, T)

        for i in eachindex(compound.atoms.aid)
            for j in eachindex(conformers)
                Atom(mol,
                    compound.atoms.aid[i],
                    haskey(compound.atoms, :element)
                        ? ElementType(Int(compound.atoms.element[i]))
                        : Elements.Unknown;
                    atom_type = haskey(compound.atoms, :label)
                        ? compound.atoms.label[i].value # does the label contain the atom type?
                        : "",
                    r = T.(conformers[j][i]),
                    formal_charge = haskey(compound.atoms, :charge)
                        ? Int(compound.atoms.charge[i])
                        : 0,
                    frame_id = j
                )
            end
        end
    end
end

function _parse_bonds!(mol::Molecule, compound::JSON3.Object, ::Type{T} = Float32) where {T <: Real}
    aidx = Dict(a.number => a.idx for a in atoms(parent_system(mol)))
    if haskey(compound, :bonds)
        for i in eachindex(compound.bonds.aid1)
            aid1 = compound.bonds.aid1[i]
            aid2 = compound.bonds.aid2[i]
            order = Int(compound.bonds.order[i])
            properties = Properties()

            # check for style PCDrawAnnotations
            annotations = []
            for j in eachindex(compound.coords)

                for k in eachindex(compound.coords[j].conformers)
                    conformer = compound.coords[j].conformers[k]
                    if haskey(conformer, :style)
                        if conformer.style.aid1[i] == aid1 && conformer.style.aid2[i] == aid2
                           push!(annotations, string(_PCBondAnnotation.T(conformer.style.annotation[i])))
                        end
                    end
                end
            end
            if length(annotations) != 0
                properties[:PCBondAnnotation_for_conformer] = annotations
            end

            Bond(mol,
                aidx[aid1],
                aidx[aid2],
                (order <= 4) ? BondOrderType(order) : BondOrder.Unknown;
                properties = properties
            )
        end
    end
end

# TODO:
# - modelling of properties
# - currently: urn.label as key and as value PCInfoDataValue
function _parse_props!(mol::Molecule, compound::JSON3.Object)
    if haskey(compound, :props)
        for i in eachindex(compound.props)
            d = compound.props[i]
            mol.properties[Symbol(d.urn.label)] = d
        end
    end
end

"""
    load_pubchem_json(fname::String, ::Type{T} = Float32) -> System{T}

Read a PubChem JSON file.

!!! note
    Conformers are stored as frames.
"""
function load_pubchem_json(fname::String, ::Type{T} = Float32) where {T <: Real}
    # TODO: full conversion, adding all properties

    pb = JSON3.read(read(fname, String))
    sys = System{T}()
    for compound in pb.PC_Compounds
        # for now, use the file name as the name for the molecule
        mol = Molecule(sys; name = basename(fname) * "_" * string(compound.id.id.cid))
        _parse_atoms!(mol, compound, T)
        _parse_bonds!(mol, compound, T)
        _parse_props!(mol, compound)
    end
    sys
end
