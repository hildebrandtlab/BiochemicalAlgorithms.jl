export load_pubchem_json

using BiochemicalAlgorithms: Molecule, Atom, BondOrder, BondOrderType, Bond, Elements, ElementType, Vector3

using StructTypes
using JSON3

mutable struct PCAtomInt
    aid::Int64
    value::Int64

    PCAtomInt() = new()
end
StructTypes.StructType(::Type{PCAtomInt}) = StructTypes.Mutable()

mutable struct PCAtomString
    aid::Int64
    value::String

    PCAtomString() = new()
end
StructTypes.StructType(::Type{PCAtomString}) = StructTypes.Mutable()

@enum(PCUrnDataType,
    STRING = 1,    # String                             [maps to a VisibleString]
    STRINGLIST,    # List of Strings                    [maps to VisibleString list]
    INT,           # 32-Bit Signed Integer              [maps to an INTEGER]
    INTVEC,        # Vector of 32-Bit Signed Integer    [maps to INTEGER vector]
    UINT,          # 32-Bit Unsigned Integer            [maps to an INTEGER]
    UINTVEC,       # Vector of 32-Bit Unsigned Integer  [maps to INTEGER vector]
    DOUBLE,        # 64-Bit Float                       [maps to a REAL]
    DOUBLEVEC,     # Vector of Double                   [maps to REAL vector]
    BOOL,          # Boolean or Binary value            [maps to a BOOLEAN]
    BOOLVEC,       # Boolean Vector                     [maps to BOOLEAN vector]
    UINT64,        # 64-Bit Unsigned Integer (Hex form) [maps to a VisibleString]
    BINARY,        # Binary Data Blob                   [maps to an OCTET STRING]
    URL,           # URL                                [maps to a VisibleString]
    UNICODE,       # UniCode String                     [maps to a VisibleString]
    DATE,          # ISO8601 Date                       [maps to a Date]
    FINGERPRINT,   # Binary Fingerprint (Gzip'ped bit   [maps to an OCTET STRING]
                   #   list w/ 4-Byte prefix denoting bit list length)
    UNKNOWN = 255    # Unknown Data Type               [maps to a set of VisibleString]
)
StructTypes.StructType(::Type{PCUrnDataType}) = StructTypes.NumberType()
StructTypes.numbertype(::Type{PCUrnDataType}) = UInt8

mutable struct PCUrn
    label::String                          # Generic Name or Label for Display  [e.g., "Log P"]
    name::Union{Nothing,String}            # Qualified Name  [e.g., "XlogP"]
    datatype::Union{Nothing,PCUrnDataType} # Specific Data Type of Value  [e.g., binary]
    parameters::Union{Nothing,String}      # Implementation Parameter  [e.g., "metal=0"]
    implementation::Union{Nothing,String}  # Implementation Name  [e.g., "E_XlogP"]
    version::Union{Nothing,String}         # Implementation Version  [e.g., "3.317"]
    software::Union{Nothing,String}        # Implementation Software  [e.g., "Cactvs"]
    source::Union{Nothing,String}          # Implementation Organization  [e.g., "xemistry.com"]
    release::Union{Nothing,String}         # NCBI Implementation Release  [e.g., "10.25.2005"]

    PCUrn() = new()
end
StructTypes.StructType(::Type{PCUrn}) = StructTypes.Mutable()

mutable struct PCDateStd
    year::UInt32                   # full year (including 1900)
    month::Union{Nothing,UInt8}    # month (1-12)
    day::Union{Nothing,UInt8}      # day of month (1-31)
    season::Union{Nothing,String}  # for "spring", "may-june", etc
    hour::Union{Nothing,UInt8}     # hour of day (0-23)
    minute::Union{Nothing,UInt8}   # minute of hour (0-59)
    second::Union{Nothing,UInt8}   # second of minute (0-59)

    PCDateStd() = new()
end
StructTypes.StructType(::Type{PCDateStd}) = StructTypes.Mutable()

mutable struct PCDate
    str::Union{Nothing,String}    # for those unparsed dates
    std::Union{Nothing,PCDateStd} # use this if you can

    PCDate() = new()
end
StructTypes.StructType(::Type{PCDate}) = StructTypes.Mutable()

mutable struct PCInfoDataValue
    bval::Bool            # Boolean or Binary
    bvec::Vector{Bool}    # Boolean Vector
    ival::Int64           # Integer (signed or unsigned)
    ivec::Vector{Int64}   # Integer Vector
    fval::Float64         # Float or Double
    fvec::Vector{Float64} # Double Vector
    sval::String          # String
    slist::Vector{String} # List of Strings
    date::PCDate          # Date
    binary::Vector{UInt8} # Binary Data
    bitlist::Vector{Bool} # Bit List (specialized version of Boolean vector)

    PCInfoDataValue() = new()
end
StructTypes.StructType(::Type{PCInfoDataValue}) = StructTypes.Mutable()

mutable struct PCInfoData
    urn::PCUrn #   Universal Resource Name  [for Value Qualification]
    value::PCInfoDataValue
  
    PCInfoData() = new()
end
StructTypes.StructType(::Type{PCInfoData}) = StructTypes.Mutable()

mutable struct PCCount
    heavy_atom::UInt32        # Total count of non-Hydrogen (Heavy) Atoms

    # StereoChemistry Counts
    atom_chiral::UInt32       # Total count of (SP3) Chiral Atoms
    atom_chiral_def::UInt32   # Total count of Defined (SP3) Chiral Atoms
    atom_chiral_undef::UInt32 # Total count of Undefined (SP3) Chiral Atoms
    bond_chiral::UInt32       # Total count of (SP2) Chiral Bonds
    bond_chiral_def::UInt32   # Total count of (SP2) Defined Chiral Bonds
    bond_chiral_undef::UInt32 # Total count of (SP2) Undefined Chiral Bonds

    # Isotopic Counts
    isotope_atom::UInt32      # Total count of Atoms with Isotopic Information

   # Discrete Structure Counts
    covalent_unit::UInt32     # Total count of covalently_bonded units in the record
    tautomers::UInt32         # Number of possible tautomers (Max. 999)

    PCCount() = new()
end
StructTypes.StructType(::Type{PCCount}) = StructTypes.Mutable()

# we need to map the json names, which contain -, to _-using symbols
count_fieldnames = []
StructTypes.foreachfield((i, name, FT, v) -> push!(count_fieldnames, (name, Symbol(replace(String(name), "_" => "-")))), PCCount())
StructTypes.names(::Type{PCCount}) = Tuple(count_fieldnames)

@enum(PCBondAnnotation,
    BA_CROSSED = 1,    # Double Bond that can be both Cis/Trans
    BA_DASHED,       # Hydrogen-Bond (3D Only?)
    BA_WAVY,         # Unknown Stereochemistry
    BA_DOTTED,       # Complex/Fractional
    BA_WEDGE_UP,     # Above-Plane
    BA_WEDGE_DOWN,   # Below-Plane
    BA_ARROW,        # Dative
    BA_AROMATIC,     # Aromatic
    BA_RESONANCE,    # Resonance
    BA_BOLD,         # Fat Bond (Non-Specific User Interpreted Information)
    BA_FISCHER,      # Interpret Bond Stereo using Fischer Conventions
    BA_CLOSECONTACT, # Identification of Atom-Atom Close Contacts (3D Only)
    BA_UNKNOWN = 255   # Unspecified or Unknown Atom-Atom Annotation
)

StructTypes.StructType(::Type{PCBondAnnotation}) = StructTypes.NumberType()
StructTypes.numbertype(::Type{PCBondAnnotation}) = UInt8

mutable struct PCDrawAnnotations
    annotation::Vector{PCBondAnnotation}
    aid1::Vector{Int64}
    aid2::Vector{Int64}

    PCDrawAnnotations() = new()
end
StructTypes.StructType(::Type{PCDrawAnnotations}) = StructTypes.Mutable()

mutable struct PCConformer
    x::Vector{Float32}
    y::Vector{Float32}
    z::Union{Nothing,Vector{Float32}}

    style::Union{Nothing,PCDrawAnnotations}
    data::Union{Nothing,Vector{PCInfoData}} # Data Associated with this Conformer

    PCConformer() = new()
end
StructTypes.StructType(::Type{PCConformer}) = StructTypes.Mutable()

@enum(PCCoordinateType,
        CT_2D = 1,           # 2D Coordinates
        CT_3D,               # 3D Coordinates (should also indicate units, below)
        CT_SUBMITTED,        # Depositor Provided Coordinates
        CT_EXPERIMENTAL,     # Experimentally Determined Coordinates
        CT_COMPUTED,         # Computed Coordinates
        CT_STANDARDIZED,     # Standardized Coordinates
        CT_AUGMENTED,        # Hybrid Original with Computed Coordinates (e.g., explicit H)
        CT_ALIGNED,          # Template used to align drawing
        CT_COMPACT,          # Drawing uses shorthand forms (e.g., COOH, OCH3, Et, etc.)
        CT_UNITS_ANGSTROMS,  # (3D) Coordinate units are Angstroms
        CT_UNITS_NANOMETERS, # (3D) Coordinate units are nanometers
        CT_UNITS_PIXEL,      # (2D) Coordinate units are pixels
        CT_UNITS_POINTS,     # (2D) Coordinate units are points
        CT_UNITS_STDBONDS,   # (2D) Coordinate units are standard bond lengths (1.0)
        CT_UNITS_UNKNOWN = 255 # Coordinate units are unknown or unspecified
)
StructTypes.StructType(::Type{PCCoordinateType}) = StructTypes.NumberType()
StructTypes.numbertype(::Type{PCCoordinateType}) = UInt8

mutable struct PCCoordinates
    type::Vector{PCCoordinateType} # Coordinate Type Information (vector)
    aid::Vector{Int64}             # Conformer Atom IDs (vector) -- to be kept synchronized with Conformers
  
    conformers::Union{Nothing,Vector{PCConformer}}  # Conformers for this Coordinate Set
    atomlabels::Union{Nothing,Vector{PCAtomString}} # Atom labels for Conformer Set
    data::Union{Nothing,Vector{PCInfoData}}         # Data Associated with these Coordinates

    PCCoordinates() = new()
end
StructTypes.StructType(::Type{PCCoordinates}) = StructTypes.Mutable()

@enum(PCStereoGroupType,
  SGT_ABSOLUTE = 1,    # Absolute configuration is known
  SGT_OR,            # Relative configuration is known (absolute configuration is unknown)
  SGT_AND,           # Mixture of stereoisomers
  SGT_UNKNOWN = 255    # Unknown configuration type
)
StructTypes.StructType(::Type{PCStereoGroupType}) = StructTypes.NumberType()
StructTypes.numbertype(::Type{PCStereoGroupType}) = UInt8

mutable struct PCStereoGroup
    type::PCStereoGroupType 
    aid::Vector{Int64}       # Atom Identifiers of atoms in this group

    PCStereoGroup() = new()
end
StructTypes.StructType(::Type{PCStereoGroup}) = StructTypes.Mutable()

mutable struct PCStereoPentagonalBiPyramid
    center::Int64   # Atom ID of Atom Center
    top::Int64      # Atom ID of Atom In-Plane and at the Top
    bottom::Int64   # Atom ID of Atom In-Plane and at the Bottom
    left::Int64     # Atom ID of Atom In-Plane and at the Left
    labove::Int64   # Atom ID of Atom Above the Plane on the Left
    lbelow::Int64   # Atom ID of Atom Below the Plane on the Left
    rabove::Int64   # Atom ID of Atom Above the Plane on the Right
    rbelow::Int64   # Atom ID of Atom Below the Plane on the Right

    PCStereoPentagonalBiPyramid() = new()
end
StructTypes.StructType(::Type{PCStereoPentagonalBiPyramid}) = StructTypes.Mutable()

mutable struct PCStereoTShape
    center::Int64   # Atom ID of Atom Center
    top::Int64      # Atom ID of Atom In-Plane and at the Top
    bottom::Int64   # Atom ID of Atom In-Plane and at the Bottom
    above::Int64    # Atom ID of Atom Above the Plane

    PCStereoTShape() = new()
end
StructTypes.StructType(::Type{PCStereoTShape}) = StructTypes.Mutable()

mutable struct PCStereoTrigonalBiPyramid
    center::Int64   # Atom ID of Atom Center
    above::Int64    # Atom ID of Atom Above the Plane
    below::Int64    # Atom ID of Atom Below the Plane
    top::Int64      # Atom ID of Atom In-Plane and at the Top
    bottom::Int64   # Atom ID of Atom In-Plane and at the Bottom
    right::Int64    # Atom ID of Atom In-Plane and to the Right

    PCStereoTrigonalBiPyramid() = new()
end
StructTypes.StructType(::Type{PCStereoTrigonalBiPyramid}) = StructTypes.Mutable()

mutable struct PCStereoOctahedral
    center::Int64   # Atom ID of Atom Center
    top::Int64      # Atom ID of Atom In-Plane and at the Top
    bottom::Int64   # Atom ID of Atom In-Plane and at the Bottom
    labove::Int64   # Atom ID of Atom Above the Plane on the Left
    lbelow::Int64   # Atom ID of Atom Below the Plane on the Left
    rabove::Int64   # Atom ID of Atom Above the Plane on the Right
    rbelow::Int64   # Atom ID of Atom Below the Plane on the Right

    PCStereoOctahedral() = new()
end
StructTypes.StructType(::Type{PCStereoOctahedral}) = StructTypes.Mutable()

@enum PCStereoCenterSquarePlanarType CSP_USHAPE = 1 CSP_ZSHAPE CSP_XSHAPE CSP_ANY CSP_UNKNOWN = 255
StructTypes.StructType(::Type{PCStereoCenterSquarePlanarType}) = StructTypes.NumberType()
StructTypes.numbertype(::Type{PCStereoCenterSquarePlanarType}) = UInt8

mutable struct PCStereoSquarePlanar
    center::Int64       # Atom ID of Atom Center
    lbelow::Int64       # Atom ID of Left Below Plane Atom
    rbelow::Int64       # Atom ID of Right Below Plane Atom
    labove::Int64       # Atom ID of Left Above Plane Atom
    rabove::Int64       # Atom ID of Right Above Plane Atom

    parity::Union{Nothing,PCStereoCenterSquarePlanarType} # StereoCenter Type

    PCStereoSquarePlanar() = new()
end
StructTypes.StructType(::Type{PCStereoSquarePlanar}) = StructTypes.Mutable()

@enum PCStereoCenterPlanarDesignation CP_SAME = 1 CP_OPPOSITE CP_ANY CP_UNKNOWN = 255
StructTypes.StructType(::Type{PCStereoCenterPlanarDesignation}) = StructTypes.NumberType()
StructTypes.numbertype(::Type{PCStereoCenterPlanarDesignation}) = UInt8

@enum PCStereoCenterPlanarType CP_PLANAR = 1 CP_CUMULENIC
StructTypes.StructType(::Type{PCStereoCenterPlanarType}) = StructTypes.NumberType()
StructTypes.numbertype(::Type{PCStereoCenterPlanarType}) = UInt8

mutable struct PCStereoPlanar
    left::Int64       # Atom ID of Left Double Bond Atom
    ltop::Int64       # Atom ID of Top Atom attached to the Left Double Bond Atom
    lbottom::Int64    # Atom ID of Bottom Atom attached to the Left Double Bond Atom
    right::Int64      # Atom ID of Top Atom attached to the Right Double Bond Atom
    rtop::Int64       # Atom Identifier of Atom Below the Plane
    rbottom::Int64    # Atom ID of Bottom Atom attached to the Right Double Bond Atom

    parity::Union{Nothing,PCStereoCenterPlanarDesignation} # StereoCenter Designation
    type::Union{Nothing,PCStereoCenterPlanarType}          # Type of StereoCenter, SP2 Planar, if not specified

    PCStereoPlanar() = new()
end
StructTypes.StructType(::Type{PCStereoPlanar}) = StructTypes.Mutable()

@enum PCStereoCenterTetrahedralDesignation TP_CLOCKWISE = 1 TP_COUNTERCLOCKWISE TP_ANY TP_UNKNOWN = 255
StructTypes.StructType(::Type{PCStereoCenterTetrahedralDesignation}) = StructTypes.NumberType()
StructTypes.numbertype(::Type{PCStereoCenterTetrahedralDesignation}) = UInt8

@enum PCStereoCenterTetrahedralType TP_TETRAHEDRAL = 1 TP_CUMULENIC TP_BIARYL
StructTypes.StructType(::Type{PCStereoCenterTetrahedralType}) = StructTypes.NumberType()
StructTypes.numbertype(::Type{PCStereoCenterTetrahedralType}) = UInt8

mutable struct PCStereoTetrahedral
    center::Int64 # Atom Identifier of Atom Center
    above::Int64  # Atom Identifier of Atom Above the Plane
    top::Int64    # Atom Identifier of Atom In-Plane and at the Top
    bottom::Int64 # Atom Identifier of Atom In-Plane and at the Bottom
    below::Int64  # Atom Identifier of Atom Below the Plane

    parity::Union{Nothing,PCStereoCenterTetrahedralDesignation} # StereoCenter Designation
    type::Union{Nothing,PCStereoCenterTetrahedralType}          # Type of StereoCenter, Tetrahedral, if not specified

    PCStereoTetrahedral() = new()
end
StructTypes.StructType(::Type{PCStereoTetrahedral}) = StructTypes.Mutable()

mutable struct PCStereoCenter
    tetrahedral::Union{Nothing,PCStereoTetrahedral}        # Tetrahedral (SP3) StereoCenter
    planar::Union{Nothing,PCStereoPlanar}                  # Planar (SP2) StereoCenter
    squareplanar::Union{Nothing,PCStereoSquarePlanar}      # Square Planar (SP4) StereoCenter
    octahedral::Union{Nothing,PCStereoOctahedral}          # Octahedral (OC-6) / Square Pyramid (SPY-5) StereoCenters
    bipyramid::Union{Nothing,PCStereoTrigonalBiPyramid}    # Trigonal BiPyramid (TBPY-4 and TBPY-5) StereoCenters
    tshape::Union{Nothing,PCStereoTShape}                  # T-Shaped (TS-3) StereoCenters
    pentagonal::Union{Nothing,PCStereoPentagonalBiPyramid} # Pentagonal BiPyramid (PBPY-7) StereoCenters

    PCStereoCenter() = new()
end
StructTypes.StructType(::Type{PCStereoCenter}) = StructTypes.Mutable()

@enum PCBondType PCBOND_SINGLE = 1 PCBOND_DOUBLE PCBOND_TRIPLE PCBOND_QUADRUPLE PCBOND_DATIVE PCBOND_COMPLEX PCBOND_IONIC PCBOND_BOND_UNKNOWN = 255
StructTypes.StructType(::Type{PCBondType}) = StructTypes.NumberType()
StructTypes.numbertype(::Type{PCBondType}) = UInt8

mutable struct PCBonds
    aid1::Vector{Int64}
    aid2::Vector{Int64}
    order::Vector{PCBondType}

    PCBonds() = new()
end
StructTypes.StructType(::Type{PCBonds}) = StructTypes.Mutable()

mutable struct PCMMDBSource
    mmdb_id::Int64
    molecule_id::Int64
    molecule_name::Vector{String}
    residue_id::Union{Nothing,Int64}
    residue_name::Union{Nothing,String}
    atom_id::Union{Nothing,Int64}
    atom_name::Union{Nothing,String}

    PCMMDBSource() = new(0, 0, [], nothing, nothing, nothing, nothing)
end
StructTypes.StructType(::Type{PCMMDBSource}) = StructTypes.Mutable()

# we need to map the json names, which contain -, to _-using symbols
fieldnames = []
StructTypes.foreachfield((i, name, FT, v) -> push!(fieldnames, (name, Symbol(replace(String(name), "_" => "-")))), PCMMDBSource())
StructTypes.names(::Type{PCMMDBSource}) = Tuple(fieldnames)

mutable struct PCAtomSource
    aid::Int64
    source::PCMMDBSource

    PCAtomSource() = new()
end
StructTypes.StructType(::Type{PCAtomSource}) = StructTypes.Mutable()

@enum PCAtomRadicalType SINGLET = 1 DOUBLET TRIPLET QUARTET QUINTET HEXTET HEPTET OCTET NONE = 255
StructTypes.StructType(::Type{PCAtomRadicalType}) = StructTypes.NumberType()
StructTypes.numbertype(::Type{PCAtomRadicalType}) = UInt8

mutable struct PCAtomRadical
    aid::Int64
    value::PCAtomRadicalType

    PCAtomRadical() = new()
end
StructTypes.StructType(::Type{PCAtomRadical}) = StructTypes.Mutable()

# Type definitions used by JSON3, derived from https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/asn_spec/pcsubstance.asn.html
@enum(PCElement, 
      H = 1, HE, LI, BE, B, C, N, O, F, NE, NA, MG, AL, SI, P, S, CL, AR, K, CA, SC,
      TI, V, CR, MN, FE, CO, NI, CU, ZN, GA, GE, AS, SE, BR, KR, RB, SR, Y, ZR, NB,
      MO, TC, RU, RH, PD, AG, CD, IN, SN, SB, TE, I, XE, CS, BA, LA, CE, PR, ND, PM,
      SM, EU, GD, TB, DY, HO, ER, TM, YB, LU, HF, TA, W, RE, OS, IR, PT, AU, HG, TL,
      PB, BI, PO, AT, RN, FR, RA, AC, TH, PA, U, NP, PU, AM, CM, BK, CF, ES, FM, MD,
      NO, LR, RF, DB, SG, BH, HS, MT, DS, RG, CN, NH, FL, MC, LV, TS, OG,
      LP = 252, R = 253, D = 254, A = 255)
StructTypes.StructType(::Type{PCElement}) = StructTypes.NumberType()
StructTypes.numbertype(::Type{PCElement}) = UInt8

mutable struct PCAtoms
    aid::Vector{Int64}
    element::Vector{PCElement}

    label::Union{Nothing,Vector{PCAtomString}}
    isotope::Union{Nothing,Vector{PCAtomInt}}
    charge::Union{Nothing,Vector{PCAtomInt}}
    radical::Union{Nothing,Vector{PCAtomRadical}}
    source::Union{Nothing,Vector{PCAtomSource}}
    comment::Union{Nothing,Vector{PCAtomString}}

    PCAtoms() = new([], [], nothing, nothing, nothing, nothing, nothing, nothing)
end
StructTypes.StructType(::Type{PCAtoms}) = StructTypes.Mutable()

mutable struct PCCid
    cid::Int64

    PCCid() = new()
end
StructTypes.StructType(::Type{PCCid}) = StructTypes.Mutable()

mutable struct PCSid
    sid::Int64

    PCSid() = new()
end
StructTypes.StructType(::Type{PCSid}) = StructTypes.Mutable()

mutable struct PCXid
    xid::Int64

    PCXid() = new()
end
StructTypes.StructType(::Type{PCXid}) = StructTypes.Mutable()

@enum PCCompoundTypeID CTI_DEPOSITED = 0 CTI_STANDARDIZED CTI_COMPONENT CTI_NEUTRALIZED CTI_MIXTURE CTI_TAUTOMER CTI_PKASTATE CTI_UNKNOWN = 255
StructTypes.StructType(::Type{PCCompoundTypeID}) = StructTypes.NumberType()
StructTypes.numbertype(::Type{PCCompoundTypeID}) = UInt8

mutable struct PCCompoundType
    type::Union{Nothing,PCCompoundTypeID}
    id::Union{Nothing,PCCid,PCSid,PCXid}

    PCCompoundType() = new(CTI_UNKNOWN, nothing)
end
StructTypes.StructType(::Type{PCCompoundType}) = StructTypes.Mutable()

mutable struct PCCompound
    id::PCCompoundType
    atoms::Union{Nothing,PCAtoms}
    bonds::Union{Nothing,PCBonds}
    stereo::Union{Nothing,Vector{PCStereoCenter}}      # StereoCenter Descriptions
    coords::Union{Nothing,Vector{PCCoordinates}}       # 2D/3D Coordinate Sets of Compound
    charge::Union{Nothing,Int32}                       # Provided Total Formal Charge  (Signed Integer)
    props::Union{Nothing,Vector{PCInfoData}}           # Derived (computed) Properties
    stereogroups::Union{Nothing,Vector{PCStereoGroup}} # Relative stereochemistry groups
    count::Union{Nothing,PCCount}                      # Counts of various properties
    vbalt::Union{Nothing,Vector{PCCompound}}           # Alternate Valence-Bond Forms

    PCCompound() = new(PCCompoundType(), nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing)
end
StructTypes.StructType(::Type{PCCompound}) = StructTypes.Mutable()

mutable struct PCResult
    PC_Compounds::Vector{PCCompound}

    PCResult() = new()
end
StructTypes.StructType(::Type{PCResult}) = StructTypes.Mutable()

function convert_coordinates(pb_coords_vec::Vector{PCCoordinates})
    if isnothing(pb_coords_vec)
        return []
    end

    result = []

    for pb_coords in pb_coords_vec
        # first, figure out what unit to use for the Coordinates
        # we interpret pixels, points, and stdbonds as Angstroms
        scale = CT_UNITS_NANOMETERS in pb_coords.type ? 10. : 1.

        # do we have 3d information, or only 2d?
        is_3d = CT_3D in pb_coords.type

        for c in pb_coords.conformers
            converted = Array{Vector3}(undef, length(c.x))

            for i in 1:length(c.x)
                converted[i] = Vector3(c.x[i], c.y[i], (is_3d && !isnothing(c.z)) ? c.z[i] : 0.0) * scale
            end

            push!(result, converted)
        end
    end

    result
end

# TODO: 
#   - full conversion, adding all properties

# NOTE: conformers are stored as frames
function load_pubchem_json(fname::String, T=Float32)
    pb = JSON3.read(read(fname, String), PCResult)

    # for now, use the file name as the name for the molecule
    mol = Molecule(fname)

    for compound in pb.PC_Compounds

        if !isnothing(compound.atoms) && !isnothing(compound.coords)
            conformers = convert_coordinates(compound.coords)

            for i in 1:length(compound.atoms.aid)
                for j in 1:length(conformers)
                    # Note: the atom will be assigned an id in add_atom!
                    atom = (number=compound.atoms.aid[i],
                            name="",
                            element = isnothing(compound.atoms.element) 
                                ? Elements.Unknown 
                                : ElementType(Int(compound.atoms.element[i])),
                            atomtype = isnothing(compound.atoms.label)
                                ? ""
                                : compound.atoms.label[i].value, # does the label contain the atom type?
                            r = T.(conformers[j][i]),
                            v = Vector3(T(0.), T(0.), T(0.)),
                            F = Vector3(T(0.), T(0.), T(0.)),
                            has_velocity = false,
                            has_force = false,
                            frame_id = j
                    )

                    push!(mol, atom)
                end
            end
        end

        if !isnothing(compound.bonds)
            for i in 1:length(compound.bonds.aid1)
                order = Int(compound.bonds.order[i])

                b = (a1 = compound.bonds.aid1[i], 
                     a2 = compound.bonds.aid2[i],
                     order = (order <= 4) ? BondOrderType(order) : BondOrder.Unknown
                    )

                push!(mol, b)
            end
        end
    end
    
    mol
end
