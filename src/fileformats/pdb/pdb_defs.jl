using AutoHashEquals
using DataStructures

using BiochemicalAlgorithms: 
    Molecule,
    Chain,
    Fragment

@enum(PDBRecordType,
    RECORD_TYPE__UNKNOWN = 0,
    RECORD_TYPE__ANISOU,
    RECORD_TYPE__ATOM,
    RECORD_TYPE__AUTHOR,
    RECORD_TYPE__CAVEAT,
    RECORD_TYPE__CISPEP,
    RECORD_TYPE__COMPND,
    RECORD_TYPE__CONECT,
    RECORD_TYPE__CON06,
    RECORD_TYPE__CON061,
    RECORD_TYPE__CON062,
    RECORD_TYPE__CON063,
    RECORD_TYPE__CON064,
    RECORD_TYPE__CRYST1,
    RECORD_TYPE__DBREF,
    RECORD_TYPE__END,
    RECORD_TYPE__ENDMDL,
    RECORD_TYPE__EXPDTA,
    RECORD_TYPE__FORMUL,
    RECORD_TYPE__FTNOTE,
    RECORD_TYPE__HEADER,
    RECORD_TYPE__HELIX,
    RECORD_TYPE__HET,
    RECORD_TYPE__HETATM,
    RECORD_TYPE__HETNAM,
    RECORD_TYPE__HETSYN,
    RECORD_TYPE__HYDBND,
    RECORD_TYPE__JRNL,
    RECORD_TYPE__KEYWDS,
    RECORD_TYPE__LINK,
    RECORD_TYPE__MASTER,
    RECORD_TYPE__MODEL,
    RECORD_TYPE__MODRES,
    RECORD_TYPE__MTRIX1,
    RECORD_TYPE__MTRIX2,
    RECORD_TYPE__MTRIX3,
    RECORD_TYPE__OBSLTE,
    RECORD_TYPE__ORIGX1,
    RECORD_TYPE__ORIGX2,
    RECORD_TYPE__ORIGX3,
    RECORD_TYPE__REMARK,
    RECORD_TYPE__REVDAT,
    RECORD_TYPE__SCALE1,
    RECORD_TYPE__SCALE2,
    RECORD_TYPE__SCALE3,
    RECORD_TYPE__SEQADV,
    RECORD_TYPE__SEQRES,
    RECORD_TYPE__SHEET,
    RECORD_TYPE__SIGATM,
    RECORD_TYPE__SIGUIJ,
    RECORD_TYPE__SITE,
    RECORD_TYPE__SLTBRG,
    RECORD_TYPE__SOURCE,
    RECORD_TYPE__SPRSDE,
    RECORD_TYPE__SSBOND,
    RECORD_TYPE__TER,
    RECORD_TYPE__TITLE,
    RECORD_TYPE__TURN,
    RECORD_TYPE__TVECT,
        
    NUMBER_OF_REGISTERED_RECORD_TYPES,

    ALL_RECORD_TYPES
)

@auto_hash_equals struct RecordTypeFormat
    record_type::PDBRecordType
    tag::String
    format_string::String
end

@auto_hash_equals struct PDBRecord
    type::PDBRecordType
    data::Tuple
end

@auto_hash_equals struct UniqueResidueID
    name::String
    chain_id::String
    number::Int
    insertion_code::String
end

@auto_hash_equals struct SSBondRecord
    number::Int
    first::UniqueResidueID
    second::UniqueResidueID
end

@auto_hash_equals struct HelixRecord
    number::Int
    name::String
    initial_residue::UniqueResidueID
    terminal_residue::UniqueResidueID
    helix_class::Int
    comment::String
end

@auto_hash_equals struct SheetRecord
    number::Int
    name::String
    initial_residue::UniqueResidueID
    terminal_residue::UniqueResidueID
    sense_of_strand::Bool
end

@auto_hash_equals struct TurnRecord
    number::Int
    name::String
    initial_residue::UniqueResidueID
    terminal_residue::UniqueResidueID
    comment::String
end

@auto_hash_equals mutable struct PDBInfo{T}
    name::String
    deposition_date::String
    id::String

    title::String

    records::Deque{PDBRecord}

    ssbonds::Deque{SSBondRecord}

    secondary_structures::Deque{Union{HelixRecord, SheetRecord, TurnRecord}}

    selected_model::Int
    current_model::Int

    current_chain::Union{Chain{T}, Nothing}
    current_residue::Union{Fragment{T}, Nothing}

    alternate_location_identifier::String

    function PDBInfo{T}(selected_model=-1) where {T}
        new("", "", "", "", Deque{PDBRecord}(), Deque{SSBondRecord}(), 
            Deque{Union{HelixRecord, SheetRecord, TurnRecord}}(), selected_model, 
            1, nothing, nothing, "A")
    end
end


const FORMAT_UNKNOWN          = ""
const FORMAT_ANISOU           = "%5ld %-4.4s%c%3.3s %c%4ld%c %7ld%7ld%7ld%7ld%7ld%7ld  %4.4s%2.2s%2.2s"
const FORMAT_ATOM             = "%5ld %-4.4s%c%3.3s %c%4ld%c   %8.3f%8.3f%8.3f%6.2f%6.2f      %4.4s%2.2s%2.2s"
const FORMAT_ATOM_PARTIAL_CRG = "%5ld %-4.4s%c%3.3s %c%4ld%c   %8.3f%8.3f%8.3f%6.2f%6.2f      %4.4s%4.4s"
const FORMAT_AUTHOR           = "  %2ld%-60.60s"
const FORMAT_CAVEAT           = "  %2ld %4.4s    %51.51s"
const FORMAT_CISPEP           = " %3ld %3.3s %c %4ld%c   %3.3s %c %4ld%c       %3ld       %6f"
const FORMAT_COMPND           = "  %2ld%-60.60s"
const FORMAT_CON06            = "%5ld%5ld%5ld%5ld%5ld"
const FORMAT_CON061           = "%5ld"
const FORMAT_CON062           = "%5ld%5ld"
const FORMAT_CON063           = "%5ld%5ld%5ld"
const FORMAT_CON064           = "%5ld%5ld%5ld%5ld"
const FORMAT_CONECT           = "%5ld%5ld%5ld%5ld%5ld%5ld%5ld%5ld%5ld%5ld%5ld"
const FORMAT_CRYST1           = "%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11.11s%4ld"
const FORMAT_DBREF            = " %4.4s %c %4ld%c %4ld%c %6.6s %8.8s %12.12s %5ld%c %5ld%c"
const FORMAT_END              = ""
const FORMAT_ENDMDL           = ""
const FORMAT_EXPDTA           = "  %2ld%60.60s"
const FORMAT_FORMUL           = "  %2ld  %3.3s %2ld%c%51.51s"
const FORMAT_FTNOTE           = " %3ld %59.59s"
const FORMAT_HEADER           = "    %-40.40s%9.9s   %4.4s"
const FORMAT_HELIX            = " %3ld %3.3s %3.3s %c %4ld%c %3.3s %c %4ld%c%2ld%30.30s %5ld"
const FORMAT_HET              = " %3.3s  %c%4ld%c  %5ld  %43.43s"
const FORMAT_HETATM           = "%5ld %-4.4s%c%3.3s %c%4ld%c   %8.3f%8.3f%8.3f%6.2f%6.2f      %4.4s%2.2s%2.2s"
const FORMAT_HETNAM           = "  %2ld %3.3s %55.55s"
const FORMAT_HETSYN           = "  %2ld %3.3s %55.55s"
const FORMAT_HYDBND           = "      %4.4s%c%3.3s %c%5ld%c %4.4s%c %c%5ld%c %4.4s%c%3.3s %c%5ld%c%6ld %6ld"
const FORMAT_JRNL             = "      %58.58s"
const FORMAT_KEYWDS           = "  %2ld%60.60s"
const FORMAT_LINK             = "      %4.4s%c%3.3s %c%4ld%c               %4.4s%c%3.3s %c%4ld%c  %6ld %6ld"
const FORMAT_MASTER           = "    %5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d"
const FORMAT_MODEL            = "    %4ld"
const FORMAT_MODRES           = " %3.3s %3.3s %c %4ld%c %3.3s  %41.41s"
const FORMAT_MTRIX1           = " %3ld%10.6f%10.6f%10.6f     %10.5f    %1ld"
const FORMAT_MTRIX2           = " %3ld%10.6f%10.6f%10.6f     %10.5f    %1ld"
const FORMAT_MTRIX3           = " %3ld%10.6f%10.6f%10.6f     %10.5f    %1ld"
const FORMAT_OBSLTE           = "  %2ld %9.9s %4.4s      %4.4s %4.4s %4.4s %4.4s %4.4s %4.4s %4.4s %4.4s"
const FORMAT_ORIGX1           = "    %10f%10f%10f     %10f"
const FORMAT_ORIGX2           = "    %10f%10f%10f     %10f"
const FORMAT_ORIGX3           = "    %10f%10f%10f     %10f"
const FORMAT_REMARK           = " %3ld %-59.59s"
const FORMAT_REVDAT           = " %3ld%2ld %9.9s %5.5s   %1ld       %6.6s %6.6s %6.6s %6.6s"
const FORMAT_SCALE1           = "    %10f%10f%10f     %10f"
const FORMAT_SCALE2           = "    %10f%10f%10f     %10f"
const FORMAT_SCALE3           = "    %10f%10f%10f     %10f"
const FORMAT_SEQADV           = " %4.4s %4.4s %c %4ld%c %4.4s %9.9s %3.3s %5ld %21.21s"
const FORMAT_SEQRES           = "  %2ld %c %4ld  %3.3s %3.3s %3.3s %3.3s %3.3s %3.3s %3.3s %3.3s %3.3s %3.3s %3.3s %3.3s %3.3s"
const FORMAT_SHEET            = " %3ld %3.3s%2ld %3.3s %c%4ld%c %3.3s %c%4ld%c%2ld %-4.4s%3.3s %c%4ld%c %-4.4s%3.3s %c%4ld%c"
const FORMAT_SIGATM           = "%5ld %4.4s%c%3.3s %c%4ld%c   %8f%8f%8f%6f%6f      %4.4s%2.2s%2.2s"
const FORMAT_SIGUIJ           = "%5ld %-4.4s%c%3.3s %c%4ld%c %7ld%7ld%7ld%7ld%7ld%7ld  %4.4s%2.2s%2.2s"
const FORMAT_SITE             = " %3ld %3.3s %2ld %3.3s %c%4ld%c %3.3s %c%4ld%c %3.3s %c%4ld%c %3.3s %c%4ld%c"
const FORMAT_SLTBRG           = "      %4.4s%c%3.3s %c%4ld%c               %4.4s%c%3.3s %c%4ld%c  %6ld%6ld"
const FORMAT_SOURCE           = "  %2ld%-60.60s"
const FORMAT_SPRSDE           = "  %2ld %9.9s %4.4s %4.4s %4.4s %4.4s %4.4s %4.4s %4.4s %4.4s %4.4s "
const FORMAT_SSBOND           = " %3ld %3.3s %c %4ld%c   %3.3s %c %4ld%c                       %6ld %6ld"
const FORMAT_TER              = "%5ld      %3.3s %c%4ld%c"
const FORMAT_TITLE            = "  %2ld%60.60s"
const FORMAT_TURN             = " %3ld %3.3s %3.3s %c%4ld%c %3.3s %c%4ld%c    %-30.30s"
const FORMAT_TVECT            = " %3ld%10f%10f%10f%30.30s"

const RECORD_TAG_UNKNOWN = "     "
const RECORD_TAG_ANISOU  = "ANISOU"
const RECORD_TAG_ATOM    = "ATOM  "
const RECORD_TAG_AUTHOR  = "AUTHOR"
const RECORD_TAG_CAVEAT  = "CAVEAT"
const RECORD_TAG_CISPEP  = "CISPEP"
const RECORD_TAG_COMPND  = "COMPND"
const RECORD_TAG_CON061  = "CONECT"
const RECORD_TAG_CON062  = "CONECT"
const RECORD_TAG_CON063  = "CONECT"
const RECORD_TAG_CON064  = "CONECT"
const RECORD_TAG_CON06   = "CONECT"
const RECORD_TAG_CONECT  = "CONECT"
const RECORD_TAG_CRYST1  = "CRYST1"
const RECORD_TAG_DBREF   = "DBREF "
const RECORD_TAG_END     = "END   "
const RECORD_TAG_ENDMDL  = "ENDMDL"
const RECORD_TAG_EXPDTA  = "EXPDTA"
const RECORD_TAG_FORMUL  = "FORMUL"
const RECORD_TAG_FTNOTE  = "FTNOTE"
const RECORD_TAG_HEADER  = "HEADER"
const RECORD_TAG_HELIX   = "HELIX "
const RECORD_TAG_HET     = "HET   "
const RECORD_TAG_HETATM  = "HETATM"
const RECORD_TAG_HETNAM  = "HETNAM"
const RECORD_TAG_HETSYN  = "HETSYN"
const RECORD_TAG_HYDBND  = "HYDBND"
const RECORD_TAG_JRNL    = "JRNL  "
const RECORD_TAG_KEYWDS  = "KEYWDS"
const RECORD_TAG_LINK    = "LINK  "
const RECORD_TAG_MASTER  = "MASTER"
const RECORD_TAG_MODEL   = "MODEL "
const RECORD_TAG_MODRES  = "MODRES"
const RECORD_TAG_MTRIX1  = "MTRIX1"
const RECORD_TAG_MTRIX2  = "MTRIX2"
const RECORD_TAG_MTRIX3  = "MTRIX3"
const RECORD_TAG_OBSLTE  = "OBSLTE"
const RECORD_TAG_ORIGX1  = "ORIGX1"
const RECORD_TAG_ORIGX2  = "ORIGX2"
const RECORD_TAG_ORIGX3  = "ORIGX3"
const RECORD_TAG_REMARK  = "REMARK"
const RECORD_TAG_REVDAT  = "REVDAT"
const RECORD_TAG_SCALE1  = "SCALE1"
const RECORD_TAG_SCALE2  = "SCALE2"
const RECORD_TAG_SCALE3  = "SCALE3"
const RECORD_TAG_SEQADV  = "SEQADV"
const RECORD_TAG_SEQRES  = "SEQRES"
const RECORD_TAG_SHEET   = "SHEET "
const RECORD_TAG_SIGATM  = "SIGATM"
const RECORD_TAG_SIGUIJ  = "SIGUIJ"
const RECORD_TAG_SITE    = "SITE  "
const RECORD_TAG_SLTBRG  = "SLTBRG"
const RECORD_TAG_SOURCE  = "SOURCE"
const RECORD_TAG_SPRSDE  = "SPRSDE"
const RECORD_TAG_SSBOND  = "SSBOND"
const RECORD_TAG_TER     = "TER   "
const RECORD_TAG_TITLE   = "TITLE "
const RECORD_TAG_TURN    = "TURN  "
const RECORD_TAG_TVECT   = "TVECT "

const RECORD_TYPE_FORMAT = [
    RecordTypeFormat(RECORD_TYPE__UNKNOWN, RECORD_TAG_UNKNOWN, FORMAT_UNKNOWN),
    RecordTypeFormat(RECORD_TYPE__ANISOU,  RECORD_TAG_ANISOU,  FORMAT_ANISOU),
    RecordTypeFormat(RECORD_TYPE__ATOM,    RECORD_TAG_ATOM,    FORMAT_ATOM),
    RecordTypeFormat(RECORD_TYPE__AUTHOR,  RECORD_TAG_AUTHOR,  FORMAT_AUTHOR),
    RecordTypeFormat(RECORD_TYPE__CAVEAT,  RECORD_TAG_CAVEAT,  FORMAT_CAVEAT),
    RecordTypeFormat(RECORD_TYPE__CISPEP,  RECORD_TAG_CISPEP,  FORMAT_CISPEP),
    RecordTypeFormat(RECORD_TYPE__COMPND,  RECORD_TAG_COMPND,  FORMAT_COMPND),
    RecordTypeFormat(RECORD_TYPE__CON06,   RECORD_TAG_CON06,   FORMAT_CON06),
    RecordTypeFormat(RECORD_TYPE__CON061,  RECORD_TAG_CON061,  FORMAT_CON061),
    RecordTypeFormat(RECORD_TYPE__CON062,  RECORD_TAG_CON062,  FORMAT_CON062),
    RecordTypeFormat(RECORD_TYPE__CON063,  RECORD_TAG_CON063,  FORMAT_CON063),
    RecordTypeFormat(RECORD_TYPE__CON064,  RECORD_TAG_CON064,  FORMAT_CON064),
    RecordTypeFormat(RECORD_TYPE__CONECT,  RECORD_TAG_CONECT,  FORMAT_CONECT),
    RecordTypeFormat(RECORD_TYPE__CRYST1,  RECORD_TAG_CRYST1,  FORMAT_CRYST1),
    RecordTypeFormat(RECORD_TYPE__DBREF,   RECORD_TAG_DBREF,   FORMAT_DBREF),
    RecordTypeFormat(RECORD_TYPE__END,     RECORD_TAG_END,     FORMAT_END),
    RecordTypeFormat(RECORD_TYPE__ENDMDL,  RECORD_TAG_ENDMDL,  FORMAT_ENDMDL),
    RecordTypeFormat(RECORD_TYPE__EXPDTA,  RECORD_TAG_EXPDTA,  FORMAT_EXPDTA),
    RecordTypeFormat(RECORD_TYPE__FORMUL,  RECORD_TAG_FORMUL,  FORMAT_FORMUL),
    RecordTypeFormat(RECORD_TYPE__FTNOTE,  RECORD_TAG_FTNOTE,  FORMAT_FTNOTE),
    RecordTypeFormat(RECORD_TYPE__HEADER,  RECORD_TAG_HEADER,  FORMAT_HEADER),
    RecordTypeFormat(RECORD_TYPE__HELIX,   RECORD_TAG_HELIX,   FORMAT_HELIX),
    RecordTypeFormat(RECORD_TYPE__HET,     RECORD_TAG_HET,     FORMAT_HET),
    RecordTypeFormat(RECORD_TYPE__HETATM,  RECORD_TAG_HETATM,  FORMAT_HETATM),
    RecordTypeFormat(RECORD_TYPE__HETNAM,  RECORD_TAG_HETNAM,  FORMAT_HETNAM),
    RecordTypeFormat(RECORD_TYPE__HETSYN,  RECORD_TAG_HETSYN,  FORMAT_HETSYN),
    RecordTypeFormat(RECORD_TYPE__HYDBND,  RECORD_TAG_HYDBND,  FORMAT_HYDBND),
    RecordTypeFormat(RECORD_TYPE__JRNL,    RECORD_TAG_JRNL,    FORMAT_JRNL),
    RecordTypeFormat(RECORD_TYPE__KEYWDS,  RECORD_TAG_KEYWDS,  FORMAT_KEYWDS),
    RecordTypeFormat(RECORD_TYPE__LINK,    RECORD_TAG_LINK,    FORMAT_LINK),
    RecordTypeFormat(RECORD_TYPE__MASTER,  RECORD_TAG_MASTER,  FORMAT_MASTER),
    RecordTypeFormat(RECORD_TYPE__MODEL,   RECORD_TAG_MODEL,   FORMAT_MODEL),
    RecordTypeFormat(RECORD_TYPE__MODRES,  RECORD_TAG_MODRES,  FORMAT_MODRES),
    RecordTypeFormat(RECORD_TYPE__MTRIX1,  RECORD_TAG_MTRIX1,  FORMAT_MTRIX1),
    RecordTypeFormat(RECORD_TYPE__MTRIX2,  RECORD_TAG_MTRIX2,  FORMAT_MTRIX2),
    RecordTypeFormat(RECORD_TYPE__MTRIX3,  RECORD_TAG_MTRIX3,  FORMAT_MTRIX3),
    RecordTypeFormat(RECORD_TYPE__OBSLTE,  RECORD_TAG_OBSLTE,  FORMAT_OBSLTE),
    RecordTypeFormat(RECORD_TYPE__ORIGX1,  RECORD_TAG_ORIGX1,  FORMAT_ORIGX1),
    RecordTypeFormat(RECORD_TYPE__ORIGX2,  RECORD_TAG_ORIGX2,  FORMAT_ORIGX2),
    RecordTypeFormat(RECORD_TYPE__ORIGX3,  RECORD_TAG_ORIGX3,  FORMAT_ORIGX3),
    RecordTypeFormat(RECORD_TYPE__REMARK,  RECORD_TAG_REMARK,  FORMAT_REMARK),
    RecordTypeFormat(RECORD_TYPE__REVDAT,  RECORD_TAG_REVDAT,  FORMAT_REVDAT),
    RecordTypeFormat(RECORD_TYPE__SCALE1,  RECORD_TAG_SCALE1,  FORMAT_SCALE1),
    RecordTypeFormat(RECORD_TYPE__SCALE2,  RECORD_TAG_SCALE2,  FORMAT_SCALE2),
    RecordTypeFormat(RECORD_TYPE__SCALE3,  RECORD_TAG_SCALE3,  FORMAT_SCALE3),
    RecordTypeFormat(RECORD_TYPE__SEQADV,  RECORD_TAG_SEQADV,  FORMAT_SEQADV),
    RecordTypeFormat(RECORD_TYPE__SEQRES,  RECORD_TAG_SEQRES,  FORMAT_SEQRES),
    RecordTypeFormat(RECORD_TYPE__SHEET,   RECORD_TAG_SHEET,   FORMAT_SHEET),
    RecordTypeFormat(RECORD_TYPE__SIGATM,  RECORD_TAG_SIGATM,  FORMAT_SIGATM),
    RecordTypeFormat(RECORD_TYPE__SIGUIJ,  RECORD_TAG_SIGUIJ,  FORMAT_SIGUIJ),
    RecordTypeFormat(RECORD_TYPE__SITE,    RECORD_TAG_SITE,    FORMAT_SITE),
    RecordTypeFormat(RECORD_TYPE__SLTBRG,  RECORD_TAG_SLTBRG,  FORMAT_SLTBRG),
    RecordTypeFormat(RECORD_TYPE__SOURCE,  RECORD_TAG_SOURCE,  FORMAT_SOURCE),
    RecordTypeFormat(RECORD_TYPE__SPRSDE,  RECORD_TAG_SPRSDE,  FORMAT_SPRSDE),
    RecordTypeFormat(RECORD_TYPE__SSBOND,  RECORD_TAG_SSBOND,  FORMAT_SSBOND),
    RecordTypeFormat(RECORD_TYPE__TER,     RECORD_TAG_TER,     FORMAT_TER),
    RecordTypeFormat(RECORD_TYPE__TITLE,   RECORD_TAG_TITLE,   FORMAT_TITLE),
    RecordTypeFormat(RECORD_TYPE__TURN,    RECORD_TAG_TURN,    FORMAT_TURN),
    RecordTypeFormat(RECORD_TYPE__TVECT,   RECORD_TAG_TVECT,   FORMAT_TVECT)
]

# TODO: Better idea for just the tag
const RECORD_MAP = Dict(
    rtf.tag => rtf for rtf in RECORD_TYPE_FORMAT
)
