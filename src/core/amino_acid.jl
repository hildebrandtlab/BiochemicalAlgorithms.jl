
using BioSymbols

export AminoAcid, name, three_letter_code, one_letter_code

const AminoAcid = BioSymbols.AminoAcid

const AminoAcidDetails = @NamedTuple begin
    name::String
    three_letter_code::String
    one_letter_code::String
end

const AminoAcidProperties = Dict{AminoAcid, AminoAcidDetails}(
    AA_A => (name="Alanine",        three_letter_code="ALA", one_letter_code="A"),
    AA_R => (name="Arginine",       three_letter_code="ARG", one_letter_code="R"),
    AA_N => (name="Asparagine",     three_letter_code="ASN", one_letter_code="N"),
    AA_D => (name="Aspartate",      three_letter_code="ASP", one_letter_code="D"),
    AA_C => (name="Cysteine",       three_letter_code="CYS", one_letter_code="C"),
    AA_Q => (name="Glutamine",      three_letter_code="GLN", one_letter_code="Q"),
    AA_E => (name="Glutamate",      three_letter_code="GLU", one_letter_code="E"),
    AA_G => (name="Glycine",        three_letter_code="GLY", one_letter_code="G"),
    AA_H => (name="Histidine",      three_letter_code="HIS", one_letter_code="H"),
    AA_I => (name="Isoleucine",     three_letter_code="ILE", one_letter_code="I"),
    AA_L => (name="Leucine",        three_letter_code="LEU", one_letter_code="L"),
    AA_K => (name="Lysine",         three_letter_code="LYS", one_letter_code="K"),
    AA_M => (name="Methionine",     three_letter_code="MET", one_letter_code="M"),
    AA_F => (name="Phenylalanine",  three_letter_code="CYS", one_letter_code="F"),
    AA_P => (name="Proline",        three_letter_code="PRO", one_letter_code="P"),
    AA_S => (name="Serine",         three_letter_code="SER", one_letter_code="S"),
    AA_T => (name="Threonine",      three_letter_code="THR", one_letter_code="T"),
    AA_W => (name="Tryptophan",     three_letter_code="TRP", one_letter_code="W"),
    AA_Y => (name="Tyrosine",       three_letter_code="TYR", one_letter_code="Y"),
    AA_V => (name="Valine",         three_letter_code="VAL", one_letter_code="V"),
    AA_O => (name="Pyrrolysine",    three_letter_code="PYL", one_letter_code="O"),
    AA_U => (name="Selenocysteine", three_letter_code="SEC", one_letter_code="U"),

    AA_B => (name="Asparagine or Aspartate",
            three_letter_code="ASX", 
            one_letter_code="B"),

    AA_J => (name="Leucine or Isoleucine",
            three_letter_code="XLE",
            one_letter_code="J"),

    AA_Z => (name="Glutamine or Glutamate",
            three_letter_code="GLX",
            one_letter_code="Z"),

    AA_X => (name="Any / Unknown", three_letter_code="XAA", one_letter_code="X"),
    
    AA_Term => (name="Termination codon", three_letter_code="TER", one_letter_code="*"),
    
    AA_Gap => (name="Gap", three_letter_code="GAP", one_letter_code="---")
 )

### Functions

name(aa::AminoAcid) = AminoAcidProperties[aa].name
three_letter_code(aa::AminoAcid) = AminoAcidProperties[aa].three_letter_code
one_letter_code(aa::AminoAcid) = AminoAcidProperties[aa].one_letter_code