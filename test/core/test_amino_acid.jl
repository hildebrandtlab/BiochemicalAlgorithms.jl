using BioSymbols

@testset "AminoAcid" begin
    # can we break the constructor of BioSymbols:
    @test_throws InexactError AminoAcid(' ') # call BioSymbols.AminoAcid(::Char)
    @test_throws MethodError AminoAcid()
    @test_throws MethodError AminoAcid(0)

    # now generate amino acids:
    ala = AA_A
    ala2 = AminoAcid('A')  
    
    @test ala isa AminoAcid
    @test name(ala) == name(ala2) == "Alanine"
    @test one_letter_code(ala) == "A"
    @test three_letter_code(ala) == "ALA"
    
    # test all possibilities
    n = ["Alanine", "Arginine", "Asparagine", "Aspartate",
        "Cysteine", "Glutamine", "Glutamate", "Glycine", "Histidine",
        "Isoleucine", "Leucine", "Lysine", "Methionine", "Phenylalanine",
        "Proline", "Serine", "Threonine", "Tryptophan", "Tyrosine",
        "Valine", "Pyrrolysine", "Selenocysteine", "Asparagine or Aspartate", 
        "Leucine or Isoleucine", "Glutamine or Glutamate", "Any / Unknown", 
        "Termination codon", "Gap"]

    t = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
        "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "CYS", 
        "PRO", "SER", "THR", "TRP", "TYR", "VAL", "PYL", 
        "SEC", "ASX", "XLE", "GLX", "XAA", "TER", "GAP"]

    o = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", 
        "I", "L", "K", "M", "F", "P", "S", "T", "W", 
        "Y", "V", "O", "U", "B", "J", "Z", "X", "*", "---"]


    for i in eachindex(n)
        a = AminoAcid(o[i][1])
        @test name(a) == n[i]
        @test three_letter_code(a) == t[i]
    end
end
