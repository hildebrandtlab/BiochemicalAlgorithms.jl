
using EnumX

export BondDef, BondDefs

@enumx BondDefs begin
    sb = 1
    db = 2
    tp = 3
    Unknown = 100
end

const BondDef = BondDefs.T

Base.parse(BondDef, s) = getproperty(BondDef, Symbol(s))