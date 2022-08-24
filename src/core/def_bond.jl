
using EnumX

export DefBond, DefBonds

@enumx DefBonds begin
    sb = 1
    db = 2
    tb = 3
    Unknown = 100
end

const DefBond = DefBonds.T

Base.parse(DefBond, x) = getproperty(DefBond, Symbol(x))