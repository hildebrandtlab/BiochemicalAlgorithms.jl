
using EnumX

export DefBond, DefBonds

@enumx DefBonds begin
    sb = 1
    db = 2
    tb = 3
    qd = 4
    qt = 5
    am = 20     # mol2 bond type, amide
    ar = 21     # mol2 bond type, aromatic
    du = 23     # mol2 bond type, dummy
    un = 100    # mol2 bond type, unknown
    nc = 101    # mol2 bond type, not connected
end

const DefBond = DefBonds.T

Base.parse(DefBonds, x) = getproperty(DefBonds, Symbol(x))