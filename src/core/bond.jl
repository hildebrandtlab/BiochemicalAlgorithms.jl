export Bond, BondOrder, BondOrderType

using EnumX

@enumx BondOrder begin
    Single = 1
    Double = 2
    Triple = 3
    Quadruple = 4
    Quintuple = 5
    Amide = 20      # mol2 format bond type
    Aromatic = 21   # mol2 format bond type
    Dummy = 23      # mol2 format bond type
    Unknown = 100
    NotConnected = 101  # mol2 format bond type
end

const BondOrderType = BondOrder.T

Base.parse(BondOrder, x) = getproperty(BondOrder, Symbol(x))

const Bond = @NamedTuple begin
    a1::Int
    a2::Int
    order::BondOrderType
end