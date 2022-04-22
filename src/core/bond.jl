export Bond

using EnumX

@enumx BondOrder begin
    Single = 1
    Double = 2
    Triple = 3
    Quadruple = 4
    Unknown = 100
end

const BondOrderType = BondOrder.T

const Bond = @NamedTuple begin
    a1::Int
    a2::Int
    order::BondOrderType
end