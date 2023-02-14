using EnumX

export BondOrder, BondOrderType

@enumx BondOrder begin
    Single = 1
    Double = 2
    Triple = 3
    Quadruple = 4
    Unknown = 100
end

const BondOrderType = BondOrder.T
