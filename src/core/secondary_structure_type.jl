export
    SecondaryStructureElement,
    SecondaryStructureType

@enumx SecondaryStructureElement begin
    Helix   = 1
    Coil    = 2
    Strand  = 3
    Turn    = 4
    Unknown = 5
end

const SecondaryStructureType = SecondaryStructureElement.T
