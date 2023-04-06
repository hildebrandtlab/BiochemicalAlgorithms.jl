module Constants

using Unitful
using PhysicalConstants

export ε₀, ε₀_nounit, N_A_nounit, e₀, e₀_nounit

const ε₀ = PhysicalConstants.CODATA2018.VacuumElectricPermittivity
const ε₀_nounit = ustrip(ε₀)

const N_A = PhysicalConstants.CODATA2018.AvogadroConstant
const N_A_nounit = ustrip(N_A)

const e₀ = PhysicalConstants.CODATA2018.ElementaryCharge
const e₀_nounit = ustrip(e₀)

end