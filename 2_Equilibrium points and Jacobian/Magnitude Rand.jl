using Flows
using ToySystems
using ToySystems.NineModeSystemEq
using Optim
using PyPlot; pygui(true)
include("9msctrl.jl")
include("Initial condition variation.jl")
include("Perturbation time.jl")
include("Contour Plots.jl")
import LinearAlgebra: norm

function mgscale(u0::AbstractVector, scale)
    u = (1/norm(u0)).*u0
    u_final = scale.*u
    return u_final
end