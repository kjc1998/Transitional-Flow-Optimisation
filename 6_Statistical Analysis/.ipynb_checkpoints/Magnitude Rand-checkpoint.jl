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

function mgscale(u0::AbstractVector, scale, ISRANGE = false)
    u = (1/norm(u0)).*u0
    if ISRANGE==false
        u_final = scale.*u
    else
        u_final = (rand()*scale).*u
    end
    return u_final
end