using Flows
using ToySystems
using ToySystems.NineModeSystemEq
using Optim
include("9msctrl.jl")
include("Initial condition variation.jl")
include("Perturbation time.jl")
include("Contour Plots.jl")
include("Magnitude Rand.jl")
include("Objective function Monitor.jl")
import LinearAlgebra: norm
using PyPlot; pygui(true)
using PyCall
np = pyimport("numpy")
sc = pyimport("scipy")
plt = pyimport("matplotlib.pyplot");

function Inipointnum(n, scale = false, ISRANGE = false, d = zeros(9))
    out = []
    if scale == false
        for i in 1:n
            out = push!(out, rand(9) + d)
        end
    else
        for i in 1:n
            out = push!(out, mgscale(rand(9),scale,ISRANGE) + d)
        end
    end
    return out
end