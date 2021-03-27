using Flows
using ToySystems
using ToySystems.NineModeSystemEq
using Optim
using PyPlot; pygui(true)
include("9msctrl.jl")
include("Initial condition variation.jl")
include("Perturbation time.jl")
import LinearAlgebra: norm

function Conplot(x::AbstractVector,y::AbstractVector,z)
    contourf(x, y, Z)

# Alternatively, you can manually set the levels
# and the norm:
# lev_exp = np.arange(np.floor(np.log10(z.min())-1),
#                    np.ceil(np.log10(z.max())+1))
# levs = np.power(10, lev_exp)
# cs = ax.contourf(X, Y, z, levs, norm=colors.LogNorm())

#     cbar = fig.colorbar(cs)

#     plt.show()
end