using Flows
using ToySystems
using ToySystems.NineModeSystemEq
using Optim
using PyPlot; pygui(true)
include("9msctrl.jl")

# Re = 500
# u0 = rand(9)
# dadt = zeros(9)
# q    = zeros(9)
# dqdt = zeros(9)
# T     = 1000
# timestep = 0.5
# ts = 0:timestep:T
# store = RAMStorage(zeros(9))
# B = [0,0,0,0,1,0,0,0,0]
# α = 0.2

objfun(u::AbstractVector) = (u[1] - 1)^2 + sum(u[i]^2 for i = 2:9)
function objfun(t, u, dudt, I, dIdt)
    return dIdt[1] = (u[1] - 1)^2 + sum(u[i]^2 for i = 2:9)
end
function objfun(xq::Coupled)
    return objfun(xq[1])
end


function Objective_con_WCTRL(Re, n, u0, T, timestep)
    f = NineModeSystem(Re)
    ϕ = flow(couple(f,objfun),RK4(couple(zeros(9),zeros(1))), TimeStepConstant(timestep))
    mon = Monitor(couple(zeros(9), zeros(1)), objfun)
    
    I0 = Float64[0.0]
    Tm = round(T/n)
    t0 = 0
    tm = Tm
    for i in 1:n
        ϕ(couple(u0, I0), (t0,tm),mon)
        tm = tm+Tm
        t0 = tm-Tm+timestep
    end
    plot(0:timestep:T, samples(mon))
    return samples(mon), times(mon)
end