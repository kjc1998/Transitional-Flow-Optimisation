using Flows
using ToySystems
using ToySystems.NineModeSystemEq
using Optim
using PyPlot; pygui(true)
include("9msctrl.jl")

function AdjointControl(u0, Re, T1, Tstep, Tend, α = 0)
    
    ts = T1:Tstep:Tend
    timestep = Tstep
    T = Tend
    
    function minus_energy_grad(t, a, dadt, v, dvdt)
        dvdt[1] += -2*(a[1] - 1)
        for i = 2:9
            dvdt[i] += -2*a[i]
        end
        return dvdt
    end
    
    objfun(u::AbstractVector) = (u[1] - 1)^2 + sum(u[i]^2 for i = 2:9)
    function objfun(t, u, dudt, I, dIdt)
        return dIdt[1] = (u[1] - 1)^2 + sum(u[i]^2 for i = 2:9)
    end
    function objfun(xq::Coupled)
        return objfun(xq[1])
    end
    
    function Objective_fun_CTRL(η::AbstractVector)
        f = NineModeSystemCTRL(Re, ts, η)
        ϕ = flow(couple(f,objfun),RK4(couple(zeros(9),zeros(1))), TimeStepConstant(timestep))
        I0 = Float64[0.0]
        u0_local = copy(u0)
        ϕ(couple(u0_local, I0), extrema(ts))
        return I0[1]/(T-T1) + α.*η2int(η, timestep)
    end
    
    function Objective_fun_Gradient!(G::AbstractVector,η::AbstractVector)
        f = NineModeSystemCTRL(Re, ts, η)
        ϕ = flow(f,RK4(zeros(9)), TimeStepConstant(timestep))
        u0_local = copy(u0)
        ϕ(u0_local, (T1, T), reset!(store))
    
        h = NineModeSystemLin(Re, true, minus_energy_grad)
        ψ_adj = flow(h, RK4(zeros(9), ContinuousMode(true)), TimeStepFromStorage(timestep))
        mon = Monitor(zeros(9), q -> -B[5]*q[5])
        ψ_adj(zeros(9), store, (T, T1), reset!(mon))
        G .= samples(mon)[end:-1:1] .+ 2.0.*α.*η #willautomatically update, vector only!
    end
    
    result = optimize(Objective_fun_CTRL, Objective_fun_Gradient!, zeros(length(ts)),LBFGS(),Optim.Options(f_tol = 1e-4,
                             store_trace = false,
                             show_trace = false))
    
    return result.minimizer
end