using Flows
using ToySystems
using ToySystems.NineModeSystemEq
using Optim
import LinearAlgebra: norm



# T     = 20000 #or 1000
# timestepdifference = 0.5



objective(u::AbstractVector) = (u[1]-1)^2 + sum(u[i]^2 for i = 2:9)
function objective(t, u, dudt, I, dIdt)
    return dIdt[1] = (u[1]-1)^2 + sum(u[i]^2 for i = 2:9)
end
function objective(xq::Coupled)
    return objective(xq[1])
end

function P(u0,Re, T, timestepdifference)
    u0_local = copy(u0)
    value = 0
    f = NineModeSystem(Re)
    ϕ = flow(couple(f, objective), RK4(couple(zeros(9), zeros(1))), TimeStepConstant(timestepdifference))
    I0 = Float64[0.0]
    mon = Monitor(couple(zeros(9), zeros(1)), objective)
    ϕ(couple(u0_local, I0), (0, T), reset!(mon))
    a1 = times(mon)
    a2 = samples(mon)
    for i in 1:length(a2)
        if sqrt(a2[i])>0.001 ##if r = 0.001, take the sqrt
            value = length(a2)
        else
            value = i
            break
        end
    end
#     println(a1[i])
#     println(i)
#     plot(a1, a2)
#     return samples(mon)[i]
    return a1[value]
end