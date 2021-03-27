objfun(u::AbstractVector) = (u[1] - 1)^2 + sum(u[i]^2 for i = 2:9)
function objfun(t, u, dudt, I, dIdt)
    return dIdt[1] = (u[1] - 1)^2 + sum(u[i]^2 for i = 2:9)
end
function objfun(xq::Coupled)
    return objfun(xq[1])
end