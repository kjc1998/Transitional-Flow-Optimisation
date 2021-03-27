# this file
using ToySystems.NineModeSystemEq


function Matching_Points(Input::Real,x::AbstractVector, y::AbstractVector)
    Up = searchsortedfirst(x, Input)
    if Input == x[Up]
        Output = y[Up]
    else
        Low = Up - 1
        Output = y[Low] + ((Input - x[Low])/(x[Up] - x[Low]))*(y[Up]-y[Low])
    end
    return Output
end

struct NineModeSystemCTRL
    sys::NineModeSystem
    ts::AbstractVector
    ηs::AbstractVector
    NineModeSystemCTRL(Re::Real, ts::AbstractVector, ηs::AbstractVector) = new(NineModeSystem(Re), ts, ηs)
end

function (eq::NineModeSystemCTRL)(t::Real, u::AbstractVector, dudt::AbstractVector)
    eq.sys(t, u, dudt)
    η = Matching_Points(t,eq.ts, eq.ηs)
    dudt[5] += η
    return dudt
end

function η2int(η::AbstractVector,ts::Real)
    output = 0
    η2 = η.^2
    for i in 1:length(η2)
        if η2[i] != η2[end]
            output = output + 0.5*(η2[i]+η2[i+1])*ts
        else
            nothing
        end
    end
    return output
end
#notebook
# include("9msctrl.jl")

# ts = 0:1:500
# ηs = sin.(ts)

# f = NineModeSystemCTRL(100, ts, ηs)

# u    = rand(9)
# dudt = zeros(9)
# f(0.4, u, dudt)

# F = flow(f, RK4(zeros(9)), TimeStepConstant(0.01))

# F(u, (0, 10))