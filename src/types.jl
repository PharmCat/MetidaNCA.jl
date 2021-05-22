
# Elimination data
struct KelData
    s::Vector{Int}
    e::Vector{Int}
    a::Vector{Float64}
    b::Vector{Float64}
    r::Vector{Float64}
    ar::Vector{Float64}
    function KelData(s, e, a, b, r, ar)::KelData
        new(s, e, a, b, r, ar)::KelData
    end
    function KelData()::KelData
        KelData(Int[], Int[], Float64[], Float64[], Float64[], Float64[])
    end
end

# Elimination settings
mutable struct ElimRange
    kelstart::Int
    kelend::Int
    kelexcl::Vector{Int}
    function ElimRange(kelstart, kelend, kelexcl)::ElimRange
        if kelstart < 0 throw(ArgumentError("Kel start point < 0")) end
        if kelend   < 0 throw(ArgumentError("Kel endpoint < 0")) end
        if any(x -> x < 0, kelexcl) throw(ArgumentError("Exclude point < 0")) end
        new(kelstart, kelend, kelexcl)::ElimRange
    end
    function ElimRange(kelstart, kelend)
        ElimRange(kelstart, kelend, Vector{Int}(undef, 0))
    end
    function ElimRange()
        ElimRange(0, 0, Vector{Int}(undef, 0))
    end
end

# Dose settings
struct DoseTime{D <: Number, T <: Number, TAU <: Number}
    dose::D
    time::T
    tau::TAU
    function DoseTime(dose, time, tau)
        new{typeof(dose), typeof(time), typeof(tau)}(dose, time, tau)::DoseTime
    end
    function DoseTime(;dose = NaN, time = 0, tau = NaN)
        DoseTime(dose, time, tau)
    end
    function DoseTime(dose)
        DoseTime(dose, 0, NaN)
    end
    function DoseTime(dose, time)
        DoseTime(dose, time, NaN)
    end
end

# PK subject
mutable struct PKSubject{T <: Number, O <: Number} <: AbstractSubject
    time::Vector{T}
    obs::Vector{O}
    kelauto::Bool
    kelrange::ElimRange
    dosetime::DoseTime
    keldata::KelData
    id::Dict
    function PKSubject(time::Vector, conc::Vector, kelauto::Bool, kelrange::ElimRange, dosetime::DoseTime, keldata::KelData, sort = Dict())
        new{eltype(time), eltype(conc)}(time, conc, kelauto, kelrange, dosetime, keldata, sort)::PKSubject
    end
    function PKSubject(time::Vector, conc::Vector, kelauto::Bool, kelrange, dosetime; sort = Dict())
        PKSubject(time, conc, kelauto, kelrange, dosetime, KelData(), sort)
    end
    function PKSubject(time::Vector, conc::Vector, sort::Dict)
        PKSubject(time, conc, true, ElimRange(), DoseTime(NaN, 0), KelData(), sort)
    end
    function PKSubject(time::Vector, conc::Vector; sort = Dict())
        PKSubject(time, conc, true, ElimRange(), DoseTime(NaN, 0), KelData(), sort)
    end
end
