
# Elimination data
struct KelData{T1<:Number,T2<:Number,T3<:Number,T4<:Number,T5<:Number,T6<:Number}
    s::Array{T1, 1}
    e::Array{T2, 1}
    a::Array{T3, 1}
    b::Array{T4, 1}
    r::Array{T5, 1}
    ar::Array{T6, 1}
    function KelData(s::Vector{T1}, e::Vector{T2}, a::Vector{T3}, b::Vector{T4}, r::Vector{T5}, ar::Vector{T6})::KelData where T1 <: Number where T2 <: Number where T3 <: Number where T4 <: Number where T5 <: Number where T6 <: Number
        new{T1,T2,T3,T4,T5,T6}(s, e, a, b, r, ar)::KelData
    end
    function KelData()::KelData
        KelData(Float64[], Float64[], Float64[], Float64[], Float64[], Float64[])
    end
end

function Base.push!(keldata::KelData, s, e, a, b, r, ar)
    push!(keldata.s, s)
    push!(keldata.e, e)
    push!(keldata.a, a)
    push!(keldata.b, b)
    push!(keldata.r, r)
    push!(keldata.ar, ar)
end

function  Base.length(keldata::KelData)
    return length(keldata.a)
end

# Elimination settings
mutable struct ElimRange{Symbol}
    kelstart::Int
    kelend::Int
    kelexcl::Vector{Int}
    function ElimRange(kelstart, kelend, kelexcl)::ElimRange
        if kelstart > kelend throw(ArgumentError("Kel start > kel end")) end
        if kelstart < 0 throw(ArgumentError("Kel start point < 0")) end
        if kelend   < 0 throw(ArgumentError("Kel endpoint < 0")) end
        if any(x -> x < 0, kelexcl) throw(ArgumentError("Exclude point < 0")) end
        if kelstart in kelexcl || kelend in kelexcl throw(ArgumentError("Kel start or kel end in exclusion")) end
        new{:point}(kelstart, kelend, kelexcl)::ElimRange
    end
    function ElimRange(kelstart, kelend)
        ElimRange(kelstart, kelend, Vector{Int}(undef, 0))
    end
    function ElimRange(;kelstart = 0, kelend = 0, kelexcl = Int[])
        ElimRange(kelstart, kelend, kelexcl)
    end
end

# Dose settings
struct DoseTime{D <: Number, T <: Number, TAU <: Number}
    dose::D
    time::T
    tau::TAU
    function DoseTime(dose::D, time::T, tau::TAU) where D where T where TAU
        if time < zero(T) throw(ArgumentError("Dose time can't be less zero!")) end
        new{D, T, TAU}(dose, time, tau)::DoseTime
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
    function PKSubject(time::Vector{T}, conc::Vector{O}, kelauto::Bool, kelrange::ElimRange, dosetime::DoseTime, keldata::KelData, sort = Dict()) where T <: Number where O <: Number
        new{T, O}(time, conc, kelauto, kelrange, dosetime, keldata, sort)::PKSubject
    end
    function PKSubject(time::Vector, conc::Vector, kelauto::Bool, kelrange::ElimRange, dosetime::DoseTime, sort::Dict)
        PKSubject(time, conc, kelauto, kelrange, dosetime, KelData(), sort)
    end
    function PKSubject(time::Vector, conc::Vector, sort::Dict)
        PKSubject(time, conc, true, ElimRange(), DoseTime(NaN, 0), KelData(), sort)
    end
    function PKSubject(time::Vector, conc::Vector, kelauto::Bool, kelrange, dosetime; sort = Dict())
        PKSubject(time, conc, kelauto, kelrange, dosetime, KelData(), sort)
    end
    function PKSubject(time::Vector, conc::Vector; sort = Dict())
        PKSubject(time, conc, true, ElimRange(), DoseTime(NaN, 0), KelData(), sort)
    end
end

function Base.length(obj::T) where T <: AbstractSubject
    length(obj.time)
end


struct NCAResult{T} <: AbstractSubjectResult{T}
    subject::T
    method::Symbol
    result::Dict
    id::Dict
    function NCAResult(subject::T, method, result, id) where T <: AbstractSubject
        new{T}(subject, method, result, id)
    end
    function NCAResult(subject::T, method, result) where T <: AbstractSubject
        NCAResult(subject, method, result, Dict())
    end
end

"""

Rule for PK subject.

STEP 1 (NaN step): replace all NaN values with nan reyword value (if nan !== NaN);
STEP 2 (LLOQ step): replace values below lloq with btmax value if this value befor Tmax or with atmax if this value after Tmax (if lloq !== NaN);
STEP 3 (remove NaN): rm == true, then remove all NaN values;
"""
struct LimitRule{T<:Real}
    lloq::T
    btmax::Float64
    atmax::Float64
    nan::Float64
    rm::Bool
    function LimitRule(lloq::T, btmax, atmax, nan, rm::Bool) where T <: Real
        new{T}(lloq, btmax, atmax, nan, rm)::LimitRule
    end
    function LimitRule(;lloq = NaN, btmax = NaN, atmax = NaN, nan = NaN, rm::Bool = false)
        LimitRule(lloq, btmax, atmax, nan, rm)
    end
end
