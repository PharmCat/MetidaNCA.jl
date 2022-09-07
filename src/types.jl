
# Elimination data
struct KelData{S<:Number, E<:Number}
    s::Vector{S}
    e::Vector{E}
    a::Vector{Float64}
    b::Vector{Float64}
    r::Vector{Float64}
    ar::Vector{Float64}
    n::Vector{Int}
    function KelData(s::Vector{S}, e::Vector{E}, a::Vector{A}, b::Vector{B}, r, ar, n)::KelData where S <: Number where E <: Number where A where B
        new{S, E}(s, e, a, b, r, ar, n)::KelData
    end
    function KelData()::KelData
        KelData(Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Int[])
    end
end

function resize!(keldata::KelData)
    resize!(keldata, 0)
end
function resize!(keldata::KelData, i::Int)
    resize!(keldata.s, i)
    resize!(keldata.e, i)
    resize!(keldata.a, i)
    resize!(keldata.b, i)
    resize!(keldata.r, i)
    resize!(keldata.ar, i)
    resize!(keldata.n, i)
    keldata
end

function Base.push!(keldata::KelData, s, e, a, b, r, ar, n)
    push!(keldata.s, s)
    push!(keldata.e, e)
    push!(keldata.a, a)
    push!(keldata.b, b)
    push!(keldata.r, r)
    push!(keldata.ar, ar)
    push!(keldata.n, n)
end

function  Base.length(keldata::KelData)
    return length(keldata.a)
end

# Elimination settings
"""
    ElimRange(kelstart::Int, kelend::Int, kelexcl::Vector{Int})::ElimRange

Elimination settings for PK subject.

* `kelstart` - start point;
* `kelend` - end point;
* `kelexcl` - excluded points.
"""
mutable struct ElimRange{Symbol}
    kelstart::Int
    kelend::Int
    kelexcl::Vector{Int}
    function ElimRange(kelstart::Int, kelend::Int, kelexcl::Vector{Int}; time = false)::ElimRange
        if kelstart > kelend throw(ArgumentError("Kel start > kel end")) end
        if kelstart < 0 throw(ArgumentError("Kel start point < 0")) end
        if kelend   < 0 throw(ArgumentError("Kel endpoint < 0")) end
        if any(x -> x < 0, kelexcl) throw(ArgumentError("Exclude point < 0")) end
        if kelstart in kelexcl || kelend in kelexcl throw(ArgumentError("Kel start or kel end in exclusion")) end
        new{:point}(kelstart, kelend, kelexcl)::ElimRange
    end
    function ElimRange(;kelstart = 0, kelend = 0, kelexcl = Int[])
        ElimRange(kelstart, kelend, kelexcl)
    end
end

# Dose settings
"""
    DoseTime(dose::D, time::T, tau::TAU) where D <: Number where T <: Number where TAU <: Number

Dose settings.

* `dose` - dose;
* `time` - dose time;
* `tau` - tau (τ);

Dose time set 0 by default.
"""
struct DoseTime{D <: Number, T <: Number, TAU <: Number}
    dose::D
    time::T
    tau::TAU
    function DoseTime(dose::D, time::T, tau::TAU) where D <: Number where T <: Number where TAU <: Number
        if time < zero(T) throw(ArgumentError("Dose time can't be less zero!")) end
        new{D, T, TAU}(dose, time, tau)::DoseTime
    end
    function DoseTime(;dose = NaN, time = 0.0, tau = NaN)
        DoseTime(dose, time, tau)
    end
    #=
    function DoseTime(dose)
        DoseTime(dose, 0, NaN)
    end
    function DoseTime(dose, time)
        DoseTime(dose, time, NaN)
    end
    =#
end

# PK subject
mutable struct PKSubject{T <: Number, O <: Union{Number, Missing}, V <: Any} <: AbstractSubject
    time::Vector{T}
    obs::Vector{O}
    kelauto::Bool
    kelrange::ElimRange
    dosetime::DoseTime
    keldata::KelData
    id::Dict{Symbol, V}
    function PKSubject(time::Vector{T}, conc::Vector{O}, kelauto::Bool, kelrange::ElimRange, dosetime::DoseTime, keldata::KelData, sort::Dict{Symbol, V} = Dict{Symbol, Any}()) where T <: Number where O <: Union{Number, Missing} where V
        new{T, O, V}(time, conc, kelauto, kelrange, dosetime, keldata, sort)::PKSubject
    end
    function PKSubject(time::Vector{T}, conc::Vector{O}, kelauto::Bool, kelrange::ElimRange, dosetime::DoseTime, sort::Dict{Symbol, V}) where T where O where V
        PKSubject(time, conc, kelauto, kelrange, dosetime, KelData(T[], T[], Float64[], Float64[], Float64[], Float64[], Int[]), sort)
    end
    #=
    function PKSubject(time::Vector, conc::Vector, sort::Dict)
        PKSubject(time, conc, true, ElimRange(), DoseTime(NaN, 0), KelData(), sort)
    end
    function PKSubject(time::Vector, conc::Vector, kelauto::Bool, kelrange, dosetime; sort = Dict())
        PKSubject(time, conc, kelauto, kelrange, dosetime, KelData(), sort)
    end
    function PKSubject(time::Vector, conc::Vector; sort = Dict())
        PKSubject(time, conc, true, ElimRange(), DoseTime(NaN, 0), KelData(), sort)
    end
    =#
end

function Base.length(obj::T) where T <: AbstractSubject
    length(obj.time)
end


struct NCAResult{T, U} <: AbstractSubjectResult{T}
    data::T
    options::Dict{Symbol}
    result::Dict{Symbol, U}
    function NCAResult(subject::T, options, result::Dict{Symbol, U}) where T <: AbstractSubject where U
        new{T, U}(subject, options, result)
    end
    #=
    function NCAResult(subject::T, method, result) where T <: AbstractSubject
        NCAResult(subject, method, result, Dict())
    end
    =#
end

"""
    LimitRule(lloq::T, btmax, atmax, nan, rm::Bool) where T <: Real

    LimitRule(;lloq = NaN, btmax = NaN, atmax = NaN, nan = NaN, rm::Bool = false)

* `lloq` - LLOQ - low limit of quantification;
* `btmax` - value for points before Tmax;
* `atmat` - values for points after Tmax;
* `nan` - values for replacing `NaN`;
* `rm` - if `true`, removee all `NaN` points.

Rule for PK subject.

* STEP 1 (NaN step): replace all `NaN` and `missing` values with nan keyword value (if `nan` not NaN);
* STEP 2 (LLOQ step): replace values below `lloq` with `btmax` value if this value befor Tmax or with atmax if this value after Tmax (if `lloq` not NaN);
* STEP 3 (remove NaN): `rm` == true, then remove all `NaN` and `missing` values.
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
    function LimitRule(;nan = NaN, lloq = NaN, btmax = NaN, atmax = NaN,  rm::Bool = false)
        LimitRule(lloq, btmax, atmax, nan, rm)
    end
end
#=
function isapplicable(lr::LimitRule)
    !isnan(lr.lloq) || !isnan(lr.nan) || lr.rm ? true : false
end
=#
#Urine PK subject
mutable struct UPKSubject{T <: Tuple{Number, Number}, O <: Union{Number, Missing}, VOL <: Union{Number, Missing}, V <: Any} <: AbstractSubject
    time::Vector{T}
    obs::Vector{O}
    vol::Vector{VOL}
    kelauto::Bool
    kelrange::ElimRange
    dosetime::DoseTime
    keldata::KelData
    id::Dict{Symbol, V}
    function UPKSubject(time::Vector{T}, conc::Vector{O}, vol::Vector{VOL}, kelauto::Bool, kelrange::ElimRange, dosetime::DoseTime, keldata::KelData, id::Dict{Symbol, V} = Dict{Symbol, Any}()) where T <: Tuple{Number, Number} where O <: Union{Number, Missing} where VOL <: Union{Number, Missing} where V
        new{T, O, VOL, V}(time, conc, vol, kelauto, kelrange, dosetime, keldata, id)
    end
    function UPKSubject(time::AbstractVector{Tuple{S,E}}, conc::Vector, vol::Vector, kelauto::Bool, kelrange::ElimRange, dosetime::DoseTime, id::Dict{Symbol, V}) where V where S where E
        ttype = promote_type(S, E)

        UPKSubject(time, conc, vol, kelauto, kelrange, dosetime, KelData(ttype[], ttype[], Float64[], Float64[], Float64[], Float64[], Int[]), id)
    end
end

# PD subject
mutable struct PDSubject{T <: Number, O <: Union{Number, Missing}, V <: Any} <: AbstractSubject
    time::Vector{T}
    obs::Vector{O}
    bl::Float64
    th::Float64
    id::Dict{Symbol, V}
    function PDSubject(time::Vector{T}, conc::Vector{O}, bl, th, sort::Dict{Symbol, V} = Dict{Symbol, Any}()) where T <: Number where O <: Union{Number, Missing} where V
        new{T, O, V}(time, conc, bl, th, sort)::PDSubject
    end
    function PDSubject(time::Vector, conc::Vector, bl, th, sort::Dict{Symbol, V}) where V
        PDSubject(time, conc, bl, th, sort)
    end
end

function gettime(subj::T) where T <: AbstractSubject
    getfield(subj, :time)
end
function getobs(subj::T) where T <: AbstractSubject
    getfield(subj, :obs)
end


struct NCAOptions
    adm::Symbol
    calcm::Symbol
    intpm::Symbol
    verbose::Int
    warn::Bool
    io::IO
    modify!::Function
end

#=
function gettime(subj::PKSubject)
end
function getobs(subj::PKSubject)
end
function gettime(subj::UPKSubject)
end
function getobs(subj::UPKSubject)
end
=#
