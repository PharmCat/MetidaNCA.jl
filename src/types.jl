
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
* `kelexcl` - excluded points (not time value).
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

abstract type AbstractCovariate end

struct ConstantCovariate{T} <: AbstractCovariate
    val::T
end
struct TimeVaryingCovariate{T} <: AbstractCovariate
    val::T
end
function makecovariate(v::AbstractVector)
    if length(v) == 1
        return ConstantCovariate(first(v))
    elseif all(x-> first(v) == x, v)
        return ConstantCovariate(first(v))
    else
        return TimeVaryingCovariate(v)
    end
end
value(v::AbstractCovariate) = v.val
Base.getindex(v::ConstantCovariate, ::Any) = v.val
Base.getindex(v::TimeVaryingCovariate, i) = v.val[i]
Base.length(::ConstantCovariate) = 1
Base.length(v::TimeVaryingCovariate) = length(value(v))

# Dose settings
"""
    DoseTime(dose::D, time::T, tau::TAU) where D <: Number where T <: Number where TAU <: Number

Dose settings.

* `dose` - dose;
* `time` - dose time;
* `tau` - tau (τ);

Dose time set 0 by default.
"""
struct DoseTime{D <: Number, T <: Number, TAU <: Number, R <: Number, RO}
    dose::D
    time::T
    tau::TAU
    rate::R
    route::RO
    function DoseTime(dose::D, time::T, tau::TAU, rate::R, route::RO) where D <: Number where T <: Number where TAU <: Number where R <: Number where RO
        if time < zero(T) throw(ArgumentError("Dose time can't be less zero!")) end
        new{D, T, TAU, R, RO}(dose, time, tau, rate, route)::DoseTime
    end
    function DoseTime(dose, time, tau)
        if isa(time, Int) time_ = float(time) else time_ = time end
        if isa(dose, Int) dose_ = float(dose) else dose_ = dose end
        if isa(tau, Int) tau_   = float(tau)  else tau_  = tau  end
        DoseTime(dose, time, tau, 0.0, nothing)
    end
    function DoseTime(;dose = NaN, time = 0.0, tau = NaN, rate = 0.0, route = nothing)
        if isa(time, Int) time_ = float(time) else time_ = time end
        if isa(dose, Int) dose_ = float(dose) else dose_ = dose end
        if isa(tau, Int) tau_   = float(tau)  else tau_  = tau  end
        DoseTime(dose, time, tau, rate, route)
    end
end

function Base.convert(::Type{<: DoseTime{D,T,TAU,R,RO}}, x::DT) where D <: Number where T <: Number where TAU <: Number where R <: Number where RO where DT <: DoseTime
    DoseTime(one(D)*x.dose, one(T)*x.time, one(TAU)*x.tau, one(R)*x.rate, x.route)
end

function Base.convert(::Type{<: Union{DoseTime, Vector{DoseTime}}}, x::Vector{DT}) where DT <: DoseTime
    return Vector{DoseTime}(x)
end

#Base.first(x::DoseTime) = x

#function checkdosetime(dt::DoseTime)
#    true
#end
#function checkdosetime(::Nothing)
#    true
#end

function checkdosetime(dt::Vector{<:DoseTime})
    if length(dt) == 0 
        error("DoseTime can't be empty.") 
    elseif length(dt) == 1 
        return true
    else
        for i = 2:length(dt)
            if dt[i].time < dt[i-1].time return false end
        end
        return true
    end
end

# PK subject
"""
    PKSubject(time::Vector{T}, conc::Vector{O}, kelauto::Bool, kelrange::ElimRange, dosetime::DoseTime, keldata::KelData, id::Dict{Symbol, V} = Dict{Symbol, Any}()) where T <: Number where O <: Union{Number, Missing} where V

Pharmacokinetic subject.

Fields:

* time::Vector{T} - time values;
* obs::Vector{O} - observations;
* kelauto::Bool 
* kelrange::ElimRange
* dosetime::DoseTime
* keldata::KelData
* id::Dict{Symbol, V}
* ncaresobs::Symbol

"""
mutable struct PKSubject{T <: Number, O, C <: Any, V <: Any} <: AbstractSubject
    time::Vector{T}
    obs::O
    covars::C
    kelauto::Bool
    kelrange::ElimRange
    dosetime::Vector{DoseTime}
    keldata::KelData
    id::Dict{Symbol, V}
    ncaresobs::Symbol
    function PKSubject(time::Vector{T}, conc::O, covars::C, kelauto::Bool, kelrange::ElimRange, dosetime, keldata::KelData, id::Dict{Symbol, V} = Dict{Symbol, Any}())  where T <: Number where O where C  where V
        if !checkdosetime(dosetime) error("DoseTime Vector should be sorted.") end
        new{T, O, C, V}(time, conc, covars, kelauto, kelrange, dosetime, keldata, id, NCARESOBS)::PKSubject
    end
    function PKSubject(time::Vector{T}, conc, kelauto::Bool, kelrange::ElimRange, dosetime, id)  where T 
        PKSubject(time, conc, nothing, kelauto, kelrange, dosetime, KelData(T[], T[], Float64[], Float64[], Float64[], Float64[], Int[]), id)
    end
    function PKSubject(time::Vector{T}, conc, covars, kelauto::Bool, kelrange::ElimRange, dosetime, id)  where T 
        PKSubject(time, conc, covars, kelauto, kelrange, dosetime, KelData(T[], T[], Float64[], Float64[], Float64[], Float64[], Int[]), id)
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


struct NCAOptions{P <: Union{AbstractVector, Nothing}}
    adm::Symbol
    calcm::Symbol
    intpm::Symbol
    partials::P
    prtext::Symbol
    verbose::Int
    warn::Bool
    io::IO
    modify!::Function
    function NCAOptions(adm, calcm, intpm, partials::P, prtext, verbose, warn, io, modify!) where P
        new{P}(adm, calcm, intpm, partials, prtext, verbose, warn, io, modify!)
    end
    function NCAOptions(;adm = :ev, calcm = :lint, intpm = :calcm, partials = nothing, prtext = :err, verbose  = 0, warn = true, io = stdout, modify! = identity)
        NCAOptions(adm, calcm, intpm, partials, prtext, verbose, warn, io, modify!) 
    end
end

"""
    NCAResult(subject::T, options, result::Dict{Symbol, U}) where T <: AbstractSubject where U

NCA resulst.

Fields:

* data::T
* options::Dict{Symbol}
* result::Dict{Symbol, U}
"""
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

See also: [`applylimitrule!`](@ref)
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
    dosetime::DoseTime
    id::Dict{Symbol, V}
    function PDSubject(time::Vector{T}, conc::Vector{O}, bl, th, dosetime::DoseTime, id::Dict{Symbol, V} = Dict{Symbol, Any}()) where T <: Number where O <: Union{Number, Missing} where V
        new{T, O, V}(time, conc, bl, th, dosetime, id)::PDSubject
    end
    function PDSubject(time::Vector, conc::Vector, bl, th, id::Dict{Symbol, V}) where V
        PDSubject(time, conc, bl, th, DoseTime(), id)
    end
end

function gettime(subj::T) where T <: AbstractSubject
    getfield(subj, :time)
end
function getobs_(obsvals::AbstractVector, ::Nothing)
    obsvals
end 
function getobs_(obsvals, ::Nothing)
    first(obsvals)
end 
function getobs_(obsvals, obs::Symbol)
    obsvals[obs]
end 
function getobs(subj::T, obs::Union{Symbol, Nothing} = nothing) where T <: AbstractSubject
    getobs_(getfield(subj, :obs), obs)
end

function getfirstkey(subj::PKSubject)
    getfirstkey_(subj.obs)
end
function getfirstkey_(obsvals::AbstractVector)
    nothing
end
function getfirstkey_(obsvals)
    first(keys(obsvals))
end


function ismultobs(subj::PKSubject)
   ismultobs_(subj.obs)
end
function ismultobs_(obsvals::AbstractVector)
    false
end
function ismultobs_(obsvals)
    true
end

function obsnumber(subj::PKSubject)
   obsnumber_(subj.obs)
end
function obsnumber_(obsvals::AbstractVector)
    1
end
function obsnumber_(obsvals)
    length(obsvals)
end



struct NCAUnits{T, O, D, V}
    time::T
    obs::O
    dose::D
    vol::V
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
