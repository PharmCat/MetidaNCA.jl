# Fills dictionary d with indexes of individual values

function nonunique(v)
    cm = StatsBase.countmap(v)
    return [x for x in unique(v) if cm[x] > 1]
end



parse_gkw(s::String) = [Symbol(s)]
parse_gkw(s::Symbol) = [s]
parse_gkw(s::AbstractVector{<:AbstractString}) = Symbol.(s)
parse_gkw(s::AbstractVector{Symbol}) = s


checkblq(blq::AbstractFloat)  = isnan(blq) ? NaN : error("Wrong BLQ value: $blq")
checkblq(blq::AbstractString) = blq == "NaN" ? NaN : blq
checkblq(::Nothing)           = nothing
checkblq(blq)                 = error("Wrong BLQ value: $blq")

obsvalue(x::Number, ::AbstractFloat)             = float(x)
obsvalue(x::Number, ::Any)                       = (y = float(x); isnan(y) ? missing : y)
obsvalue(s::AbstractString, blq::AbstractString) = s == blq ? NaN : obsvalue(tryparse(Float64, s), blq)
obsvalue(s::AbstractString, blq)                 = obsvalue(tryparse(Float64, s), blq)
obsvalue(::Any, ::Any)                           = missing

# Returns a new vector, `vec` is not changed.
function checkobsvec(vec::AbstractVector, warn, blq)
    b = checkblq(blq)
    T = nonmissingtype(eltype(vec))
    if warn
        T !== Union{} && T <: Integer && @warn "Concentration values transformed to float."
        T <: Number || @warn "Some concentration values maybe not a number, try to fix."
    end
    out = [obsvalue(x, b) for x in vec]
    if warn && !(T <: Number)
        for (x, y) in zip(vec, out)
            ismissing(y) && !ismissing(x) && !(x isa Number) && @warn "Value $x parsed as missing"
        end
    end
    return out
end


function checktimevec(vec::AbstractVector)
    timevals = identity.(vec)
    eltype(timevals) <: Number || error("Some time values not a number ($(eltype(timevals))))!")
    return timevals
end

function checkvalues!(timevals_sp, obsvals; warn = true, blq = NaN)
    timevals = checktimevec(timevals_sp)
    
    obsvals = checkobsvec.(obsvals, warn, blq)
    
    timevals, obsvals
end

function makedosetimevec(::Nothing, v)
    [DoseTime(NaN, v, NaN)] 
end
function makedosetimevec(x::DoseTime, ::Any)
    [x] 
end
function makedosetimevec(x::AbstractVector{<: DoseTime}, ::Any) 
    x
end

"""
    pkimport(data, time, obs, sort;
        elimrange = ElimRange(),
        dosetime = DoseTime(),
        limitrule::Union{Nothing, LimitRule} = nothing,
        warn = true,
        kwargs...)

Import PK data from table `data`.

* `time` - time column;
* `obs` - concentration column;
* `sort` - subject sorting columns.

keywords:
* `blq` - BLQ value, can be only NaN or string values, by default  blq = NaN;
* `elimrange` - set elimination range settings;
* `dosetime` - set dose and dose time, by default dosetime = 0, dose is `NaN`;
* `limitrule` - apply limitrule to subject;
* `warn` - false for warnings supress;
* checktime (true) - check uniqueness of time values; 
* nutcfunk (last) - function used to chose concentration value for identical time (if checktime == `true`).

!!! note

    BLQ replacement:
    1. `NaN` - internal representation of BLQ value
    2. If `blq` is string and values is string - replace "NaN" -> "missing" then "blq" -> NaN
    2. If `blq` is string and values is not string - replace "NaN" -> "missing" then "blq" -> NaN
    


See also: [`ElimRange`](@ref), [`DoseTime`](@ref), [`LimitRule`](@ref).
"""
function pkimport(data, time, obs, sort;
    blq = NaN,
    covars = [],
    elimrange = ElimRange(), 
    dosetime = nothing, 
    warn = true, 
    nutcfunk = last,
    units::Union{Nothing, NCAUnits} = nothing)

    checkblq(blq)
    
    sort = parse_gkw(sort)
    obs  = parse_gkw(obs)

    if length(covars) > 0
        covars = parse_gkw(covars)
    end

    Tables.istable(data) || error("Data not a table!")

    cols   = Tables.columns(data)
    cdata  = Tuple(Tables.getcolumn(cols, y) for y in sort)
    d      = Dict{Tuple{eltype.(cdata)...}, Vector{Int}}()
    indsdict!(d, cdata)

    timec = Tables.getcolumn(data, time)
    #concc = Tables.getcolumn(data, obs)

    dosetime_ = makedosetimevec(dosetime, zero(eltype(timec)))

    if !checkdosetime(dosetime_)
        warn && @warn "DoseTime sorted..."
        sort!(dosetime_, by = x -> x.time)
    end

    any(isnanormissing, timec) && error("Some time values is NaN or Missing!")

    sdata = Vector{PKSubject}(undef, length(d))
    i = one(Int)
    korder = sort!(collect(keys(d)), by = k -> first(d[k]))
    @inbounds for k in korder
        v = d[k]
        timevals = view(timec, v)

        if  !allunique(timevals)
            nuv      = nonunique(timevals) # non unique values
            error("Subject: $k, Not all time values is unique ($nuv)!")
        end

        permute!(v, sortperm(timevals))
        timevals_spv = view(timec, v)
        obsvals_spv  = [view(Tables.getcolumn(data, y), v) for y in obs]
        

        timevals_sp, obsvals_sp = checkvalues!(timevals_spv, obsvals_spv; warn = warn, blq = blq)

        if length(covars) > 0
            covars_v = (; zip(covars, [makecovariate(Tables.getcolumn(data, y)[v]) for y in covars])...)
        else
            covars_v = nothing
        end
 
        sdata[i] = PKSubject(timevals_sp, (; zip(obs, obsvals_sp)...), covars_v, elimrange,  dosetime_, Dict(sort .=> k))
        i += one(Int)
    end
    ds = DataSet(identity.(sdata))

    if !(isnothing(units)) ds.metadata[:units] = units end
    return ds
end
"""
    pkimport(data, time, obs;
        warn = true,
        kwargs...)

Import PK data from tabular data `data`, `time` - time column, `obs` - concentration column.
"""
function pkimport(data, time, obs; 
    blq = NaN, elimrange = ElimRange(), dosetime = nothing, id = Dict{Symbol, Any}(), warn = true)
    
    checkblq(blq)
    
    obs   = parse_gkw(obs)

    obsvals = [copy(Tables.getcolumn(data, y)) for y in obs]

    timevals, obsvals = checkvalues!(Tables.getcolumn(data, time), obsvals; warn = warn, blq = blq)
    
    dosetime_ = makedosetimevec(dosetime, zero(eltype(timevals)))

    if length(timevals) > 1
        for i = 1:length(timevals)-1
            if timevals[i] > timevals[i+1]
                spt = sortperm(timevals)
                for o in obsvals
                    permute!(o, spt)
                end
                permute!(timevals, spt)
                break
            end
        end
    end
    if !checkdosetime(dosetime_)
        warn && @warn "DoseTime sorted..."
        sort!(dosetime_, by = x -> x.time)
    end

    pks = PKSubject(timevals, (; zip(obs, obsvals)...), elimrange,  dosetime_, id)

    pks
end
"""
    pkimport(time, obs;
        elimrange = ElimRange(),
        dosetime = DoseTime(),
        id = Dict{Symbol, Any}(),
        limitrule::Union{Nothing, LimitRule} = nothing,
        warn = true,
        kwargs...)

Import PK data from time vector `time` and concentration vector `obs`.

"""
function pkimport(time, obs; 
    blq = NaN, obsname = NCARESOBS, elimrange = ElimRange(), dosetime = nothing, id = Dict{Symbol, Any}(), warn = true)
    
    checkblq(blq)

    timevals = checktimevec(time)
   
    obsvals  = checkobsvec(copy(obs), warn, blq)

    dosetime_ = makedosetimevec(dosetime, zero(eltype(timevals)))

    if !checkdosetime(dosetime_)
        warn && @warn "DoseTime sorted..."
        sort!(dosetime_, by = x -> x.time)
    end

    if length(timevals) > 1
        for i = 1:length(timevals)-1
            if timevals[i] > timevals[i+1]
                spt = sortperm(timevals)
                permute!(obsvals, spt)
                permute!(timevals, spt)
                break
            end
        end
    end
    pks = PKSubject(timevals, NamedTuple{(obsname,)}((obsvals,)), elimrange, dosetime_, id)

    pks
end

function pkimport(data; time, conc, sort = nothing, kwargs...)
    if isnothing(sort)
        return pkimport(data, time, conc; kwargs...)
    end
    return pkimport(data, time, conc, sort; kwargs...)
end

"""
    upkimport(data, stime, etime, conc, vol, sort; elimrange = ElimRange(), dosetime = DoseTime())

Import urine PK data from table `data`.

* `stime` - start time column;
* `etime` - end time column;
* `conc` - concentration column;
* `vol` - volume column;
* `sort` - subject sorting columns.
"""
function upkimport(data, stime, etime, conc, vol, sort; kwargs...)
    
    sort = parse_gkw(sort)

    cols   = Tables.columns(data)
    cdata  = Tuple(Tables.getcolumn(cols, y) for y in sort)
    d      = Dict{Tuple{eltype.(cdata)...}, Vector{Int}}()
    indsdict!(d, cdata)

    tnames = Symbol.(Tables.columns(data)) 
    isa(stime, Symbol) || stime in tnames || error("column Start Time ($stime) not found")
    isa(etime, Symbol) || etime in tnames || error("column End Time ($etime) not found")
    isa(conc, Symbol) || conc in tnames || error("column Concentration ($conc) not found")
    isa(vol, Symbol) || vol in tnames || error("column Volume ($vol) not found")

    stimec = Tables.getcolumn(data, stime)
    etimec = Tables.getcolumn(data, etime)
    concc = Tables.getcolumn(data, conc)
    volc = Tables.getcolumn(data, vol)

    any(isnanormissing, stimec) && error("Some Start Time values is NaN or Missing!")
    any(isnanormissing, etimec) && error("Some End Time values is NaN or Missing!")

    sdata = Vector{UPKSubject}(undef, length(d))
    i = one(Int)
    korder = sort!(collect(keys(d)), by = k -> first(d[k]))
    @inbounds for k in korder
        v = d[k]
        stimevals = view(stimec, v)
        etimevals = view(etimec, v)
        concvals  = view(concc, v)
        volvals   = view(volc, v)

        sdata[i] = upkimport(stimevals, etimevals, concvals, volvals; id = Dict(sort .=> k), kwargs...)
        i += one(Int)
    end
    return DataSet(identity.(sdata))
end
"""
    upkimport(data, stime, etime, conc, vol;  elimrange = ElimRange(), dosetime = DoseTime())

Import single urine PK data from table `data`.

* `stime` - start time column;
* `etime` - end time column;
* `conc` - concentration column;
* `vol` - volume column.
"""
function upkimport(data, stime, etime, conc, vol; kwargs...)
    upkimport(Tables.getcolumn(data, stime), Tables.getcolumn(data, etime), Tables.getcolumn(data, conc), Tables.getcolumn(data, vol); kwargs...)
end
"""
    upkimport(stime, etime, conc, vol; elimrange = ElimRange(), dosetime = DoseTime())

Import urine PK data from time vectors:
* `stime` - start times;
* `etime` - end times;
* `conc` - concentrations;
* `vol` - volumes.
"""
function upkimport(stime, etime, conc, vol;
    blq = NaN, 
    elimrange = ElimRange(), 
    dosetime = nothing, 
    id = Dict{Symbol, Any}(), 
    warn = true)

    any(isnanormissing, stime) && error("Some Start Time values is NaN or Missing!")
    any(isnanormissing, etime) && error("Some End Time values is NaN or Missing!")
    timeranges = collect(zip(stime, etime))
    sp = sortperm(stime)
    timevals_sp = timeranges[sp]
    concvals_sp = conc[sp]
    volvals_sp  = vol[sp]

    #if isnothing(dosetime) dosetime = DoseTime(NaN, zero(promote_type(eltype(stimec), eltype(etimec))), NaN) end

    time_type = promote_type(typeof(zero(eltype(stime))), typeof(zero(eltype(etime))))
    zerotime  = zero(time_type)
    if isnothing(dosetime)
        dosetime = DoseTime(NaN, zerotime, NaN*zerotime)
    else
        if !(time_type <: typeof(dosetime.time)) && !(time_type <: Real)
            warn && @warn "Type of dose time can be wrong... try to fix it"
            dosetime = DoseTime(dosetime.dose, dosetime.time*oneunit(time_type), dosetime.tau)
        end
    end


    if length(timevals_sp) > 1
        for c = 2:length(timevals_sp)
            timevals_sp[c][1] == timevals_sp[c-1][2] || error("Start time ($(timevals_sp[c][1])) for observation $c not equal End time ($(timevals_sp[c-1][2])) for observation $(c-1)!")
        end
    end
    UPKSubject(timevals_sp, concvals_sp, volvals_sp, elimrange,  dosetime, id)
end

"""
    pdimport(data, time, obs, sort;
        bl = 0,
        th = 0,
        limitrule::Union{Nothing, LimitRule} = nothing,
        warn = true)

Import pharmackodynamic data from table:

* `data` - data table;
* `time` - observation time;
* `obs` - observation value;
* `sort` - sorting columns.

Keywords:

* `bl` - baseline;
* `th` - threshold;
* `limitrule` - limit rule;
* `warn` - warning supress if `false`.

"""
function pdimport(data, time, obs, sort; 
    blq = nothing, 
    bl = 0, 
    th = 0, 
    dosetime::Union{Nothing, DoseTime} = nothing, 
    warn = true)
    
    sort = parse_gkw(sort)

    Tables.istable(data) || error("Data not a table!")

    cols   = Tables.columns(data)
    cdata  = Tuple(Tables.getcolumn(cols, y) for y in sort)
    d      = Dict{Tuple{eltype.(cdata)...}, Vector{Int}}()
    indsdict!(d, cdata)

    timec = Tables.getcolumn(data, time)
    obsc  = Tables.getcolumn(data, obs)

    if isnothing(dosetime) dosetime = DoseTime(NaN, zero(eltype(timec)), NaN) end

    any(isnanormissing, timec) && error("Some time values is NaN or Missing!")

    sdata = Vector{PDSubject}(undef, length(d))
    i = one(Int)
    korder = sort!(collect(keys(d)), by = k -> first(d[k]))
    @inbounds for k in korder
        v = d[k]
        timevals = timec[v]
        obsvals  = obsc[v]

        if  !allunique(timevals)
            nuv      = nonunique(timevals) # non unique values
            error("Subject: $k, Not all time values is unique ($nuv)!")
        end
        sp = sortperm(timevals)
        timevals_spv = checktimevec(view(timevals, sp))
        obsvals_spv  = checkobsvec(obsvals[sp], warn, blq)
        sdata[i]     = PDSubject(timevals_spv, obsvals_spv, bl, th, dosetime, Dict(sort .=> k))
        i += one(Int)
    end
    ds = DataSet(identity.(sdata))

    ds
end

"""
    pdimport(data, time, obs;
        warn = true,
        kwargs...)

Import PD data from tabular data `data`, `time` - time column, `obs` - observations column.
"""
function pdimport(data, time, obs; blq = nothing, bl = 0, th = 0, id = Dict{Symbol, Any}(), warn = true)

    timevals = checktimevec(Tables.getcolumn(data, time))
    obsvals  = checkobsvec(copy(Tables.getcolumn(data, obs)), warn, blq)

    spt = sortperm(timevals)
    permute!(obsvals, spt)
    permute!(timevals, spt)

    pds = PDSubject(timevals, obsvals, bl, th,  id)

    pds
end
"""
    pdimport(time, obs;
        bl = 0,
        th = 0,
        id = Dict{Symbol, Any}(),
        warn = true)

Import PD data from time vector `time` and observations vector `obs`.
"""
function pdimport(time, obs; blq = nothing, bl = 0, th = 0, id = Dict{Symbol, Any}(), warn = true)
    timevals = checktimevec(time)
    obsvals  = checkobsvec(copy(obs), warn, blq)

    spt = sortperm(timevals)
    permute!(obsvals, spt)
    permute!(timevals, spt)

    pds = PDSubject(timevals, obsvals, bl, th,  id)

    pds
end
