# Заполняет словарь d индексами индивидуальных значений

nonunique(v) = [k for (k, v) in StatsBase.countmap(v) if v > 1]

function floatparse(data, warn)
    if !isa(data, AbstractFloat) && !ismissing(data)
        if isa(data, AbstractString)
            tp =  tryparse(Float64, data)
            if isnothing(tp)
                warn && @warn "Value $data parsed as `NaN`"
                return NaN
            else
                return tp
            end
        else
            try
                return float(data)
            catch
                return NaN
            end
        end
    elseif ismissing(data)
        return NaN
    else
        return data
    end
end


#=
function floatparse(data::AbstractVector)
    v = Vector{Float64}(undef, length(data))
    @inbounds for i = 1:length(data)
        if !isa(data[i], AbstractFloat) && !ismissing(data[i])
            if isa(data[i], AbstractString)
                try
                    v[i] = parse(Float64, data[i])
                catch
                    v[i] = NaN
                end
            else
                try
                    v[i] = float(data[i])
                catch
                    v[i] = NaN
                end
            end
        elseif ismissing(data[i])
            v[i] = NaN
        else
            v[i] = data[i]
        end
    end
    identity.(v)
end
=#
#=
Check element type of time column
Check element type of observation / concentration column
return new arrays
=#
function checkvalues(timevals_sp, concvals_sp; warn = true)

    timevals_sp_ = identity.(timevals_sp)
    eltype(timevals_sp_) <: Number || error("Some time values not a number!")

    if !(eltype(concvals_sp) <: Union{Number, Missing})
        warn && @warn "Some concentration values maybe not a number, try to fix."
        concvals_sp_ = floatparse.(concvals_sp, warn)
    else
        concvals_sp_ = identity.(concvals_sp)
    end
    timevals_sp_, concvals_sp_
end

"""
    pkimport(data, time, conc, sort;
        kelauto = true,
        elimrange = ElimRange(),
        dosetime = DoseTime(),
        limitrule::Union{Nothing, LimitRule} = nothing,
        warn = true,
        kwargs...)

Import PK data from table `data`.

* `time` - time column;
* `conc` - concentration column;
* `sort` - subject sorting columns.

keywords:

* `kelauto` - if `true` auto range settings, if `false` used `kelstart`/`kelend` from `elimrange`;
* `elimrange` - set elimination range settings;
* `dosetime` - set dose and dose time, by default dosetime = 0, dose is `NaN`;
* `limitrule` - apply limitrule to subject;
* `warn` - false for warnings supress.

!!! note

    If time column have non-unique values - last pair time-concentration will be used.

"""
function pkimport(data, time, conc, sort; kelauto = true,  elimrange = ElimRange(), dosetime = DoseTime(), limitrule::Union{Nothing, LimitRule} = nothing, warn = true, kwargs...)
    if isa(sort, String) sort = [Symbol(sort)] end
    if isa(sort, Symbol) sort = [sort] end

    Tables.istable(data) || error("Data not a table!")

    cols   = Tables.columns(data)
    cdata  = Tuple(Tables.getcolumn(cols, y) for y in sort)
    d      = Dict{Tuple{eltype.(cdata)...}, Vector{Int}}()
    indsdict!(d, cdata)

    timec = Tables.getcolumn(data, time)
    concc = Tables.getcolumn(data, conc)

    any(isnanormissing, timec) && error("Some time values is NaN or Missing!")

    sdata = Vector{PKSubject}(undef, length(d))
    i = one(Int)
    @inbounds for (k, v) in d
        timevals = view(timec, v)
        concvals = view(concc, v)
        if !allunique(timevals)

            nuv = nonunique(timevals)
            warn && @warn "Not all time values is unique ($nuv), last observation used! ($k)"

            nuvinds = findall(x -> x == first(nuv), timevals)
            resize!(nuvinds, length(nuvinds) - 1)
            if length(nuv) > 1
                for cnt = 2:length(nuv)
                    nuvinds_ = findall(x -> x == nuv[cnt], timevals)
                    resize!(nuvinds_, length(nuvinds_) - 1)
                    append!(nuvinds, nuvinds_)
                end
            end
            sort!(nuvinds)
            deleteat!(v, nuvinds)
            deleteat!(v, nuvinds)
            timevals = view(timec, v)
            concvals = view(concc, v)
        end
        sp = sortperm(timevals)
        timevals_spv = view(timevals, sp)
        concvals_spv = view(concvals, sp)
        timevals_sp, concvals_sp = checkvalues(timevals_spv, concvals_spv; warn = warn)
        sdata[i] = PKSubject(timevals_sp, concvals_sp, kelauto, elimrange,  dosetime, Dict(sort .=> k))
        i += one(Int)
    end
    ds = DataSet(identity.(sdata))
    if !isnothing(limitrule)
        applylimitrule!(ds, limitrule)
    end
    ds
end
"""
    pkimport(data, time, conc;
        warn = true,
        kwargs...)

Import PK data from tabular data `data`, `time` - time column, `conc` - concentration column.
"""
function pkimport(data, time, conc; warn = true, kwargs...)
    timevals_sp, concvals_sp = checkvalues(Tables.getcolumn(data, time), Tables.getcolumn(data, conc); warn = warn)
    pkimport(timevals_sp, concvals_sp; warn = warn, kwargs...)
end
"""
    pkimport(time, conc;
        kelauto = true,
        elimrange = ElimRange(),
        dosetime = DoseTime(),
        id = Dict{Symbol, Any}(),
        limitrule::Union{Nothing, LimitRule} = nothing,
        warn = true,
        kwargs...)

Import PK data from time vector `time` and concentration vector `conc`.
"""
function pkimport(time, conc; kelauto = true,  elimrange = ElimRange(), dosetime = DoseTime(), id = Dict{Symbol, Any}(), limitrule::Union{Nothing, LimitRule} = nothing, warn = true, kwargs...)
    timevals_sp, concvals_sp = checkvalues(time, conc, warn = warn)
    pks = PKSubject(timevals_sp, concvals_sp, kelauto, elimrange,  dosetime, id)
    if !isnothing(limitrule)
        applylimitrule!(pks, limitrule)
    end
    pks
end


"""
    upkimport(data, stime, etime, conc, vol, sort; kelauto = true,  elimrange = ElimRange(), dosetime = DoseTime())

Import urine PK data from table `data`.

* `stime` - start time column;
* `etime` - end time column;
* `conc` - concentration column;
* `vol` - volume column;
* `sort` - subject sorting columns.
"""
function upkimport(data, stime, etime, conc, vol, sort; kelauto = true,  elimrange = ElimRange(), dosetime = DoseTime())
    if isa(sort, String) sort = [Symbol(sort)] end
    if isa(sort, Symbol) sort = [sort] end

    cols   = Tables.columns(data)
    cdata  = Tuple(Tables.getcolumn(cols, y) for y in sort)
    d      = Dict{Tuple{eltype.(cdata)...}, Vector{Int}}()
    indsdict!(d, cdata)

    tnames = Symbol.(names(data))
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
    @inbounds for (k, v) in d
        stimevals = view(stimec, v)
        etimevals = view(etimec, v)
        concvals  = view(concc, v)
        volvals   = view(volc, v)

        sdata[i] = upkimport(stimevals, etimevals, concvals, volvals; kelauto = kelauto,  elimrange = elimrange, dosetime = dosetime, id = Dict(sort .=> k))
        i += one(Int)
    end
    return DataSet(identity.(sdata))
end
"""
    upkimport(data, stime, etime, conc, vol; kelauto = true,  elimrange = ElimRange(), dosetime = DoseTime())

Import single urine PK data from table `data`.

* `stime` - start time column;
* `etime` - end time column;
* `conc` - concentration column;
* `vol` - volume column.
"""
function upkimport(data, stime, etime, conc, vol; kelauto = true,  elimrange = ElimRange(), dosetime = DoseTime())
    upkimport(Tables.getcolumn(data, stime), Tables.getcolumn(data, etime), Tables.getcolumn(data, conc), Tables.getcolumn(data, vol); kelauto = kelauto,  elimrange = elimrange, dosetime = dosetime)
end
"""
    upkimport(stime, etime, conc, vol; kelauto = true,  elimrange = ElimRange(), dosetime = DoseTime())

Import urine PK data from time vectors:
* `stime` - start times;
* `etime` - end times;
* `conc` - concentrations;
* `vol` - volumes.
"""
function upkimport(stime, etime, conc, vol; kelauto = true,  elimrange = ElimRange(), dosetime = DoseTime(), id = Dict{Symbol, Any}())
    any(isnanormissing, stime) && error("Some Start Time values is NaN or Missing!")
    any(isnanormissing, etime) && error("Some End Time values is NaN or Missing!")
    timeranges = collect(zip(stime, etime))
    sp = sortperm(stime)
    timevals_sp = timeranges[sp]
    concvals_sp = conc[sp]
    volvals_sp  = vol[sp]
    if length(timevals_sp) > 1
        for c = 2:length(timevals_sp)
            timevals_sp[c][1] == timevals_sp[c-1][2] || error("Start time ($(timevals_sp[c][1])) for observation $c not equal End time ($(timevals_sp[c-1][2])) for observation $(c-1)!")
        end
    end
    UPKSubject(timevals_sp, concvals_sp, volvals_sp, kelauto, elimrange,  dosetime, id)
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
function pdimport(data, time, obs, sort; bl = 0, th = 0, limitrule::Union{Nothing, LimitRule} = nothing, warn = true)
    if isa(sort, String) sort = [Symbol(sort)] end
    if isa(sort, Symbol) sort = [sort] end

    Tables.istable(data) || error("Data not a table!")

    cols   = Tables.columns(data)
    cdata  = Tuple(Tables.getcolumn(cols, y) for y in sort)
    d      = Dict{Tuple{eltype.(cdata)...}, Vector{Int}}()
    indsdict!(d, cdata)

    timec = Tables.getcolumn(data, time)
    obsc  = Tables.getcolumn(data, obs)

    any(isnanormissing, timec) && error("Some time values is NaN or Missing!")

    sdata = Vector{PDSubject}(undef, length(d))
    i = one(Int)
    @inbounds for (k, v) in d
        timevals = timec[v]
        obsvals  = obsc[v]

        if !allunique(timevals)

            nuv = nonunique(timevals)
            warn && @warn "Not all time values is unique ($nuv), last observation used! ($k)"

            nuvinds = findall(x -> x == first(nuv), timevals)
            resize!(nuvinds, length(nuvinds) - 1)
            if length(nuv) > 1
                for cnt = 2:length(nuv)
                    nuvinds_ = findall(x -> x == nuv[cnt], timevals)
                    resize!(nuvinds_, length(nuvinds_) - 1)
                    append!(nuvinds, nuvinds_)
                end
            end
            sort!(nuvinds)
            deleteat!(v, nuvinds)
            deleteat!(v, nuvinds)
            timevals = view(timec, v)
            obsvals = view(obsc, v)
        end
        sp = sortperm(timevals)
        timevals_spv = view(timevals, sp)
        obsvals_spv = view(obsvals, sp)
        timevals_sp, obsvals_sp = checkvalues(timevals_spv, obsvals_spv, warn = warn)
        sdata[i] = PDSubject(timevals_sp, obsvals_sp, bl, th, Dict(sort .=> k))
        i += one(Int)
    end
    ds = DataSet(identity.(sdata))
    if !isnothing(limitrule)
        applylimitrule!(ds, limitrule)
    end
    ds
end

"""
    pdimport(data, time, obs;
        warn = true,
        kwargs...)

Import PD data from tabular data `data`, `time` - time column, `obs` - observations column.
"""
function pdimport(data, time, obs; warn = true,  kwargs...)
    timevals_sp, obsvals_sp = checkvalues(Tables.getcolumn(data, time), Tables.getcolumn(data, obs), warn = warn)
    pdimport(timevals_sp, obsvals_sp; kwargs...)
end
"""
    pdimport(time, obs;
        bl = 0,
        th = 0,
        id = Dict{Symbol, Any}(),
        warn = true)

Import PD data from time vector `time` and observations vector `obs`.
"""
function pdimport(time, obs; bl = 0, th = 0, id = Dict{Symbol, Any}(), warn = true)
    timevals_sp, obsvals_sp = checkvalues(copy(time), copy(obs), warn = warn)
    PDSubject(timevals_sp, obsvals_sp, bl, th,  id)
end
