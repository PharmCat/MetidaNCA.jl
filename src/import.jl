# Заполняет словарь d индексами индивидуальных значений
#=
function indsdict!(d::Dict{T}, cdata::Tuple) where T
    @inbounds for (i, element) in enumerate(zip(cdata...))
        ind = ht_keyindex(d, element)
        if ind > 0
            push!(d.vals[ind], i)
        else
            v = Vector{Int}(undef, 1)
            v[1] = i
            d[element] = v
        end
    end
    d
end
=#
nonunique(v) = [k for (k, v) in StatsBase.countmap(v) if v > 1]

function floatparse(data)
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

function checkvalues(timevals_sp, concvals_sp)
    if !(eltype(timevals_sp) <: Number)
        timevals_sp = identity.(timevals_sp)
        eltype(timevals_sp) <: Number || error("Some time values not a number!")
    end
    if !(eltype(concvals_sp) <: Union{Number, Missing})
        concvals_sp = identity.(concvals_sp)
        if !(eltype(concvals_sp) <: Union{Number, Missing})
            @warn "Some concentration values not a number, try to fix."
            concvals_sp = floatparse(concvals_sp)
        end
    end
    timevals_sp, concvals_sp
end

"""
    pkimport(data, time, conc, sort; kelauto = true,  elimrange = ElimRange(), dosetime = DoseTime())

Import PK data from table `data`.

* `time` - time column;
* `conc` - concentration column;
* `sort` - subject sorting columns.

keywords:

* `kelauto` - if `true` auto range settings, if `false` used `kelstart`/`kelend` from `elimrange`;
* `elimrange` - set elimination range settings;
* `dosetime` - set dose and dose time, by default dosetime = 0, dose is `NaN`.

!!! note

    If time column have non-unique values - last pair time-concentration will be used.

"""
function pkimport(data, time, conc, sort; kelauto = true,  elimrange = ElimRange(), dosetime = DoseTime())
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
        timevals = timec[v]
        concvals = concc[v]
        if !allunique(timevals)
            @warn "Not all time values is unique, last observation used! ($k)"
            nuv = nonunique(timevals)
            nuvinds = findall(x -> x == first(nuv), timevals)[1:end-1]
            if length(nuv) > 1
                for cnt = 2:length(nuv)
                    append!(nuvinds, findall(x -> x == nuv[cnt], timevals)[1:end-1])
                end
            end
            sort!(nuvinds)
            deleteat!(timevals, nuvinds)
            deleteat!(concvals, nuvinds)
        end
        sp = sortperm(timevals)
        timevals_sp = timevals[sp]
        concvals_sp = concvals[sp]
        timevals_sp, concvals_sp = checkvalues(timevals_sp, concvals_sp)
        sdata[i] = PKSubject(timevals_sp, concvals_sp, kelauto, elimrange,  dosetime, Dict(sort .=> k))
        i += one(Int)
    end
    return DataSet(identity.(sdata))
end
"""
    pkimport(data, time, conc; kelauto = true,  elimrange = ElimRange(), dosetime = DoseTime())

Import PK data from tabular data `data`, `time` - time column, `conc` - concentration column.
"""
function pkimport(data, time, conc; kelauto = true,  elimrange = ElimRange(), dosetime = DoseTime())
    timevals_sp, concvals_sp = checkvalues(copy(Tables.getcolumn(data, time)), copy(Tables.getcolumn(data, conc)))
    pkimport(timevals_sp, concvals_sp; kelauto = kelauto,  elimrange = elimrange, dosetime = dosetime)
end
"""
    pkimport(time, conc; kelauto = true,  elimrange = ElimRange(), dosetime = DoseTime())

Import PK data from time vector `time` and concentration vector `conc`.
"""
function pkimport(time, conc; kelauto = true,  elimrange = ElimRange(), dosetime = DoseTime())
    timevals_sp, concvals_sp = checkvalues(copy(time), copy(conc))
    PKSubject(timevals_sp, concvals_sp, kelauto, elimrange,  dosetime, Dict{Symbol, Any}())
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
        stimevals = stimec[v]
        etimevals = etimec[v]
        concvals  = concc[v]
        volvals   = volc[v]

        timeranges = collect(zip(stimevals, etimevals))
        sp = sortperm(stimevals)

        timevals_sp = timeranges[sp]
        concvals_sp = concvals[sp]
        volvals_sp  = volvals[sp]

        if length(timevals_sp) > 1
            for c = 2:length(timevals_sp)
                timevals_sp[c][1] == timevals_sp[c-1][2] || error("Start time ($(timevals_sp[c][1])) for observation $c not equal End time ($(timevals_sp[c-1][2])) for observation $(c-1)!")
            end
        end
        sdata[i] = UPKSubject(timevals_sp, concvals_sp, volvals_sp, kelauto, elimrange,  dosetime, Dict(sort .=> k))
        i += one(Int)
    end
    return DataSet(identity.(sdata))
end
