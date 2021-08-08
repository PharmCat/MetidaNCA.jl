
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

"""
function pkimport(data, time, conc, sort; kelauto = true,  elimrange = ElimRange(), dosetime = DoseTime())
    if isa(sort, String) sort = [Symbol(sort)] end
    if isa(sort, Symbol) sort = [sort] end
    
    cols   = Tables.columns(data)
    cdata  = Tuple(Tables.getcolumn(cols, y) for y in sort)
    d      = Dict{Tuple{eltype.(cdata)...}, Vector{Int}}()
    indsdict!(d, cdata)

    timec = Tables.getcolumn(data, time)
    concc = Tables.getcolumn(data, conc)
    sdata = Vector{PKSubject}(undef, length(d))
    i = one(Int)
    @inbounds for (k, v) in d
        timevals = timec[v]
        concvals = concc[v]
        if !allunique(timevals) @warn "Not all time values is unique!" end
        sp = sortperm(timevals)
        sdata[i] = PKSubject(timevals[sp], concvals[sp], kelauto, elimrange,  dosetime, Dict(sort .=> k))
        i += one(Int)
    end
    return DataSet(identity.(sdata))
end
"""
    pkimport(data, time, conc; kelauto = true,  elimrange = ElimRange(), dosetime = DoseTime())

Import PK data from tabular data `data`, `time` - time column, `conc` - concentration column.
"""
function pkimport(data, time, conc; kelauto = true,  elimrange = ElimRange(), dosetime = DoseTime())
    pkimport(copy(Tables.getcolumn(data, time)), copy(Tables.getcolumn(data, conc)); kelauto = kelauto,  elimrange = elimrange, dosetime = dosetime)
end
"""
    pkimport(time, conc; kelauto = true,  elimrange = ElimRange(), dosetime = DoseTime())

Import PK data from time vector `time` and concentration vector `conc`.
"""
function pkimport(time, conc; kelauto = true,  elimrange = ElimRange(), dosetime = DoseTime())
    PKSubject(copy(time), copy(conc), kelauto, elimrange,  dosetime, Dict{Symbol, Any}())
end
