


function pkimport(data, time, conc, sort)
    cols   = Tables.columns(data)
    cdata  = Tuple(Tables.getcolumn(cols, y) for y in Tuple(sort))
    ctype  = Tuple{eltype.(cdata)...}
    d      = Dict{ctype, Set{Int}}()
    @inbounds for (i, element) in enumerate(zip(cdata...))
        ind = ht_keyindex(d, element)
        if ind > 0
            push!(d.vals[ind], i)
        else
            d[element] = Set{Int}(i)
        end
    end
    sdata = Vector{PKSubject}(undef, length(d))
    timec = Tables.getcolumn(data, time)
    concc = Tables.getcolumn(data, conc)
    i = one(Int)
    @inbounds for (k, v) in d
        cv = collect(v)
        timevals = timec[cv]
        concvals = concc[cv]
        if !allunique(timevals) @warn "Not all time values is unique!" end
        sp = sortperm(timevals)
        sdata[i] = PKSubject(timevals[sp], concvals[sp], Dict(sort .=> k))
        i += one(Int)
    end
    return DataSet(sdata)
end
