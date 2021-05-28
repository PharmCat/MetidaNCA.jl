
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
function pkimport(data, time, conc, sort)
    cols   = Tables.columns(data)
    cdata  = Tuple(Tables.getcolumn(cols, y) for y in sort)
    d      = Dict{Tuple{eltype.(cdata)...}, Vector{Int}}()
    indsdict!(d, cdata)
    sdata = Vector{PKSubject}(undef, length(d))
    timec = Tables.getcolumn(data, time)
    concc = Tables.getcolumn(data, conc)
    i = one(Int)
    @inbounds for (k, v) in d
        timevals = timec[v]
        concvals = concc[v]
        if !allunique(timevals) @warn "Not all time values is unique!" end
        sp = sortperm(timevals)
        sdata[i] = PKSubject(timevals[sp], concvals[sp], Dict(sort .=> k))
        i += one(Int)
    end
    return DataSet(sdata)
end
