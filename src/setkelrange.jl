#Subject
"""
    setkelrange!(data::T, range::ElimRange{:point}; kelauto = false) where T <: PKSubject

Set `range` for subject `data`. Set `kelauto` if possible.
"""
function setkelrange!(data::T, range::ElimRange{:point}; kelauto = false) where T <: PKSubject
    if range.kelend > length(data) throw(ArgumentError("Kel endpoint out of range")) end
    data.kelrange = range
    setkelauto!(data, kelauto)
    data
end
#DS ind Int
"""
    setdosetime!(data::DataSet{T}, dosetime::DoseTime, ind::Int) where T <: PKSubject
"""
function setkelrange!(data::DataSet{T}, range::ElimRange{:point}, ind::Int; kelauto = false) where T <: PKSubject
    setkelrange!(data[ind], range; kelauto = kelauto)
    data
end
#DS iter Int
"""
    setkelrange!(data::DataSet{T}, range::ElimRange{:point}, inds::Union{Vector{Int}, UnitRange{Int}, Tuple{Vararg{Int}}}; kelauto = false)
"""
function setkelrange!(data::DataSet{T}, range::ElimRange{:point}, inds::Union{Vector{Int}, UnitRange{Int}, Tuple{Vararg{Int}}}; kelauto = false) where T <: PKSubject
    for i in inds
        setkelrange!(data[i], range; kelauto = kelauto)
    end
    data
end
#DS all
"""
    setkelrange!(data::DataSet{T}, range::ElimRange{:point}; kelauto = false) where T <: PKSubject
"""
function setkelrange!(data::DataSet{T}, range::ElimRange{:point}; kelauto = false) where T <: PKSubject
    for i = 1:length(data)
        setkelrange!(data[i], range; kelauto = kelauto)
    end
    data
end
#DS Dict
"""
    setkelrange!(data::DataSet{T}, range::ElimRange{:point}, sort::Dict; kelauto = false) where T <: PKSubject
"""
function  setkelrange!(data::DataSet{T}, range::ElimRange{:point}, sort::Dict; kelauto = false) where T <: PKSubject
    for i = 1:length(data)
        if sort ⊆ data[i].id setkelrange!(data[i], range; kelauto = kelauto) end
    end
    data
end
#GET subj
"""
    getkelrange(data::T) where T <: PKSubject
"""
function getkelrange(data::T) where T <: PKSubject
    data.kelrange
end
