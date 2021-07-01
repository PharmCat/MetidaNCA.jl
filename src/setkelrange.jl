"""
    setkelrange!(data::T, range::ElimRange{:point}; kelauto = false) where T <: PKSubject
"""
function setkelrange!(data::T, range::ElimRange{:point}; kelauto = false) where T <: PKSubject
    if range.kelend > length(data) throw(ArgumentError("Kel endpoint out of range")) end
    if range.kelend > range.kelstart > 0 setkelauto!(data, kelauto) end
    data.kelrange = range
    data
end
"""
    setdosetime!(data::DataSet{T}, dosetime::DoseTime, ind::Int) where T <: PKSubject
"""
function setkelrange!(data::DataSet{T}, range::ElimRange{:point}, ind::Int; kelauto = false) where T <: PKSubject
    setkelrange!(data[ind], range)
    data
end
"""
    getkelrange(data::T) where T <: PKSubject
"""
function getkelrange(data::T) where T <: PKSubject
    data.kelrange
end
