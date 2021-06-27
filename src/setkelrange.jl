

function setkelrange!(data::T, range::ElimRange{:point}; kelauto = false) where T <: PKSubject
    if range.kelend > length(data) throw(ArgumentError("Kel endpoint out of range")) end
    if range.kelend > range.kelstart > 0 setkelauto!(data, kelauto) end
    data.kelrange = range
end
