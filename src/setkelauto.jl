
"""
    setkelauto!(data::T, kelauto::Bool) where T <: PKSubject

Set range for elimination parameters calculation for subject.

data     - PK subject;
kelauto  - value.
"""
function setkelauto!(data::T, kelauto::Bool) where T <: PKSubject
    if data.kelrange.kelend  > 0 && data.kelrange.kelstart > 0 data.kelauto = kelauto end
    data
end
"""
    setkelauto!(data::DataSet{T}, kelauto::Bool, ind::Int) where T <: PKSubject
"""
function setkelauto!(data::DataSet{T}, kelauto::Bool, ind::Int) where T <: PKSubject
    setkelauto!(data[ind], kelauto)
    data
end
"""
    getkelauto!(data::T) where T <: PKSubject
"""
function getkelauto(data::T) where T <: PKSubject
    data.kelauto
end
