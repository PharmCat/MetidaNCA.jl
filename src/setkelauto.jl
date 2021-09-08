#Subject
"""
    setkelauto!(data::T, kelauto::Bool) where T <: PKSubject

Set range for elimination parameters calculation for subject.

* `data`     - PK subject;
* `kelauto`  - value.
"""
function setkelauto!(data::T, kelauto::Bool) where T <: PKSubject
    if !kelauto
        if !(data.kelrange.kelend > data.kelrange.kelstart > 0) error("Start point: $(data.kelrange.kelstart) end point: $(data.kelrange.kelend), check that data.kelrange.kelend > data.kelrange.kelstart > 0") end
    end
    data.kelauto = kelauto
    data
end

#DS ind Int
"""
    setkelauto!(data::DataSet{T}, kelauto::Bool, ind::Int) where T <: PKSubject
"""
function setkelauto!(data::DataSet{T}, kelauto::Bool, ind::Int) where T <: PKSubject
    setkelauto!(data[ind], kelauto)
    data
end
#DS iter Int
"""
    setkelauto!(data::DataSet{T}, kelauto::Bool, inds::Union{Vector{Int}, UnitRange{Int}, Tuple{Vararg{Int}}}) where T <: PKSubject
"""
function setkelauto!(data::DataSet{T}, kelauto::Bool, inds::Union{Vector{Int}, UnitRange{Int}, Tuple{Vararg{Int}}}) where T <: PKSubject
    for i in inds
        setkelauto!(data[i], kelauto)
    end
    data
end
#DS all
"""
    setkelauto!(data::DataSet{T}, kelauto::Bool) where T <: PKSubject
"""
function setkelauto!(data::DataSet{T}, kelauto::Bool) where T <: PKSubject
    for i = 1:length(data)
        setkelauto!(data[i], kelauto)
    end
    data
end
#DS Dict
"""
    setkelauto!(data::DataSet{T}, kelauto::Bool, sort::Dict) where T <: PKSubject
"""
function setkelauto!(data::DataSet{T}, kelauto::Bool, sort::Dict) where T <: PKSubject
    for i = 1:length(data)
        if sort ⊆ data[i].id setkelauto!(data[i], kelauto) end
    end
    data
end

#GET subj
"""
    getkelauto!(data::T) where T <: PKSubject
"""
function getkelauto(data::T) where T <: PKSubject
    data.kelauto
end
