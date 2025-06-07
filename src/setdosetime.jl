#Subject
"""
    setdosetime!(data::AbstractSubject, dosetime)

Set dose time `dosetime` for subject `data`.
"""
function setdosetime!(data::AbstractSubject, dosetime)
    data.dosetime = dosetime
    data
end


#Subject
"""
    setdosetime!(data::PKSubject, dosetime::AbstractVector{<:DoseTime})

Set dose time `dosetime` for subject `data`.
"""
function setdosetime!(data::PKSubject, dosetime::AbstractVector{<:DoseTime})
    data.dosetime = dosetime
    data
end

#Subject
"""
    setdosetime!(data::PKSubject, dosetime::DoseTime)

Set dose time `dosetime` for subject `data`.
"""
function setdosetime!(data::PKSubject, dosetime::DoseTime)
    setdosetime!(data, [dosetime])
    data
end

#DS ind Int
"""
    setdosetime!(data::DataSet{T}, dosetime::DoseTime, ind::Int) where T <: AbstractSubject

* `ind` - index in DataSet.
"""
function setdosetime!(data::DataSet{<:AbstractSubject}, dosetime::DoseTime, ind::Int)
    setdosetime!(data[ind], dosetime)
    data
end
#DS iter Int
"""
    setdosetime!(data::DataSet{T}, dosetime::DoseTime, inds::Union{Vector{Int}, UnitRange{Int}, Tuple{Vararg{Int}}}) where T <: PKSubject

* `inds` - indexes in DataSet.
"""
function setdosetime!(data::DataSet{<:AbstractSubject}, dosetime::DoseTime, inds::Union{Vector{Int}, UnitRange{Int}, Tuple{Vararg{Int}}})
    for i in inds
        setdosetime!(data[i], dosetime)
    end
    data
end
#DS all
"""
    setdosetime!(data::DataSet{T}, dosetime::DoseTime) where T <: AbstractSubject

For all subjects in DataSet.
"""
function setdosetime!(data::DataSet{<:AbstractSubject}, dosetime::DoseTime)
    for i = 1:length(data)
        setdosetime!(data[i], dosetime)
    end
    data
end
#DS Dict
"""
    setdosetime!(data::DataSet{T}, dosetime::DoseTime, sort::Dict) where T <: AbstractSubject

Set dose time `dosetime` for subjects if `sort` ⊆ subject's `id`.
"""
function setdosetime!(data::DataSet{<:AbstractSubject}, dosetime::DoseTime, sort::Dict)
    for i = 1:length(data)
        if sort ⊆ data[i].id setdosetime!(data[i], dosetime) end
    end
    data
end
#GET subj
"""
    getdosetime(data::T) where T <: PKSubject

Return dosetime.
"""
function getdosetime(data::AbstractSubject)
    data.dosetime
end


#Subject
#DS ind Int
#DS iter Int
#DS all
#DS Dict
#GET subj
