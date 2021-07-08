#Subject
"""
    setdosetime!(data::T, dosetime::DoseTime) where T <: PKSubject
"""
function setdosetime!(data::T, dosetime::DoseTime) where T <: PKSubject
    data.dosetime = dosetime
    data
end
#DS ind Int
"""
    setdosetime!(data::DataSet{T}, dosetime::DoseTime, ind::Int) where T <: PKSubject
"""
function setdosetime!(data::DataSet{T}, dosetime::DoseTime, ind::Int) where T <: PKSubject
    setdosetime!(data[ind], dosetime)
    data
end
#DS iter Int
"""
    setdosetime!(data::DataSet{T}, dosetime::DoseTime, inds::Union{Vector{Int}, UnitRange{Int}, Tuple{Vararg{Int}}}) where T <: PKSubject
"""
function setdosetime!(data::DataSet{T}, dosetime::DoseTime, inds::Union{Vector{Int}, UnitRange{Int}, Tuple{Vararg{Int}}}) where T <: PKSubject
    for i in inds
        setdosetime!(data[i], dosetime)
    end
    data
end
#DS all
"""
    setdosetime!(data::DataSet{T}, dosetime::DoseTime) where T <: PKSubject
"""
function setdosetime!(data::DataSet{T}, dosetime::DoseTime) where T <: PKSubject
    for i = 1:length(data)
        setdosetime!(data[i], dosetime)
    end
    data
end
#DS Dict
"""
    setdosetime!(data::DataSet{T}, dosetime::DoseTime, sort::Dict) where T <: PKSubject
"""
function setdosetime!(data::DataSet{T}, dosetime::DoseTime, sort::Dict) where T <: PKSubject
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
function getdosetime(data::T) where T <: PKSubject
    data.dosetime
end


#Subject
#DS ind Int
#DS iter Int
#DS all
#DS Dict
#GET subj
