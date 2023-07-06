# BL
# Subject
"""
    setbl!(data::T, bl) where T <: PDSubject

Set `baseline` for pd subject `data`.
"""
function setbl!(data::T, bl) where T <: PDSubject
    isnanormissing(bl) && error("Baseline can't be NaN or missing.")
    data.bl = bl
    data
end
#DS ind Int
"""
    setbl!(data::DataSet{T}, bl, ind::Int) where T <: PDSubject

Set baseline for subject `ind` in `data`.
"""
function setbl!(data::DataSet{T}, bl, ind::Int) where T <: PDSubject
    setbl!(data[ind], bl)
    data
end
#DS iter Int
"""
    setbl!(data::DataSet{T}, bl, inds::Union{Vector{Int}, UnitRange{Int}, Tuple{Vararg{Int}}})
 
Set baseline for all subjects in range or vector `ind` in `data`.
"""
function setbl!(data::DataSet{T}, bl, inds::Union{Vector{Int}, UnitRange{Int}, Tuple{Vararg{Int}}}) where T <: PDSubject
    for i in inds
        setbl!(data[i], bl)
    end
    data
end
#DS all
"""
    setbl!(data::DataSet{T}, bl) where T <: PDSubject

Set `baseline` for all subjects in `data`.
"""
function setbl!(data::DataSet{T}, bl) where T <: PDSubject
    for i = 1:length(data)
        setbl!(data[i], bl)
    end
    data
end
#DS Dict
"""
    setbl!(data::DataSet{T}, bl, sort::Dict) where T <: PDSubject

Set `baseline` only for subjects which `sort` ⊆ `id` is `true`.
"""
function  setbl!(data::DataSet{T}, bl, sort::Dict; kelauto = false) where T <: PDSubject
    for i = 1:length(data)
        if sort ⊆ data[i].id setbl!(data[i], bl) end
    end
    data
end
#GET subj
"""
    getbl(data::T) where T <: PDSubject
"""
function getbl(data::T) where T <: PDSubject
    data.bl
end
################################################################################
# TH
# Subject
"""
    setth!(data::T, th) where T <: PDSubject

Set `threshold` for subject `data`.
"""
function setth!(data::T, th) where T <: PDSubject
    isnanormissing(th) && error("Threshold can't be NaN or missing.")
    data.th = th
    data
end
#DS ind Int
"""
    setth!(data::DataSet{T}, th, ind::Int) where T <: PDSubject
"""
function setth!(data::DataSet{T}, th, ind::Int) where T <: PDSubject
    setth!(data[ind], th)
    data
end
#DS iter Int
"""
    setth!(data::DataSet{T}, th, inds::Union{Vector{Int}, UnitRange{Int}, Tuple{Vararg{Int}}})
"""
function setth!(data::DataSet{T}, th, inds::Union{Vector{Int}, UnitRange{Int}, Tuple{Vararg{Int}}}) where T <: PDSubject
    for i in inds
        setth!(data[i], th)
    end
    data
end
#DS all
"""
    setth!(data::DataSet{T}, th) where T <: PDSubject
"""
function setth!(data::DataSet{T}, th) where T <: PDSubject
    for i = 1:length(data)
        setth!(data[i], th)
    end
    data
end
#DS Dict
"""
    setth!(data::DataSet{T}, th, sort::Dict) where T <: PDSubject
"""
function  setth!(data::DataSet{T}, th, sort::Dict; kelauto = false) where T <: PDSubject
    for i = 1:length(data)
        if sort ⊆ data[i].id setth!(data[i], th) end
    end
    data
end
#GET subj
"""
    getth(data::T) where T <: PDSubject
"""
function getth(data::T) where T <: PDSubject
    data.th
end
