"""
    getkeldata(data::T) where T <: PKSubject
"""
function getkeldata(data::T) where T <: PKSubject
    data.keldata
end
"""
    getkeldata(data::T) where T <: PKSubject
"""
function getkeldata(data::T) where T <: NCAResult
    getkeldata(data.subject)
end
