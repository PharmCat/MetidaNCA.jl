"""
    getkeldata(data::T) where T <: PKSubject
"""
function getkeldata(data::T) where T <: PKSubject
    data.keldata
end
"""
    getkeldata(data::T) where T <: PKSubject
"""
function getkeldata(ncar::T) where T <: NCAResult
    getkeldata(ncar.data)
end
