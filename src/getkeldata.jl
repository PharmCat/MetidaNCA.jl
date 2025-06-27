function getkeldata(data::T) where T <: PKSubject
    error("Method deprecated...")
end
"""
    getkeldata(data::T) where T <: PKSubject
"""
function getkeldata(ncar::T) where T <: NCAResult
    ncar.keldata
end
