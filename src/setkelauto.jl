

"""
    setkelauto!(data::T, kelauto::Bool) where T <: PKSubject 

Set range for elimination parameters calculation for subject.

data     - PK subject;
kelauto  - value.
"""
function setkelauto!(data::T, kelauto::Bool) where T <: PKSubject
    data.kelauto = kelauto
end
