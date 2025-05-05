#=
function Base.append!(t::MetidaTable, t2::MetidaTable)
    if !(names(t) ⊆ names(t2)) error("Names for t not in t2") end
    for n in names(t)
        append!(t.table[n], t2.table[n])
    end
    t
end
=#
function MetidaBase.metida_table_(obj::DataSet{T}; obstime = false) where T <: Union{PKSubject, PDSubject}
    if obstime
        return metida_table_obstime(obj)
    else
        return metida_table_all(obj)
    end
end

function metida_table_all(obj::DataSet{T}) where T <: Union{PKSubject, PDSubject}
    idset  = Set(keys(first(obj).id))
    if length(obj) > 1
        for i = 2:length(obj)
            union!(idset,  Set(keys(obj[i].id)))
        end
    end
    mt1 = metida_table_((fill(getid(obj, 1, c), length(obj[1])) for c in idset)...; names = idset)
    mt2 = metida_table_(deepcopy(gettime(obj[1])), deepcopy(getobs(obj[1])); names = [:time, :obs])
    mtm = merge(mt1, mt2)
    if length(obj) > 1
        
        for i = 2:length(obj)
            mt1 = metida_table_((fill(getid(obj, i, c), length(obj[i])) for c in idset)...; names = idset)
            mt2 = metida_table_(gettime(obj[i]), getobs(obj[i]); names = [:time, :obs])
            amtm = merge(mt1, mt2)
            for n in keys(mtm)
                append!(mtm[n], amtm[n])
            end
        end
    end
    mtm
end

function metida_table_obstime(obj::DataSet{T}) where T <: Union{PKSubject, PDSubject}
    mtm = metida_table_(deepcopy(gettime(obj[1])), deepcopy(getobs(obj[1])); names = [:time, :obs])
    if length(obj) > 1
        for i = 2:length(obj)
            append!(mtm[:time], gettime(obj[i]))
            append!(mtm[:obs], getobs(obj[i]))
        end
    end
    mtm
end
