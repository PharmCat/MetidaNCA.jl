#=
function Base.append!(t::MetidaTable, t2::MetidaTable)
    if !(names(t) ⊆ names(t2)) error("Names for t not in t2") end
    for n in names(t)
        append!(t.table[n], t2.table[n])
    end
    t
end
=#
function MetidaBase.metida_table(obj::DataSet{T}) where T <: PKSubject
    idset  = Set(keys(first(obj).id))
    if length(obj) > 1
        for i = 2:length(obj)
            union!(idset,  Set(keys(obj[i].id)))
        end
    end
    mt1 = MetidaBase.metida_table((fill(getid(obj, 1, c), length(obj[1])) for c in idset)...; names = idset)
    mt2 = MetidaBase.metida_table(deepcopy(obj[1].time), deepcopy(obj[1].obs); names = [:time, :obs])
    mtm = MetidaTable(merge(mt1.table, mt2.table))
    if length(obj) > 1
        for i = 2:length(obj)
            mt1 = MetidaBase.metida_table((fill(getid(obj, i, c), length(obj[i])) for c in idset)...; names = idset)
            mt2 = MetidaBase.metida_table(obj[i].time, obj[i].obs; names = [:time, :obs])
            append!(mtm, MetidaTable(merge(mt1.table, mt2.table)))
        end
    end
    mtm
end
