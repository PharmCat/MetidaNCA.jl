
function Base.show(io::IO, obj::DoseTime)
    print(io, "Dose - $(obj.dose); Time - $(obj.time); Tau - $(obj.tau)")

end
function Base.show(io::IO, obj::ElimRange)
    print(io, "Elimination range: $(obj.kelstart) - $(obj.kelend) ")
    if length(obj.kelexcl) > 0
        print(io, "Exclusions: $(obj.kelexcl[1])")
        if length(obj.kelexcl) > 1 for i = 2:length(obj.kelexcl) print(io, ", $(obj.kelexcl[i])") end end
        print(io, ".")
    else
        print(io, "No exclusion.")
    end
end
function Base.show(io::IO, obj::KelData)
    println(io, "Elimination table:")
    mt = metida_table(obj.s, obj.e, obj.a, obj.b, obj.r, obj.ar; names =  Tuple(Symbol.(["Strat time", "End time", "a", "b", "r²", "Adjusted r²"])))
    PrettyTables.pretty_table(io, mt; tf = PrettyTables.tf_compact)
end

# PK Subject

function Base.show(io::IO, obj::PKSubject)
    println(io, "  Pharmacokinetic subject")
    println(io, "Observations: $(length(obj)); ", obj.dosetime)
    println(io,  obj.kelrange)
    PrettyTables.pretty_table(io, metida_table(obj.time, obj.obs; names = (:Time, :Concentration)); tf = PrettyTables.tf_compact)

end
function Base.show(io::IO, obj::DataSet{Subj}) where Subj <: PKSubject
    println(io, "DataSet: Pharmacokinetic subject")
    println(io, "Length: $(length(obj))")
    for i = 1:length(obj)
        println(io, "Subject $(i): ", obj[i].id)
    end
end

function Base.show(io::IO, obj::T) where T <: NCAResult
    println(io, "  Pharmacokinetic subject NCA result")
    PrettyTables.pretty_table(io, obj.result; header = ["Parameter", "Value"], tf = PrettyTables.tf_compact)
end

function Base.show(io::IO, obj::DataSet{Res}) where Res <: NCAResult
    println(io, "DataSet: Pharmacokinetic subject NCA result")
    println(io, "Length: $(length(obj))")
    for i = 1:length(obj)
        println(io, "Subject $(i): ", obj[i].id)
    end
end
