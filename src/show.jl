
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
    m = hcat(copy(obj.s), obj.e, obj.a, obj.b, obj.r, obj.ar)
    println(io, "Elimination table:")
    PrettyTables.pretty_table(io, m; header = ["Strat time", "End time", "a", "b", "r²", "Adjusted r²"])
end

# PK Subject

function Base.show(io::IO, obj::PKSubject)
    println(io, "  Pharmacokinetic subject")
    println(io, "Observations: $(length(obj)); ", obj.dosetime)
    println(io,  obj.kelrange)
    println(io, "Time   Concentration")
    for i = 1:length(obj)
        println(io, obj.time[i], " => ", obj.obs[i])
    end
end
function Base.show(io::IO, obj::DataSet{Subj}) where Subj <: PKSubject
    println(io, "DataSet: Pharmacokinetic subject")
    println(io, "Length: $(length(obj))")
    for i = 1:length(obj)
        println(io, "Subject $(i): ", obj[i].id)
    end
end
