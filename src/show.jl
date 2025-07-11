
function Base.show(io::IO, obj::DoseTime)
    print(io, "Dose - $(obj.dose); Time - $(obj.time); Tau - $(obj.tau)")

end
function Base.show(io::IO, obj::ElimRange)
    print(io, "Elimination range: $(obj.kelstart) - $(obj.kelend) ")
    if length(obj.kelexcl) > 0
        print(io, "Exclusions: $(obj.kelexcl[1])")
        if length(obj.kelexcl) > 1 
            for i = 2:length(obj.kelexcl) print(io, ", $(obj.kelexcl[i])") end 
        end
        print(io, ".")
    else
        print(io, "No exclusion.")
    end
end
function Base.show(io::IO, obj::KelData)
    println(io, "Elimination table:")
    header = ["Strat time", "End time", "a", "b", "r²", "Adjusted r²", "N"]
    mt = metida_table(obj.s, obj.e, obj.a, obj.b, obj.r, obj.ar, obj.n; names =  Tuple(Symbol.(header)))
    PrettyTables.pretty_table(io, mt; tf = PrettyTables.tf_compact, header = header)
end

# PK Subject

function Base.show(io::IO, obj::PKSubject)
    print(io, "  Pharmacokinetic subject")
    if isnothing(obsnames(obj))
        println(io, "")
    else
        on = obsnames(obj)
        print(io, " (Observation names: $(on[1])")
        if length(on) > 1
            for i = 2:length(on)
                print(io, ", $(on[i])")
            end
        end
        println(io, ")")
    end
    if length(obj.id) > 0
        print(io, "ID: ")
        for (k, v) in obj.id
            print(io, "$k => $v;")
        end
        println(io, "")
    end
    print(io, "Observations: $(length(obj)); ")
    if isa(obj.dosetime, Vector)
        println(io,"")
        for i in 1:length(obj.dosetime)
            println(io, "Dose $i: ", obj.dosetime[i])
        end
    else
        println(io, obj.dosetime)
    end
    println(io,  obj.kelrange)
    if ismultobs(obj)
        if obsnumber(obj) == 1
            PrettyTables.pretty_table(io, metida_table(obj.time, getobs(obj); names = (:Time, getfirstkey(obj))); tf = PrettyTables.tf_compact)
        else
            PrettyTables.pretty_table(io, metida_table(obj.time, obj.obs...; names = append!([:Time], keys(obj.obs))); tf = PrettyTables.tf_compact)
        end
    else
        PrettyTables.pretty_table(io, metida_table(obj.time, getobs(obj); names = (:Time, :Observations)); tf = PrettyTables.tf_compact)
    end
end
function Base.show(io::IO, obj::UPKSubject)
    println(io, "  Pharmacokinetic subject (urine)")
    println(io, "Observations: $(length(obj)); ", obj.dosetime)
    println(io,  obj.kelrange)
    PrettyTables.pretty_table(io, metida_table(getindex.(obj.time, 1), getindex.(obj.time, 2), obj.obs, obj.vol); tf = PrettyTables.tf_compact, header  = ["Start time", "End time", "Concentration", "Volume"])

end

function Base.show(io::IO, obj::PDSubject)
    println(io, "  Pharmacodynamics subject")
    if length(obj.id) > 0
        print(io, "ID: ")
        for (k, v) in obj.id
            print(io, "$k => $v;")
        end
        println(io, "")
    end
    println(io, "Observations: $(length(obj)); ")
    PrettyTables.pretty_table(io, metida_table(obj.time, obj.obs; names = (:Time, :Observation)); tf = PrettyTables.tf_compact)
end

function subject_type_str(subj::Type{PKS}) where PKS <: PKSubject
    "Pharmacokinetics subject"
end
function subject_type_str(subj::Type{UPKS}) where UPKS <: UPKSubject
    "Pharmacokinetics subject (urine)"
end
function subject_type_str(subj::Type{PDS}) where PDS <: PDSubject
    "Pharmacodynamics subject"
end
function Base.show(io::IO, obj::DataSet{ST}) where ST <: AbstractSubject
    println(io, "DataSet: $(subject_type_str(ST))")
    println(io, "Length: $(length(obj))")
    lo = min(length(obj), 20)
    for i = 1:lo
        print(io, "Subject $(i): ")
        if length(obj[i].id) > 0
            for (k, v) in obj[i].id
                print(io, "$k => $v, ")
            end
            println(io, "")
        else
            println(io, "-")
        end
    end
    if lo < length(obj) 
        printstyled(io, "$(length(obj) - lo) subjects omitted... \n"; color = :blue)
    end
end

function Base.show(io::IO, obj::T) where T <: NCAResult
    print(io, "  PK/PD subject NCA result")
    if (obj.obsname != NCARESOBS) && !isnothing(obj.obsname)
        printstyled(io, " (NCA obs: $(obj.obsname))\n"; color = :blue)
    else
        println(io, "")
    end
    PrettyTables.pretty_table(io, obj.result; header = ["Parameter", "Value"], tf = PrettyTables.tf_compact)
end

function Base.show(io::IO, obj::DataSet{Res}) where Res <: NCAResult
    println(io, "DataSet: PK/PD NCA result")
    println(io, "Length: $(length(obj))")
    lo = min(length(obj), 20)
    for i = 1:lo
        print(io, "Subject $(i): ")
        if length(obj[i].data.id) > 0
            for (k, v) in obj[i].data.id
                print(io, "$k => $v, ")
            end
            println(io, "")
        else
            println(io, "-")
        end
    end
    if lo < length(obj) 
        printstyled(io, "$(length(obj) - lo) subjects omitted... \n"; color = :blue)
    end
end

function Base.show(io::IO, obj::PKPlot) 
    show(io, obj.plot) 
end
function Base.display(obj::PKPlot) 
    display(obj.plot) 
end
function Base.display(m::MIME, obj::PKPlot) 
    display(m, obj.plot) 
end
function Base.display(d::AbstractDisplay, mime::AbstractString, obj::PKPlot) 
    display(d, mime, obj.plot) 
end


function Base.show(io::IO, obj::DataSet{PKPlot}) 
    println(io, "DataSet: PKPlot")
    println(io, "Length: $(length(obj))")
    lo = min(length(obj), 20)
    for i = 1:lo
        print(io, "Subject $(i) plot: ")
        if length(obj[i].id) > 0
            for (k, v) in obj[i].id
                print(io, "$k => $v, ")
            end
            println(io, "")
        else
            println(io, "-")
        end
    end
    if lo < length(obj) 
        printstyled(io, "$(length(obj) - lo) subjects omitted... \n"; color = :blue)
    end
end
