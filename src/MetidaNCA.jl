# Metida
# Copyright © 2019-2020 Vladimir Arnautov aka PharmCat <mail@pharmcat.net>

module MetidaNCA

using RecipesBase

import Base: ht_keyindex, length, length, push!

import MetidaBase: Tables, PrettyTables, AbstractIdData, AbstractSubject, DataSet, AbstractSubjectResult, AbstractResultData, getid, metida_table

export pkimport, nca!, DoseTime, ElimRange, LimitRule, setdosetime!, setkelauto!, setkelrange!, applylimitrule!, pkplot

const LOG2 = log(2)

    include("types.jl")
    include("setkelauto.jl")
    include("setkelrange.jl")
    include("setdosetime.jl")
    include("getkeldata.jl")
    include("applylimitrule.jl")
    include("show.jl")
    include("import.jl")
    include("nca.jl")
    include("plots.jl")

end # module
