# Metida
# Copyright © 2019-2020 Vladimir Arnautov aka PharmCat <mail@pharmcat.net>

module MetidaNCA

using RecipesBase

import Base: length, length, push!, resize!
import MetidaBase
import MetidaBase: Tables, StatsBase, SnoopPrecompile,
PrettyTables,
AbstractIdData,
AbstractSubject,
DataSet,
AbstractSubjectResult,
AbstractResultData,
isnanormissing,
getid,
getdata,
metida_table, metida_table_, MetidaTable,
uniqueidlist,
indsdict!,
subset

export pkimport, upkimport, pdimport, nca!, nca, DoseTime, ElimRange, LimitRule, NoPageSort,
auc_sparse,
setdosetime!, setkelauto!, setkelrange!, applylimitrule!, setbl!, setth!,
pkplot,
getkeldata, getkelauto, getkelrange, getdosetime, getbl, getth, subset
metida_table

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
    include("metidatable.jl")
    include("setblth.jl")
    include("timefilter.jl")
    include("sparse.jl")
    include("atomic.jl")
    include("precompile.jl")

end # module
