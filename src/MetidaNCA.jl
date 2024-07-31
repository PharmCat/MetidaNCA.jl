# Metida
# Copyright © 2019-2020 Vladimir Arnautov aka PharmCat <mail@pharmcat.net>

module MetidaNCA

using RecipesBase

import Base: length, length, push!, resize!
import MetidaBase
import MetidaBase: Tables, StatsBase, PrecompileTools,
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

using MetidaBase.Requires

export pkimport, upkimport, pdimport, nca!, nca, DoseTime, ElimRange, LimitRule, NoPageSort,
auc_sparse,
setdosetime!, setkelauto!, setkelrange!, applylimitrule!, setbl!, setth!,
pkplot,
getkeldata, getkelauto, getkelrange, getdosetime, getbl, getth, subset,
metida_table,
PKSubject, UPKSubject, PDSubject, NCAResult

function __init__()
    @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
        savefig = Plots.savefig
    end
end

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
