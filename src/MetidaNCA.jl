# MetidaNCA
# Copyright © 2019-2025 Vladimir Arnautov aka PharmCat <mail@pharmcat.net>
# SPDX-License-Identifier: MIT
module MetidaNCA

using RecipesBase

import RecipesBase: plot!, plot
import Statistics: mean, quantile
import Base: length, push!, resize!, ht_keyindex, convert, first, display, show
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
pkplot, vpcplot,
getkeldata, getkelauto, getkelrange, getdosetime, getbl, getth, subset,
metida_table,
PKSubject, UPKSubject, PDSubject, NCAResult

function __init__()
    @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
        savefig = Plots.savefig
        current = Plots.current

        function mergeplots!(sp1::Plots.Subplot, sp2::Plots.Subplot)
            append!(sp1.series_list, sp2.series_list)
            Plots.expand_extrema!(sp1[:xaxis], Plots.xlims(sp2))
            Plots.expand_extrema!(sp1[:yaxis], Plots.ylims(sp2))
            Plots.expand_extrema!(sp1[:zaxis], Plots.zlims(sp2))
            return sp1
        end
        
        function mergeplots!(plt, plts...)
            for (i, sp) in enumerate(plt.subplots)
                for other_plt in plts
                    if i in eachindex(other_plt.subplots)
                      mergeplots!(sp, other_plt[i])
                    end
                end
            end
            return plt
        end
        function Plots.png(plt::PKPlot, io::IO)
            Plots.png(plt.plot, io)
        end
        function Plots.png(plt::PKPlot, fn)
            Plots.png(plt.plot, fn)
        end
        function Plots.png(plt::PKPlot)
            Plots.png(plt.plot)
        end
    end
end

function mergeplots! end

function png end

const LOG2 = log(2)

const NCARESOBS = Symbol("")

    include("types.jl")
    include("setkelauto.jl")
    include("setkelrange.jl")
    include("setdosetime.jl")
    include("getkeldata.jl")
    include("applylimitrule.jl")
    include("import.jl")
    include("nca.jl")
    include("plots.jl")
    include("metidatable.jl")
    include("setblth.jl")
    include("timefilter.jl")
    include("dropnanormissing.jl")
    include("sparse.jl")
    include("atomic.jl")
    include("precompile.jl")
    include("show.jl")

end # module
