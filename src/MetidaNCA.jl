# Metida
# Copyright © 2019-2020 Vladimir Arnautov aka PharmCat <mail@pharmcat.net>

module MetidaNCA

import Base: ht_keyindex, length, length, push!
import MetidaBase: AbstractIdData, AbstractSubject, DataSet, AbstractSubjectResult, AbstractResultData
import Tables

export pkimport, nca!

const LOG2 = log(2)

    include("types.jl")
    include("show.jl")
    include("import.jl")
    include("nca.jl")

end # module
