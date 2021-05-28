# Metida
# Copyright © 2019-2020 Vladimir Arnautov aka PharmCat <mail@pharmcat.net>

module MetidaNCA

import Base: ht_keyindex, length
import MetidaBase: AbstractIdData, AbstractSubject, DataSet
using Tables

export pkimport

    include("types.jl")
    include("show.jl")
    include("import.jl")
    include("nca.jl")

end # module
