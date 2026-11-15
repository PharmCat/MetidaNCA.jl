################################################################################
# Loader of the combined reference value set.
################################################################################

using CSV, DataFrames

const REFERENCE_DIR = joinpath(@__DIR__, "csv")
const REFERENCE_VALUES_CSV    = joinpath(REFERENCE_DIR, "reference_values.csv")
const REFERENCE_SCENARIOS_CSV = joinpath(REFERENCE_DIR, "reference_scenarios.csv")

_blank(x) = ismissing(x) || (x isa AbstractString && isempty(strip(x)))
_num(x)   = _blank(x) ? nothing : (x isa Number ? float(x) : parse(Float64, strip(String(x))))
_sym(x)   = _blank(x) ? nothing : Symbol(strip(String(x)))


function reference_scenarios(path::AbstractString = REFERENCE_SCENARIOS_CSV)
    isfile(path) || error("Metadata file not found: $path")
    df = CSV.File(path) |> DataFrame
    out = Dict{String, NamedTuple}()
    for r in eachrow(df)
        origin = _blank(r.origin) ? "external" : lowercase(strip(String(r.origin)))
        origin in ("external", "internal") ||
            error("Unknown origin \"$origin\" for scenario $(r.scenario); expected external or internal")
        out[String(r.scenario)] = (
            origin   = origin,
            dataset  = String(r.dataset),
            dose     = _num(r.dose),
            dosetime = _num(r.dosetime),
            tau      = _num(r.tau),
            calcm    = _sym(r.calcm),
            adm      = _sym(r.adm),
            bl       = _num(r.bl),
            th       = _num(r.th),
            sources  = _blank(r.sources) ? String[] : String.(strip.(split(String(r.sources), '|'))),
        )
    end
    return out
end

#=
Values by scenario. By default the keys are MetidaNCA parameter names; rows
without a match are skipped.
Vectors are ordered by subject number.
=#
function reference_values(path::AbstractString = REFERENCE_VALUES_CSV; phoenix::Bool = false)
    isfile(path) || error("Data file not found: $path")
    df = CSV.File(path) |> DataFrame
    namecol = phoenix ? :phoenix_parameter : :metida_parameter
    # full list of the scenario subjects, from all its rows
    allsubj = Dict{String, Vector{Int}}()
    for r in eachrow(df)
        push!(get!(allsubj, String(r.scenario), Int[]), Int(r.subject))
    end
    for (k, v) in allsubj
        allsubj[k] = sort!(unique!(v))
    end
    raw = Dict{String, Dict{Symbol, Dict{Int, Float64}}}()
    for r in eachrow(df)
        _blank(r[namecol]) && continue
        sc  = String(r.scenario)
        par = Symbol(strip(String(r[namecol])))
        get!(get!(raw, sc, Dict{Symbol, Dict{Int, Float64}}()), par, Dict{Int, Float64}())[Int(r.subject)] = float(r.value)
    end
    out = Dict{String, Dict{Symbol, Vector{Float64}}}()
    for (sc, pars) in raw
        subj = allsubj[sc]
        d = Dict{Symbol, Vector{Float64}}()
        for (par, bysubj) in pars
            # a parameter can have no value for some subjects (for example, Swing
            # is not defined at Cmin = 0); such positions are filled with NaN so
            # that the vector length always equals the number of subjects
            d[par] = [get(bysubj, i, NaN) for i in subj]
        end
        out[sc] = d
    end
    return out
end

"""
    reference_subjects(path = REFERENCE_VALUES_CSV) -> Dict{String, Vector{Int}}
Subject numbers of each scenario, in the order of the values in the vectors
returned by `reference_values`.
"""
function reference_subjects(path::AbstractString = REFERENCE_VALUES_CSV)
    df = CSV.File(path) |> DataFrame
    out = Dict{String, Vector{Int}}()
    for r in eachrow(df)
        push!(get!(out, String(r.scenario), Int[]), Int(r.subject))
    end
    for (k, v) in out
        out[k] = sort!(unique!(v))
    end
    return out
end

"""
    check_references() -> Bool
"""
function check_references()
    sc = reference_scenarios()
    rv = reference_values()
    rp = reference_values(; phoenix = true)
    subj = reference_subjects()
    problems = String[]
    # check against the scenario list of the values file: rows with origin = internal
    present = keys(subj)
    for k in setdiff(present, keys(sc))
        push!(problems, "No metadata for scenario: $k")
    end
    for k in setdiff(keys(sc), present)
        push!(problems, "No values for scenario: $k")
    end
    gaps = String[]
    for (k, d) in rv
        n = length(subj[k])
        for (par, v) in d
            length(v) == n || push!(problems, "$k / $par: $(length(v)) values for $n subject(s)")
            m = count(isnan, v)
            m == 0 || push!(gaps, "$k / $par: missing value(s) for $m subject(s)")
        end
    end
    isempty(gaps) || @info "Missing values replaced by NaN" Missings = gaps

    isempty(problems) || error("References set not agreed:\n  " * join(problems, "\n  "))
    # correspondence, rp slice - rows with origin = internal
    nval = nrow(CSV.File(REFERENCE_VALUES_CSV) |> DataFrame)
    next = count(x -> x.origin == "external", values(sc))
    nint = count(x -> x.origin == "internal", values(sc))
    @info "Refs loaded" Jobs = length(sc) Values = nval External = next Internal = nint
    return true
end
