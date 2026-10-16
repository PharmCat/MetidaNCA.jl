################################################################################
# Acceptance criteria and divergence journal for the MetidaNCA validation report.
################################################################################

# ---------------------------------------------------------------- criteria ---
# A computed value x is accepted against a reference r when
#
#       |x - r| <= ACCEPT_RTOL * |r| + ACCEPT_ATOL
#
const ACCEPT_RTOL = 1.0e-4
const ACCEPT_ATOL = 0.0

const REQUIRED_SCENARIOS = [
    "PKG-TEST",
    "PK-LINT-EV-DT0",
    "PK-LULD-EV-DT025-TAU9",
    "PK-LINT-IV-TAU12",
    "PK-LOGT-EV-TAU12",
    "PK-LULDT-EV-TAU12",
    "URINE-LINT",
    "PD-LINT",
    "TRAP-EXP",
    "TRAP-ASC",
    "TRAP-ID",
    "TRAP-REIMP",
    "DERIVED-ID",
    "GRAPH-DATA",
    "CONCLUSION",
]

# The kind of evidence each scenario provides. This is what the coverage tables
# report, and it is the distinction a validation file has to make explicit:
#
#   external    comparison with an independent implementation;
#   internal    a value fixed by this package and carried forward. 
#   analytic    comparison with a closed-form value. 
#   derived     the parameter is a documented function of parameters that are
#               themselves qualified, and the report verifies that the function
#               holds.
#   operational the function runs and returns what it is supposed to return -
#               for graphics, that the figure carries the data that was
#               computed. Not a numeric qualification of a parameter.
const SCENARIO_KIND = Dict(
    "PKG-TEST"              => "operational",
    "PK-LINT-EV-DT0"        => "external",
    "PK-LULD-EV-DT025-TAU9" => "external",
    "PK-LINT-IV-TAU12"      => "external",
    "PK-LOGT-EV-TAU12"      => "external",
    "PK-LULDT-EV-TAU12"     => "internal",
    "URINE-LINT"            => "external",
    "PD-LINT"               => "external",
    "TRAP-EXP"              => "analytic",
    "TRAP-ASC"              => "analytic",
    "TRAP-ID"               => "analytic",
    "TRAP-REIMP"            => "analytic",
    "DERIVED-ID"            => "derived",
    "GRAPH-DATA"            => "operational",
    "CONCLUSION"            => "operational",
)

scenario_kind(s) = get(SCENARIO_KIND, s, "unclassified")

# Short titles for the per-scenario summary in the Conclusion.
const SCENARIO_TITLES = Dict(
    "PKG-TEST"              => "Package test suite",
    "PK-LINT-EV-DT0"        => "Plasma :lint, EV, dt 0, no tau, dose 100",
    "PK-LULD-EV-DT025-TAU9" => "Plasma :luld, EV, dt 0.25, tau 9, dose 100",
    "PK-LINT-IV-TAU12"      => "Plasma :lint, IV, dt 0, tau 12, dose 120",
    "PK-LOGT-EV-TAU12"      => "Plasma :logt, EV, dt 0, tau 12, dose 120",
    "PK-LULDT-EV-TAU12"     => "Plasma :luldt, EV, dt 0, tau 12, dose 120",
    "URINE-LINT"            => "Urine :lint, dose 100",
    "PD-LINT"               => "PD :lint, BL 3.0, TH 1.5",
    "TRAP-EXP"              => "Integration rules on a mono-exponential profile",
    "TRAP-ASC"              => "Integration rules on a strictly ascending profile",
    "TRAP-ID"               => "Switching-rule identities of :luldt",
    "TRAP-REIMP"            => "Integration rules against a direct reimplementation",
    "DERIVED-ID"            => "Derived parameters against their defining identities",
    "GRAPH-DATA"            => "Figures carry the computed data",
    "CONCLUSION"            => "Conclusion rendered",
)

# Every comparison, pass or fail.
const ACCEPT_LOG = joinpath(@__DIR__, "src", "validation_acceptance.tsv")
# Divergences only - the journal.
const ACCEPT_JOURNAL = joinpath(@__DIR__, "src", "validation_divergences.tsv")

const ACCEPT_LOG_HEADER = ["scenario", "parameter", "subject", "value", "reference", "difference", "status"]
const ACCEPT_JOURNAL_HEADER = ["scenario", "parameter", "subject", "value", "reference", "difference", "relative_error", "criterion"]

# How many divergences are listed in the PDF before the table is truncated and
# the reader is sent to the journal file.
const ACCEPT_JOURNAL_LIMIT = 200

# ------------------------------------------------------------------ logging ---

"""
    acceptance_init()

Remove any log and journal left over from a previous build and start new ones.
Must be called from `validation.jl` BEFORE `weave`, so that artefacts of a
previous build cannot be mistaken for the result of this one.
"""
function acceptance_init()
    mkpath(dirname(ACCEPT_LOG))
    open(ACCEPT_LOG, "w") do io
        println(io, join(ACCEPT_LOG_HEADER, '\t'))
    end
    rm(ACCEPT_JOURNAL; force = true)
    return ACCEPT_LOG
end

function _accept(value, reference)
    (isnan(value) && isnan(reference)) && return true
    (isnan(value) || isnan(reference)) && return false
    return isapprox(value, reference; rtol = ACCEPT_RTOL, atol = ACCEPT_ATOL)
end

function _log_row(scenario, parameter, subject, value, reference, difference, status)
    open(ACCEPT_LOG, "a") do io
        println(io, join((scenario, parameter, subject, value, reference, difference, status), '\t'))
    end
    return status
end

"""
    acceptance_record(scenario, parameter, subject, value, reference) -> "PASS"/"FAIL"

Compare one value against its reference, append the outcome to the log and
return the verdict for display in the report table.
"""
function acceptance_record(scenario, parameter, subject, value, reference)
    status = _accept(value, reference) ? "PASS" : "FAIL"
    _log_row(scenario, string(parameter), string(subject), value, reference, value - reference, status)
    return status
end

"""
    acceptance_mark(scenario, ok::Bool, note = "")

Record a non-numeric outcome: the test suite run, or the Conclusion chunk
confirming that it rendered.
"""
function acceptance_mark(scenario, ok::Bool, note = "")
    _log_row(scenario, note, "-", "-", "-", "-", ok ? "PASS" : "FAIL")
    return ok
end

function _read_log()
    isfile(ACCEPT_LOG) || return NamedTuple[]
    rows = NamedTuple[]
    for (i, line) in enumerate(eachline(ACCEPT_LOG))
        i == 1 && continue                       # header
        isempty(strip(line)) && continue
        f = split(line, '\t')
        length(f) == length(ACCEPT_LOG_HEADER) || continue
        push!(rows, (scenario = f[1], parameter = f[2], subject = f[3],
                     value = f[4], reference = f[5], difference = f[6], status = f[7]))
    end
    return rows
end

# A numeric comparison carries a value
_is_check(r) = r.value != "-"

# Bookkeeping marks the report makes about itself. 
const SELF_MARKS = ["CONCLUSION"]

# --------------------------------------------------------------- inspection ---
# Everything below reads the log. None of it depends on DataFrames, so it can
# be called from `validation.jl` as well as from the report.
#
# `inreport = true` is the view used by the Conclusion while it renders itself;
# `inreport = false` (the default) is the full view used by the build gate.

_rows(inreport) = inreport ? filter(r -> r.scenario ∉ SELF_MARKS, _read_log()) : _read_log()
_required(inreport) = inreport ? setdiff(REQUIRED_SCENARIOS, SELF_MARKS) : REQUIRED_SCENARIOS

"""
    acceptance_scenarios(; inreport = false) -> Vector{String}

Scenarios that reported at least one row, in the order of `REQUIRED_SCENARIOS`
first and then any additional ones.
"""
function acceptance_scenarios(; inreport = false)
    reported = unique(getfield.(_rows(inreport), :scenario))
    known = filter(s -> s in reported, _required(inreport))
    return vcat(known, sort(setdiff(reported, _required(inreport))))
end

"""
    acceptance_missing(; inreport = false) -> Vector{String}

Required scenarios absent from the log, i.e. whose chunk errored out or was
removed from the report.
"""
acceptance_missing(; inreport = false) =
    setdiff(_required(inreport), unique(getfield.(_rows(inreport), :scenario)))

"""
    acceptance_counts(; inreport = false) -> NamedTuple

`(checks, failed_checks, failed_marks, scenarios, missing)`.
"""
function acceptance_counts(; inreport = false)
    rows = _rows(inreport)
    failed = filter(r -> r.status == "FAIL", rows)
    return (checks        = count(_is_check, rows),
            failed_checks = count(_is_check, failed),
            failed_marks  = count(!_is_check, failed),
            scenarios     = length(unique(getfield.(rows, :scenario))),
            missing       = length(acceptance_missing(; inreport = inreport)))
end

"""
    acceptance_ok(; inreport = false) -> Bool

`true` when the log exists, is not empty, every required scenario reported and
every recorded outcome is a PASS.
"""
function acceptance_ok(; inreport = false)
    isfile(ACCEPT_LOG) || return false
    rows = _rows(inreport)
    isempty(rows) && return false
    isempty(acceptance_missing(; inreport = inreport)) || return false
    return !any(r -> r.status == "FAIL", rows)
end

acceptance_status(; inreport = false) = acceptance_ok(; inreport = inreport) ? "PASS" : "FAIL"

_accnum(s) = s == "-" ? nothing : tryparse(Float64, s)

# Relative error |v - r| / |r|, as a string for display. Returns "n/a" when the
# reference is zero (the criterion then reduces to the absolute tolerance) or
# when either side is not a number.
function _relerr_str(value, reference)
    v, r = _accnum(value), _accnum(reference)
    (isnothing(v) || isnothing(r)) && return "n/a"
    (isnan(v) || isnan(r)) && return "NaN"
    iszero(r) && return "n/a"
    return string(round(abs(v - r) / abs(r), sigdigits = 3))
end

# The criterion the comparison was judged against, spelled out per row
function _criterion_str(reference)
    r = _accnum(reference)
    (isnothing(r) || isnan(r)) && return "-"
    return string(round(ACCEPT_RTOL * abs(r) + ACCEPT_ATOL, sigdigits = 3))
end

"""
    acceptance_divergences(; inreport = false) -> Vector{NamedTuple}

Every recorded outcome that is outside the acceptance criteria, with its
relative error and the criterion it was judged against. This is the journal.
"""
function acceptance_divergences(; inreport = false)
    out = NamedTuple[]
    for r in filter(x -> x.status == "FAIL", _rows(inreport))
        push!(out, (scenario = r.scenario, parameter = r.parameter, subject = r.subject,
                    value = r.value, reference = r.reference, difference = r.difference,
                    relative_error = _relerr_str(r.value, r.reference),
                    criterion = _criterion_str(r.reference)))
    end
    return out
end

"""
    acceptance_write_journal() -> String

Write the divergence journal to `ACCEPT_JOURNAL` and return its path. The file
is written even when there are no divergences - an empty journal with a header
is a positive statement, a missing journal is not.

Called from the Conclusion chunk, and again from `validation.jl` in case that
chunk did not run.
"""
function acceptance_write_journal()
    mkpath(dirname(ACCEPT_JOURNAL))
    open(ACCEPT_JOURNAL, "w") do io
        println(io, join(ACCEPT_JOURNAL_HEADER, '\t'))
        for d in acceptance_divergences()
            println(io, join((d.scenario, d.parameter, d.subject, d.value, d.reference,
                              d.difference, d.relative_error, d.criterion), '\t'))
        end
        for s in acceptance_missing()
            println(io, join((s, "scenario did not report", "-", "-", "-", "-", "n/a", "-"), '\t'))
        end
    end
    return ACCEPT_JOURNAL
end

# ------------------------------------------------------------------- tables ---
# These build the tables shown in the report and need DataFrames.

const _REFERENCE_CACHE = Dict{String, Any}()

"""
    scenario_references(path, scenario, subjects) -> (Dict{Symbol,Vector{Float64}}, Vector{Symbol})

Every reference value the external set holds for one scenario, in the order of
`subjects`, together with the sorted list of parameters.
"""
function scenario_references(path, scenario, subjects)
    df = get!(() -> DataFrame(CSV.File(path)), _REFERENCE_CACHE, path)
    byparsubj = Dict{Tuple{Symbol, Int}, Float64}()
    pars = Symbol[]
    for r in eachrow(df)
        String(r.scenario) == scenario || continue
        ismissing(r.metida_parameter) && continue
        isempty(strip(String(r.metida_parameter))) && continue
        p = Symbol(strip(String(r.metida_parameter)))
        byparsubj[(p, Int(r.subject))] = float(r.value)
        p in pars || push!(pars, p)
    end
    isempty(pars) && error("The reference set holds nothing for scenario \"$scenario\"")
    out = Dict(p => [get(byparsubj, (p, Int(s)), NaN) for s in subjects] for p in pars)
    return out, sort!(pars)
end

"""
    mdtable(df) -> String

Render a data frame as the Markdown pipe table the report embeds, with every
`|` inside a cell or a column name escaped.
"""
function mdtable(df)
    esc(x) = replace(string(x), "|" => "\\|")
    out = DataFrame()
    found = String[]
    for n in names(df)
        occursin('|', n) && push!(found, "column name \"$n\"")
        col = df[!, n]
        for v in col
            occursin('|', string(v)) && push!(found, "a cell of column \"$n\"")
        end
        out[!, esc(n)] = esc.(col)
    end
    # `where` is a reserved word, so the key has to be something else.
    isempty(found) ||
        @warn "A Markdown table cell contained \"|\" and was escaped" places = unique(found)
    return sprint(show, latexify(out; latex = false))
end

"""
    acceptance_table(scenario, dsnca, refdict, list, subjects) -> DataFrame

Per-subject comparison table for a PK scenario. Values are displayed rounded to
six significant digits; the acceptance check itself is made on the unrounded
values.
"""
function acceptance_table(scenario, dsnca, refdict, list, subjects)
    dft = DataFrame("Parameter" => String[], "Subject" => String[], "Value" => Any[],
                    "Reference" => Any[], "Difference" => Any[], "Status" => String[])
    for l in list
        haskey(refdict, l) || continue
        # A parameter the package does not produce for this configuration is not
        # a divergence. `Swing` is the case to keep in mind: it is undefined when
        # Cmin is zero, so it exists for some subjects of a scenario and not for
        # others, and the reference set is blank in exactly those places.
        first(dsnca).result isa AbstractDict && !haskey(first(dsnca).result, l) && continue
        vals = dsnca[:, l]
        refs = refdict[l]
        for i in eachindex(subjects)
            ismissing(vals[i]) && continue
            v = float(vals[i])
            r = float(refs[i])
            isnan(r) && !isnan(v) && continue      # no reference for this subject
            status = acceptance_record(scenario, l, Int(subjects[i]), v, r)
            push!(dft, (string(l), string(Int(subjects[i])), round(v, sigdigits = 6),
                        round(r, sigdigits = 6), round(v - r, sigdigits = 6), status))
        end
    end
    return dft
end

"""
    acceptance_table_single(scenario, result, refdict, list) -> DataFrame

Comparison table for a single-subject result (urine, pharmacodynamics).
"""
function acceptance_table_single(scenario, result, refdict, list; index = nothing)
    dft = DataFrame("Parameter" => String[], "Value" => Any[], "Reference" => Any[],
                    "Difference" => Any[], "Status" => String[])
    for k in list
        haskey(refdict, k) || continue
        res = isnothing(index) ? result : result[index]
        res.result isa AbstractDict && !haskey(res.result, k) && continue
        v = float(isnothing(index) ? result[k] : result[index, k])
        # scenario_references returns one value per subject; a single-subject
        # scenario has one.
        rv = refdict[k]
        r = float(rv isa AbstractVector ? first(rv) : rv)
        isnan(r) && !isnan(v) && continue
        status = acceptance_record(scenario, k, "-", v, r)
        push!(dft, (string(k), round(v, sigdigits = 6), round(r, sigdigits = 6),
                    round(v - r, sigdigits = 6), status))
    end
    return dft
end

"""
    acceptance_note(scenario) -> String

One-line summary printed under each results table. Names the diverging
parameters directly, so a divergence is visible at the table it came from and
not only in the Conclusion.
"""
function acceptance_note(scenario)
    rows = filter(r -> r.scenario == scenario, _read_log())
    n = count(_is_check, rows)
    bad = filter(r -> r.status == "FAIL", rows)
    isempty(bad) && return "Comparisons: $n; outside acceptance criteria: 0. Result: PASS."
    named = unique([isempty(r.parameter) ? "-" : r.parameter for r in bad])
    shown = join(Iterators.take(named, 8), ", ")
    length(named) > 8 && (shown *= ", ...")
    return "Comparisons: $n; outside acceptance criteria: $(length(bad)) " *
           "($shown). Result: FAIL - see the divergence journal in the Conclusion."
end

"""
    acceptance_summary_table() -> DataFrame

Per-scenario summary for the Conclusion: how many comparisons each scenario
made, how many fell outside the criteria, and its verdict. A required scenario
that produced no rows at all appears as NOT REPORTED rather than being absent.
"""
function acceptance_summary_table()
    rows = _rows(true)
    dft = DataFrame("Scenario" => String[], "Description" => String[],
                    "Checks" => Any[], "Failures" => Any[], "Status" => String[])
    for s in acceptance_scenarios(; inreport = true)
        srows = filter(r -> r.scenario == s, rows)
        nfail = count(r -> r.status == "FAIL", srows)
        push!(dft, (s, get(SCENARIO_TITLES, s, "-"), count(_is_check, srows), nfail,
                    nfail == 0 ? "PASS" : "FAIL"))
    end
    for s in acceptance_missing(; inreport = true)
        push!(dft, (s, get(SCENARIO_TITLES, s, "-"), 0, "-", "NOT REPORTED"))
    end
    return dft
end

# Full-precision values go to the TSV; the PDF gets six significant digits, so
# that the journal table stays inside the page width.
function _fmt(s)
    v = _accnum(s)
    isnothing(v) && return s
    isnan(v) && return "NaN"
    return string(round(v, sigdigits = 6))
end

"""
    acceptance_divergence_table(; limit = ACCEPT_JOURNAL_LIMIT) -> DataFrame

The divergence journal as shown in the Conclusion. Empty when nothing diverged;
truncated at `limit` rows, the untruncated journal being the TSV file. The
per-row criterion is kept in the TSV only - in the report the criterion is the
single `rtol` stated under Acceptance criteria.
"""
function acceptance_divergence_table(; limit = ACCEPT_JOURNAL_LIMIT)
    dft = DataFrame("Scenario" => String[], "Parameter" => String[], "Subject" => String[],
                    "Value" => String[], "Reference" => String[], "Difference" => String[],
                    "Rel. error" => String[])
    for d in Iterators.take(acceptance_divergences(; inreport = true), limit)
        push!(dft, (d.scenario, isempty(d.parameter) ? "-" : d.parameter, d.subject,
                    _fmt(d.value), _fmt(d.reference), _fmt(d.difference), d.relative_error))
    end
    for s in acceptance_missing(; inreport = true)
        push!(dft, (s, "scenario did not report", "-", "-", "-", "-", "n/a"))
    end
    return dft
end

"""
    acceptance_journal_block() -> String

The whole "divergence journal" block of the Conclusion, ready to be printed
into the report: either the rendered journal table with its caption, or the
statement that nothing diverged. Emitted from one function because a Markdown
caption line cannot be made conditional.
"""
function acceptance_journal_block()
    dft = acceptance_divergence_table()
    total = length(acceptance_divergences(; inreport = true)) +
            length(acceptance_missing(; inreport = true))
    io = IOBuffer()
    if total == 0
        println(io, "No divergences. Every recorded comparison is within the acceptance ",
                    "criteria, the package test suite passed and every required scenario ",
                    "reported. The journal file `docs/src/validation_divergences.tsv` is ",
                    "written with its header and no rows.")
        return String(take!(io))
    end
    print(io, mdtable(dft))
    println(io)
    println(io)
    println(io, "Table: Divergence journal - every recorded outcome that failed acceptance")
    println(io)
    if total > ACCEPT_JOURNAL_LIMIT
        println(io, "Listed above: ", ACCEPT_JOURNAL_LIMIT, " of ", total,
                    " divergences. The remaining ", total - ACCEPT_JOURNAL_LIMIT,
                    " are in `docs/src/validation_divergences.tsv`.")
    else
        println(io, "Divergences: ", total, ", all listed above.")
    end
    println(io, "Values are shown to six significant digits; full precision and the ",
                "per-row criterion are in `docs/src/validation_divergences.tsv`.")
    return String(take!(io))
end

"""
    acceptance_conclusion_text() -> String

The verdict paragraph of the Conclusion, stating the outcome, the criteria and
where to look next. Written for a reader who has the report and nothing else.
"""
function acceptance_conclusion_text()
    c = acceptance_counts(; inreport = true)
    ok = acceptance_ok(; inreport = true)
    io = IOBuffer()
    if ok
        println(io, "**Result: PASS.** ", c.checks, " comparisons over ",
                c.scenarios, " scenarios, all within the acceptance criteria ",
                "(rtol = ", ACCEPT_RTOL, ", atol = ", ACCEPT_ATOL, "). ",
                "The package test suite passed and every required scenario reported. ",
                "The divergence journal below is empty.")
    else
        println(io, "**Result: FAIL.** ", c.checks, " comparisons over ",
                c.scenarios, " scenarios; ", c.failed_checks,
                " outside the acceptance criteria (rtol = ", ACCEPT_RTOL,
                ", atol = ", ACCEPT_ATOL, ").")
        c.failed_marks > 0 && println(io, "\nRecorded outcomes that failed: ",
                                      c.failed_marks,
                                      " (see the per-scenario summary - a failed test suite is reported here).")
        c.missing > 0 && println(io, "\nRequired scenarios that did not report: ",
                                 join(acceptance_missing(; inreport = true), ", "),
                                 ". A scenario reports nothing when its code chunk errored out, ",
                                 "so its results table is missing from this report as well.")
        println(io, "\nEvery divergence is listed in the journal below and in ",
                "`docs/src/validation_divergences.tsv`. ",
                "This report is not evidence of conformance in its present state.")
    end
    return String(take!(io))
end

# ------------------------------------------------------------------ coverage ---
# The coverage tables are derived from the same log as the verdict

"""
    acceptance_parameter_kinds() -> Dict{String, Vector{String}}

For every parameter that any numeric comparison in the log refers to, the kinds
of evidence that compared it (see `SCENARIO_KIND`).
"""
function acceptance_parameter_kinds()
    out = Dict{String, Set{String}}()
    for r in _rows(true)
        _is_check(r) || continue
        # An operational block (graphics) checks properties of a figure, not the
        # value of a parameter
        scenario_kind(r.scenario) == "operational" && continue
        push!(get!(out, r.parameter, Set{String}()), scenario_kind(r.scenario))
    end
    return Dict(k => sort!(collect(v)) for (k, v) in out)
end

"""
    acceptance_parameter_counts() -> Dict{String, Int}

Number of numeric comparisons per parameter.
"""
function acceptance_parameter_counts()
    out = Dict{String, Int}()
    for r in _rows(true)
        _is_check(r) || continue
        scenario_kind(r.scenario) == "operational" && continue
        out[r.parameter] = get(out, r.parameter, 0) + 1
    end
    return out
end

_mark(b) = b ? "yes" : "-"

"""
    acceptance_coverage_table(producible) -> DataFrame

Parameter-level coverage. 

Status is:

  * `qualified`       - compared with an external or an analytic reference;
  * `regression only` - compared only with a value this package itself fixed;
  * `NOT COVERED`     - not compared with anything in this report.
"""
function acceptance_coverage_table(producible)
    kinds  = acceptance_parameter_kinds()
    counts = acceptance_parameter_counts()
    names  = sort!(unique(vcat(String.(producible), collect(keys(kinds)))))
    dft = DataFrame("Parameter" => String[], "External" => String[], "Analytic" => String[],
                    "Derived" => String[], "Internal" => String[],
                    "Checks" => Int[], "Status" => String[])
    for n in names
        k = get(kinds, n, String[])
        ext, ana = "external" in k, "analytic" in k
        der, int = "derived" in k, "internal" in k
        status = (ext || ana) ? "qualified" :
                 der          ? "derived" :
                 int          ? "regression only" : "NOT COVERED"
        push!(dft, (n, _mark(ext), _mark(ana), _mark(der), _mark(int),
                    get(counts, n, 0), status))
    end
    return dft
end

"""
    acceptance_coverage_summary(producible) -> String

One paragraph under the parameter-coverage table: how many parameters are
qualified, how many rest on a regression baseline only, how many are not
covered, and which ones those are.
"""
function acceptance_coverage_summary(producible)
    dft = acceptance_coverage_table(producible)
    nq = count(==("qualified"), dft.Status)
    nd = count(==("derived"), dft.Status)
    nr = count(==("regression only"), dft.Status)
    nn = count(==("NOT COVERED"), dft.Status)
    io = IOBuffer()
    print(io, "Parameters produced or checked during this build: ", nrow(dft), ". ",
              "Qualified against an external or analytic reference: ", nq, ". ",
              "Established through a defining identity only: ", nd, ". ",
              "Resting on an internal regression baseline only: ", nr, ". ",
              "Not covered by any comparison in this report: ", nn, ".")
    if nd > 0
        print(io, "\n\nIdentity only: ",
                  join(dft.Parameter[dft.Status .== "derived"], ", "),
                  ". No other program reports these, so there is nothing to compare them ",
                  "with; what the report establishes is that each is what its definition ",
                  "says it is, computed from parameters that are themselves qualified.")
    end
    if nr > 0
        print(io, "\n\nRegression baseline only: ",
                  join(dft.Parameter[dft.Status .== "regression only"], ", "),
                  ". A baseline shows that behaviour has not changed; it is not evidence of ",
                  "correctness.")
    end
    if nn > 0
        print(io, "\n\nNot covered: ",
                  join(dft.Parameter[dft.Status .== "NOT COVERED"], ", "),
                  ". These are computed by the package and reported to the user, but no ",
                  "comparison in this report establishes their values. Where they are ",
                  "exercised by the package test suite, that is an operational check, not a ",
                  "numeric qualification.")
    end
    return String(take!(io))
end

# ---- functionality coverage -------------------------------------------------
# Declared scope of the validation, one row per capability of the package. The
# status column is derived from the log, so an item claiming coverage by a
# scenario that did not report is shown as such rather than taken on trust.
#
# `scenarios` are the scenario codes that exercise the item in THIS report;
# `elsewhere` says where the item is covered when this report does not cover it.
const COVERAGE_ITEMS = [
    (item = "pkimport - plasma data import",
     scenarios = ["PK-LINT-EV-DT0", "PK-LULD-EV-DT025-TAU9"], elsewhere = "",
     note = "Every plasma scenario imports through it"),
    (item = "upkimport - urine data import",
     scenarios = ["URINE-LINT"], elsewhere = "",
     note = ""),
    (item = "pdimport - pharmacodynamics import",
     scenarios = ["PD-LINT"], elsewhere = "",
     note = ""),
    (item = "calcm = :lint",
     scenarios = ["PK-LINT-EV-DT0", "PK-LINT-IV-TAU12", "TRAP-ASC", "TRAP-REIMP"], elsewhere = "",
     note = ""),
    (item = "calcm = :logt",
     scenarios = ["PK-LOGT-EV-TAU12", "TRAP-EXP", "TRAP-REIMP"], elsewhere = "",
     note = ""),
    (item = "calcm = :luld",
     scenarios = ["PK-LULD-EV-DT025-TAU9", "TRAP-EXP", "TRAP-REIMP"], elsewhere = "",
     note = ""),
    (item = "calcm = :luldt",
     scenarios = ["PK-LULDT-EV-TAU12", "TRAP-EXP", "TRAP-ID", "TRAP-REIMP"], elsewhere = "",
     note = "No external counterpart exists; qualified analytically, see the integration-rules chapter"),
    (item = "adm = :ev",
     scenarios = ["PK-LINT-EV-DT0", "PK-LOGT-EV-TAU12"], elsewhere = "",
     note = ""),
    (item = "adm = :iv",
     scenarios = ["PK-LINT-IV-TAU12", "TRAP-EXP"], elsewhere = "",
     note = ""),
    (item = "Infusion (DoseTime rate, IT)",
     scenarios = String[], elsewhere = "not implemented",
     note = "DoseTime.rate is stored and printed but enters no computation; MRT carries no IT/2 correction"),
    (item = "Steady state (tau): AUCtau, Ctau, Cavg, Fluc, Swing, Accind",
     scenarios = ["PK-LULD-EV-DT025-TAU9", "PK-LINT-IV-TAU12", "PK-LOGT-EV-TAU12"], elsewhere = "",
     note = ""),
    (item = "Dose time other than zero",
     scenarios = ["PK-LULD-EV-DT025-TAU9"], elsewhere = "",
     note = ""),
    (item = "Terminal phase: automatic interval search",
     scenarios = ["PK-LINT-EV-DT0", "TRAP-EXP"], elsewhere = "",
     note = "Kel, HL, Rsq, ARsq, NpLZ compared in every plasma scenario"),
    (item = "Terminal phase: user-set interval (ElimRange, setkelrange!, kelexcl)",
     scenarios = String[], elsewhere = "package test suite",
     note = "Exercised by the test suite; no reference comparison in this report"),
    (item = "Partial AUC (partials, prtext)",
     scenarios = String[], elsewhere = "package test suite",
     note = "Exercised by the test suite; no reference comparison in this report"),
    (item = "BLOQ and missing values (LimitRule, applylimitrule!)",
     scenarios = String[], elsewhere = "package test suite",
     note = "Exercised by the test suite; no reference comparison in this report"),
    (item = "Sparse data (auc_sparse)",
     scenarios = String[], elsewhere = "package test suite",
     note = "Formula checked by the test suite; standard errors are not implemented"),
    (item = "Multiple dose times (usedose)",
     scenarios = String[], elsewhere = "package test suite",
     note = ""),
    (item = "Multiple observation columns",
     scenarios = String[], elsewhere = "package test suite",
     note = ""),
    (item = "Units (Unitful quantities)",
     scenarios = String[], elsewhere = "package test suite",
     note = ""),
    (item = "Export: DataFrame, metida_table",
     scenarios = String[], elsewhere = "package test suite",
     note = "Row and subject order checked by the test suite"),
    (item = "Graphics: individual profile, linear and semi-logarithmic",
     scenarios = ["GRAPH-DATA"], elsewhere = "",
     note = "Figures reproduced below; the data they carry is compared with the computed data"),
    (item = "Graphics: data set overlay, mean profile",
     scenarios = ["GRAPH-DATA"], elsewhere = "",
     note = ""),
    (item = "Graphics: terminal-phase interval and excluded points",
     scenarios = ["GRAPH-DATA"], elsewhere = "",
     note = ""),
    (item = "Graphics: pharmacodynamics with baseline and threshold",
     scenarios = ["GRAPH-DATA"], elsewhere = "",
     note = ""),
    (item = "Graphics: urine excretion rate against interval midpoints",
     scenarios = ["GRAPH-DATA"], elsewhere = "",
     note = "The figure is checked against midtime and excretionrate, and its maximum against TmaxRate and MaxRate"),
    (item = "Minrate (urine)",
     scenarios = String[], elsewhere = "cannot be qualified",
     note = "Identically zero: nca! prepends a zero to the rate vector for the curve origin and takes the minimum after that (nca.jl, urine branch)"),
    (item = "Derived parameters (Cllast, Vzlast, Vsstau, LZ, Rsqn, AUCBTW)",
     scenarios = ["DERIVED-ID"], elsewhere = "",
     note = "No external counterpart exists; each is checked against the identity that defines it"),
    (item = "Descriptive statistics (n, mean, SD, CV, geometric mean and CV, median, range)",
     scenarios = String[], elsewhere = "not implemented",
     note = "The package computes individual parameters only"),
    (item = "Subsetting: timefilter, subset",
     scenarios = String[], elsewhere = "package test suite",
     note = ""),
]

"""
    acceptance_functionality_table() -> DataFrame

Declared scope of the validation, one row per capability, with the status of
each derived from the log: which of its scenarios reported, and with what kind
of evidence.
"""
function acceptance_functionality_table()
    reported = Set(unique(getfield.(_rows(true), :scenario)))
    dft = DataFrame("Functionality" => String[], "Scenarios" => String[],
                    "Evidence" => String[], "Note" => String[])
    for it in COVERAGE_ITEMS
        present = filter(s -> s in reported, it.scenarios)
        kinds = sort!(unique(scenario_kind.(present)))
        evidence =
            !isempty(present)        ? join(kinds, ", ") :
            !isempty(it.scenarios)   ? "DECLARED BUT DID NOT REPORT" :
            isempty(it.elsewhere)    ? "NOT VALIDATED" : it.elsewhere
        push!(dft, (it.item, isempty(present) ? "-" : join(present, ", "), evidence, it.note))
    end
    return dft
end

"""
    acceptance_functionality_summary() -> String

One paragraph under the functionality table, counting the capabilities covered
by this report, those covered elsewhere, and those covered by nothing.
"""
function acceptance_functionality_summary()
    dft = acceptance_functionality_table()
    inreport = count(r -> r != "-", dft.Scenarios)
    suite    = count(==("package test suite"), dft.Evidence)
    absent   = count(x -> x in ("not implemented", "not supported", "NOT VALIDATED",
                               "cannot be qualified"), dft.Evidence)
    declared = count(==("DECLARED BUT DID NOT REPORT"), dft.Evidence)
    io = IOBuffer()
    print(io, "Capabilities in the declared scope: ", nrow(dft),
              ". Qualified in this report: ", inreport,
              ". Exercised by the package test suite only: ", suite,
              ". Not implemented, not supported, or not validated anywhere: ", absent, ".")
    declared > 0 && print(io, " Declared as covered by a scenario that did not report: ",
                          declared, " - see the Conclusion.")
    print(io, "\n\nA capability exercised by the test suite is covered operationally: the code runs ",
              "and its result is self-consistent. That is weaker than a comparison with an ",
              "external or analytic reference, and this table keeps the two apart on purpose.")
    return String(take!(io))
end

# ------------------------------------------------------------------ verdict ---

"""
    acceptance_verdict(; io = stdout, strict = true) -> Bool

Build gate. Prints the verdict and the full divergence journal to `io`, writes
the journal file, and returns `true` when everything is within criteria.

With `strict = true` (the default, and what CI uses) a negative verdict raises,
so the documentation build fails loudly. With `strict = false` the build
completes and the report is left for inspection - the verdict is in the report
itself either way.

Called from `validation.jl` after `weave`, i.e. outside Weave, because Weave
renders an error raised inside a chunk into the document and still reports
success.
"""
function acceptance_verdict(; io = stdout, strict = true)
    msg = IOBuffer()

    if !isfile(ACCEPT_LOG)
        println(msg, "VALIDATION FAILED: acceptance log not found at $ACCEPT_LOG.")
        println(msg, "The report was not generated through validation.jl, or every chunk failed.")
    elseif isempty(_read_log())
        # The log is created by acceptance_init() before the build and appended
        # to by the chunks. A header and nothing else means every chunk that
        # should have recorded something failed, which in practice means the
        # preamble failed and took the whole document with it.
        println(msg, "VALIDATION FAILED: acceptance log at $ACCEPT_LOG has its header and no rows -")
        println(msg, "not one comparison was recorded, although the log was created, so the build did start.")
        println(msg, "")
        println(msg, "Every chunk that records a comparison failed. The usual cause is the preamble:")
        println(msg, "the document generator aborts a chunk at its first failing statement and carries")
        println(msg, "on with the next one, so a package missing from the active environment, or an")
        println(msg, "error in `include(\"acceptance.jl\")`, leaves every later chunk without the")
        println(msg, "recording functions. The report is then built out of error messages.")
        println(msg, "")
        println(msg, "Open the generated report and find the FIRST line beginning with \"Error:\" -")
        println(msg, "that is the failure; the rest are its consequences. Then:")
        println(msg, "  * a missing package: instantiate the documentation environment and build with it")
        println(msg, "        julia --project=docs -e 'using Pkg; Pkg.develop(path=\".\"); Pkg.instantiate()'")
        println(msg, "        julia --project=docs docs/make.jl")
        println(msg, "  * anything else: the message names the scenario or the function.")
        println(msg, "")
        println(msg, "Scenarios that should have reported: ", join(REQUIRED_SCENARIOS, ", "), ".")
    else
        journal = acceptance_write_journal()
        c = acceptance_counts()
        miss = acceptance_missing()
        divs = acceptance_divergences()

        # Scenario counts are reported without the report's own bookkeeping
        # marks, so that the console line and the report agree.
        nscen = acceptance_counts(; inreport = true).scenarios

        if acceptance_ok()
            println(io, "VALIDATION PASSED: ", c.checks, " checks over ", nscen,
                    " scenarios, all within criteria (rtol = ", ACCEPT_RTOL, ").")
            println(io, "Full log: ", ACCEPT_LOG)
            println(io, "Divergence journal (empty): ", journal)
            return true
        end

        println(msg, "VALIDATION FAILED.")
        println(msg, "Acceptance criteria: rtol = $ACCEPT_RTOL, atol = $ACCEPT_ATOL.")
        println(msg, "Checks performed: $(c.checks) over $nscen scenarios.")
        if c.failed_marks > 0
            for r in filter(r -> r.status == "FAIL" && !_is_check(r), _read_log())
                println(msg, "Recorded outcome failed: ", r.scenario,
                        isempty(r.parameter) ? "" : " (" * r.parameter * ")")
            end
        end
        if !isempty(miss)
            println(msg, "Scenarios that did not report (chunk errored or was removed): ",
                    join(miss, ", "))
        end
        nchecks = count(d -> d.value != "-", divs)
        if nchecks > 0
            println(msg, "Comparisons outside acceptance criteria: $nchecks")
            println(msg, "  scenario / parameter / subject : value vs reference (rel. error, criterion)")
            for d in Iterators.take(filter(d -> d.value != "-", divs), ACCEPT_JOURNAL_LIMIT)
                println(msg, "  ", d.scenario, " / ", d.parameter, " / subject ", d.subject,
                        " : ", d.value, " vs ", d.reference,
                        " (rel. err ", d.relative_error, ", criterion ", d.criterion, ")")
            end
            nchecks > ACCEPT_JOURNAL_LIMIT &&
                println(msg, "  ... and $(nchecks - ACCEPT_JOURNAL_LIMIT) more")
        end
        println(msg, "Divergence journal: $journal")
        println(msg, "Full log: $ACCEPT_LOG")
        println(msg, "The report itself carries the same journal in its Conclusion section.")
    end

    text = String(take!(msg))
    if strict
        error(text)
    else
        print(io, text)
        println(io, "VALIDATION_STRICT=false - the build is not aborted; inspect the report and the journal above.")
        return false
    end
end
