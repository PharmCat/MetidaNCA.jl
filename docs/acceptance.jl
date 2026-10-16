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

# Scenarios that MUST be present in the log for the report to be considered
# complete. A chunk that errors out produces no rows, and its absence here is
# what turns a silently missing table into a failed gate. "CONCLUSION" is the
# Conclusion chunk marking itself once it has rendered the journal, so that a
# report without a verdict is caught too.
const REQUIRED_SCENARIOS = [
    "PKG-TEST",
    "PK-LINT-EV-DT0",
    "PK-LULD-EV-DT025-TAU9",
    "PK-LINT-IV-TAU12",
    "PK-LOGT-EV-TAU12",
    "URINE-LINT",
    "PD-LINT",
    "CONCLUSION",
]

# Short titles for the per-scenario summary in the Conclusion. The method is
# named by its `calcm` value, so a row points straight at the code that produced
# it; the Results section spells the methods out in full.
const SCENARIO_TITLES = Dict(
    "PKG-TEST"              => "Package test suite",
    "PK-LINT-EV-DT0"        => "Plasma :lint, EV, dt 0, no tau, dose 100",
    "PK-LULD-EV-DT025-TAU9" => "Plasma :luld, EV, dt 0.25, tau 9, dose 100",
    "PK-LINT-IV-TAU12"      => "Plasma :lint, IV, dt 0, tau 12, dose 120",
    "PK-LOGT-EV-TAU12"      => "Plasma :logt, EV, dt 0, tau 12, dose 120",
    "URINE-LINT"            => "Urine :lint, dose 100",
    "PD-LINT"               => "PD :lint, BL 3.0, TH 1.5",
    "CONCLUSION"            => "Conclusion rendered",
)

# Every comparison, pass or fail.
const ACCEPT_LOG = joinpath(@__DIR__, "src", "validation_acceptance.tsv")
# Divergences only - the journal. Always written, even when empty.
const ACCEPT_JOURNAL = joinpath(@__DIR__, "src", "validation_divergences.tsv")

const ACCEPT_LOG_HEADER = ["scenario", "parameter", "subject", "value", "reference", "difference", "status"]
const ACCEPT_JOURNAL_HEADER = ["scenario", "parameter", "subject", "value", "reference", "difference", "relative_error", "criterion"]

# How many divergences are listed in the PDF before the table is truncated and
# the reader is sent to the journal file. The journal file is never truncated.
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

# A numeric comparison carries a value; a recorded outcome such as the test
# suite run does not. The two are counted separately so that the reason for a
# failure is never ambiguous.
_is_check(r) = r.value != "-"

# Bookkeeping marks the report makes about itself. The Conclusion records its
# mark as the last thing it does, so while the Conclusion is being rendered its
# own mark does not exist yet. The report-side view (`inreport = true`) leaves
# these out, otherwise every Conclusion would report itself as missing; the
# build gate keeps them, which is the whole point of the mark.
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

_num(s) = s == "-" ? nothing : tryparse(Float64, s)

# Relative error |v - r| / |r|, as a string for display. Returns "n/a" when the
# reference is zero (the criterion then reduces to the absolute tolerance) or
# when either side is not a number.
function _relerr_str(value, reference)
    v, r = _num(value), _num(reference)
    (isnothing(v) || isnothing(r)) && return "n/a"
    (isnan(v) || isnan(r)) && return "NaN"
    iszero(r) && return "n/a"
    return string(round(abs(v - r) / abs(r), sigdigits = 3))
end

# The criterion the comparison was judged against, spelled out per row, so the
# journal is self-contained.
function _criterion_str(reference)
    r = _num(reference)
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
        vals = dsnca[:, l]
        refs = refdict[l]
        for i in eachindex(subjects)
            v = float(vals[i])
            r = float(refs[i])
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
        v = float(isnothing(index) ? result[k] : result[index, k])
        r = float(refdict[k])
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
    v = _num(s)
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
    show(io, latexify(dft; latex = false))
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
        println(msg, "VALIDATION FAILED: acceptance log at $ACCEPT_LOG is empty - no comparison was performed.")
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
