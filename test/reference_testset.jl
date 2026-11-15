################################################################################
# Check against the combined reference set
# test/csv/reference_values.csv 
# test/csv/reference_scenarios.cs
# (docs/acceptance.jl): ref_rtol 1e-4.
################################################################################

include("references.jl")

const REF_RTOL = 1.0e-4

# Number and description of each scenario. Numbers 1..10 follow the blocks in
# tests.jl, so that a check can be matched to its block:
#
#   #5 and #9 - origin = internal: :luldt and :luld/IV, which Phoenix does not
#               support. The values are taken from test/refdicts.jl (refdict5
#               and refdict9);
#   #10       - in tests.jl this is the same scenario as #3 (120 / 0 / tau 12 /
#               IV / lint)
#
const REF_TESTS = Dict{String, Tuple{Int, String}}(
    "pk_d100_dt0_notau_lint"    => ( 1, "Linear trapezoidal, Dose 100, Dosetime 0, no tau"),
    "pk_d100_dt0.25_tau9_luld"  => ( 2, "Linear up Log down, Dose 100, Dosetime 0.25, tau 9"),
    "pk_d120_dt0_tau12_lint_iv" => ( 3, "Linear trapezoidal, IV, Dose 120, Dosetime 0.0, tau 12"),
    "pk_d120_dt0_tau12_logt"    => ( 4, "Log trapezoidal ATM, Dose 120, Dosetime 0, tau 12"),
    "pk_d120_dt0_tau12_luldt"   => ( 5, "Linear up Log down ATM, Dose 120, Dosetime 0, tau 12"),
    "pk_d120_dt0_tau12_luld"    => ( 6, "Linear up Log down, Dose 120, Dosetime 0, tau 12"),
    "pk_d100_dt2_tau10_lint"    => ( 7, "Linear trapezoidal, Dose 100, Dosetime 2.0, tau 10"),
    "pk_d100_dt0.25_tau9_luld_iv" => ( 9, "Linear up Log down, Dose 100, Dosetime 0.25, tau 9 IV"),
    # in tests.jl block #8 is titled "Linear trapezoidal" but uses :luld; the
    # description here follows the code
    "pk_d100_dt0_tau100_luld"   => ( 8, "Linear up Log down, Dose 100, Dosetime 0.0, tau 100"),
    "pk_d100_dt0_tau36_luld"    => (11, "Linear up Log down, Dose 100, Dosetime 0, tau 36"),
    "pk_d100_dt0_tau60_luld"    => (12, "Linear up Log down, Dose 100, Dosetime 0, tau 60"),
    "pk_d120_dt0_tau12_lint_ev" => (13, "Linear trapezoidal, Dose 120, Dosetime 0, tau 12"),
    "urine_d100_lint"           => (14, "Urine data; Linear-trapezoidal rule, Dose 100"),
    "pd_bl3_th1.5"              => (15, "Pharmacodynamics data; Linear-trapezoidal rule, BL 3, TH 1.5"),
)

"""
    reftitle(scenario, meta, fallback) -> String
"""
function reftitle(sc, meta, number, descr)
    mark = meta.origin == "internal" ? " [internal]" : ""
    return rpad("  #$(number) $(descr)$(mark)", 60)
end

"""
    refresult(scenario, meta, data) -> (result, bysubject::Bool)
"""
function refresult(meta, pkdata, upkdata, pddata)
    if meta.dataset == "pkdata2.csv"
        dt = isnothing(meta.tau) ?
             MetidaNCA.DoseTime(dose = meta.dose, time = meta.dosetime) :
             MetidaNCA.DoseTime(dose = meta.dose, time = meta.dosetime, tau = meta.tau)
        ds = MetidaNCA.pkimport(pkdata, :Time, :Concentration, [:Subject, :Formulation];
                                dosetime = dt)
        sort!(ds, :Subject)
        return MetidaNCA.nca!(ds; adm = meta.adm, calcm = meta.calcm), true
    elseif meta.dataset == "upkdata.csv"
        upk = MetidaNCA.upkimport(upkdata, :st, :et, :conc, :vol, :subj;
                                  dosetime = MetidaNCA.DoseTime(dose = meta.dose))
        return MetidaNCA.nca!(upk), true
    elseif meta.dataset == "pddata.csv"
        pd = MetidaNCA.pdimport(pddata, :time, :obs;
                                bl = meta.bl, th = meta.th, id = Dict(:subj => 1))
        return MetidaNCA.nca!(pd), false
    else
        error("Unknown data set in reference_scenarios.csv: $(meta.dataset)")
    end
end

#=
Lower bound of the number of comparisons with the reference set.

A reference parameter that is absent from the result is skipped (the line
`haskey(result, par) || continue` below): this skips parameters that do not
exist for the route of administration. The skip is silent, so renaming a
parameter reduces the number of comparisons while the test set stays green.
The final number of comparisons is checked against this bound.
=#
const REF_CHECKS_MIN = 3340

# verbose = true: the summary shows which scenarios were checked and how many values
@testset verbose = true "  Reference set (Extrnal Software)                        " begin

    refchecks = Ref(0)
    refskipped = Dict{String, Vector{Symbol}}()

    # check_references() returns true or throws with the list of differences; the
    # messages about skipped values go through @info, so @test_nowarn does not fit
    @test check_references()

    refsc   = reference_scenarios()
    refval  = reference_values()
    refsubj = reference_subjects()

    @test length(refsc) == length(refval)
    @test !isempty(refval)

    # numbers and descriptions come from REF_TESTS; a scenario that is not there
    # is still checked, under the number after the largest known one
    nextnum = maximum(first.(values(REF_TESTS)))
    order = Dict{String, Tuple{Int, String}}()
    for sc in sort(collect(keys(refval)))
        if haskey(REF_TESTS, sc)
            order[sc] = REF_TESTS[sc]
        else
            nextnum += 1
            order[sc] = (nextnum, sc)
            @info "Scenario without a description in REF_TESTS, added in order" scenario = sc number = nextnum
        end
    end

    for sc in sort(collect(keys(refval)); by = x -> first(order[x]))
        number, descr = order[sc]
        meta = refsc[sc]
        @testset "$(reftitle(sc, meta, number, descr))" begin
            res, bysubject = refresult(meta, pkdata2, upkdata, pddata)

            if bysubject
                @test MetidaNCA.getid(res, :, bysubject && meta.dataset == "pkdata2.csv" ?
                                            :Subject : :subj) == refsubj[sc]
            end

            checked = 0
            for (par, ref) in sort(collect(refval[sc]); by = first)
                result = bysubject ? res[1].result : res.result
                if !haskey(result, par)                  # MetidaNCA has no such parameter
                    push!(get!(refskipped, sc, Symbol[]), par)
                    continue
                end
                col = bysubject ? res[:, par] : [result[par]]
                any(ismissing, col) && continue
                vals = float.(col)
                length(vals) == length(ref) || continue
                for i in eachindex(vals)
                    # gap in the export (Swing at Cmin = 0)
                    isnan(ref[i]) && continue
                    checked += 1
                    @test isapprox(vals[i], ref[i]; rtol = REF_RTOL) ||
                          error("$sc / $par / subject $(refsubj[sc][i]): " *
                                "computed $(vals[i]), reference $(ref[i])")
                end
            end
            @test checked > 0                            # the scenario checked something
            refchecks[] += checked
        end
    end

    for sc in sort(collect(keys(refskipped)))
        @info "Reference parameters absent from the result" scenario = sc parameters = sort(unique(refskipped[sc]))
    end
    @info "Comparisons with the reference set" checks = refchecks[] minimum = REF_CHECKS_MIN
    @test refchecks[] >= REF_CHECKS_MIN
end
