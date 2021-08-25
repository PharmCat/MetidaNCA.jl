#using MetidaNCA
using Test
using DataFrames, CSV, Plots

path     = dirname(@__FILE__)
io       = IOBuffer();
pkdata2  = CSV.File(joinpath(path, "csv", "pkdata2.csv")) |> DataFrame
missingpk  = CSV.File(joinpath(path, "csv", "missingpk.csv")) |> DataFrame
aucallpk  = CSV.File(joinpath(path, "csv", "aucalltest.csv")) |> DataFrame
include("refdicts.jl")
# Cmax
# Tmax
# Cdose
# Tlag
# Clast
# AUClast
# AUMClast / AUMCtau
# AUCall
# Rsq
# Adjusted Rsq
# Kel
# HL
# LZint
# Clast_pred
# AUCinf
# AUCinf_pred
# AUMCinf
# AUMCinf_pred
# AUCpct
# MRTlast
# MRTinf / MRTtauinf
# MRTinf_pred
# Cllast
# Clinf / Cltau
# Vzlast
# Vzinf / Vztau
# Vssinf

# AUCtau
# Ctau
# Cavg
# Ctaumin
# Accind
# Fluc
# Fluctau
# Swing
# Swingtau
@testset "   Simple test                                             " begin
    # Basic dataset scenario
    ds = MetidaNCA.pkimport(pkdata2, :Time, :Concentration, [:Subject, :Formulation]; dosetime = MetidaNCA.DoseTime(dose = 100, time = 0))
    sort!(ds, :Subject)
    show(io, ds)
    dsnca = MetidaNCA.nca!(ds, adm = :ev, calcm = :lint)
    @test  MetidaNCA.getid(dsnca, :, :Subject) == collect(1:10)
    show(io, dsnca)

    mtds = MetidaNCA.metida_table(ds)
    @test size(mtds, 1) == size(pkdata2, 1)

    dsncafromds = MetidaNCA.nca(pkdata2, :Time, :Concentration, [:Subject, :Formulation])
    sort!(dsncafromds, :Subject)
    @test dsnca[:, :AUClast] == dsncafromds[:, :AUClast]

    # Plotting
    MetidaNCA.pkplot(ds; typesort = :Subject, pagesort = nothing, sort = Dict(:Formulation => "R"))
    MetidaNCA.pkplot(ds; typesort = :Formulation, pagesort = nothing, legend = true)
    pl = MetidaNCA.pkplot(ds; elim = true, ls = true)
    pl = MetidaNCA.pkplot(ds; typesort = :Subject, pagesort = :Formulation, elim = true, ls = true)
    pl = MetidaNCA.pkplot(ds; typesort = :Formulation, pagesort = :Subject, xticksn = 8, yticksn = 10)

    MetidaNCA.setdosetime!(ds, MetidaNCA.DoseTime(dose = 100, time = 0.25))
    @test first(MetidaNCA.nca!(ds)[:, :Cdose]) == 0

    # Single subject scenario
    tdat = pkdata2[1:16, :Time]
    cdat = pkdata2[1:16, :Concentration]
    ds = MetidaNCA.pkimport(tdat, cdat)
    show(io, ds)
    show(io, MetidaNCA.getdosetime(ds))
    show(io, MetidaNCA.getkelrange(ds))
    sbj = MetidaNCA.nca!(ds)
    show(io, sbj)
    show(io, MetidaNCA.getkeldata(sbj))
    ct = MetidaNCA.ctmax(ds)
    @test  sbj[:Cmax] == ct[1]
    @test  sbj[:Tmax] == ct[2]

    dsncafromds =  MetidaNCA.nca(pkdata2[1:16, :], :Time, :Concentration)
    @test  sbj[:AUClast]  ≈ dsncafromds[:AUClast]

    dsncafromds = MetidaNCA.nca(tdat, cdat)
    @test  sbj[:AUClast]  ≈ dsncafromds[:AUClast]

    # Missing NaN

    dsncafromds =  MetidaNCA.nca(missingpk, :Time, :Concentration)
    @test  sbj[:AUClast]  ≈ dsncafromds[:AUClast]

    dsncafromds =  MetidaNCA.nca(missingpk, :Time, :Concentration;
    limitrule = MetidaNCA.LimitRule(;lloq = 0, btmax = 0, atmax = NaN, nan = NaN, rm = true))
    @test  sbj[:AUClast]  ≈ dsncafromds[:AUClast]

    function newparam(data, result)
        result[:AUChalf] = result[:AUClast] / 2
    end

    dsncafromds =  MetidaNCA.nca(missingpk, :Time, :Concentration;
    limitrule = MetidaNCA.LimitRule(;lloq = 0, btmax = 0, atmax = NaN, nan = NaN, rm = true), modify! = newparam)
    dsncafromds[:AUChalf] ≈ dsncafromds[:AUClast] / 2

    #redirect_stderr(Base.DevNull())
    missingpk.ConcentrationStr = string.(missingpk.Concentration)
    @test_logs (:warn, "Some concentration values not a number, try to fix") pkiw = MetidaNCA.pkimport(missingpk, :Time, :ConcentrationStr)

end

@testset "  #1 Linear trapezoidal, Dose 100, Dosetime 0, no tau      " begin

ds = MetidaNCA.pkimport(pkdata2, :Time, :Concentration, [:Subject, :Formulation]; dosetime = MetidaNCA.DoseTime(dose = 100, time = 0))
sort!(ds, :Subject)
dsnca = MetidaNCA.nca!(ds, adm = :ev, calcm = :lint)

    # Cmax
    @test dsnca[:, :Cmax] == refdict[:Cmax]

    # Tmax
    @test dsnca[:, :Tmax] == refdict[:Tmax]

    # Cdose
    @test round.(dsnca[:, :Cdose], sigdigits = 6) == round.(refdict[:Cdose], sigdigits = 6)

    # Tlag
    #=
    @test round.(dsnca[:, :Tlag], sigdigits = 6) == round.([0
    0
    0
    0
    0.5
    0
    0
    0
    0
    0], sigdigits = 6)
    =#

    # Clast
    @test dsnca[:, :Clast] == refdict[:Clast]

    # AUClast
    @test round.(dsnca[:, :AUClast], sigdigits = 6) == round.(refdict[:AUClast], sigdigits = 6)

    # AUMClast
    @test round.(dsnca[:, :AUMClast], sigdigits = 6) == round.(refdict[:AUMClast], sigdigits = 6)

    # AUCall
    @test round.(dsnca[:, :AUCall], sigdigits = 6) == round.(refdict[:AUCall], sigdigits = 6)

    # Rsq
    @test round.(dsnca[:, :Rsq], digits = 6) == round.(refdict[:Rsq], digits = 6)

    # Adjusted Rsq
    @test round.(dsnca[:, :ARsq], digits = 6) == round.(refdict[:ARsq], digits = 6)

    # Kel
    @test round.(dsnca[:, :Kel], sigdigits = 6) == round.(refdict[:Kel], sigdigits = 6)

    # HL
    @test round.(dsnca[:, :HL], sigdigits = 5) == round.(refdict[:HL], sigdigits = 5)

    # Clast_pred
    @test round.(dsnca[:, :Clast_pred], sigdigits = 6) == round.(refdict[:Clast_pred], sigdigits = 6)

    # AUCinf
    @test round.(dsnca[:, :AUCinf], sigdigits = 6) == round.(refdict[:AUCinf], sigdigits = 6)

    # AUCinf_pred

    # AUMCinf

    # AUMCinf_pred

    # AUCpct
    @test round.(dsnca[:, :AUCpct], sigdigits = 5) == round.(refdict[:AUCpct], sigdigits = 5)

    # MRTlast
    @test round.(dsnca[:, :MRTlast], digits = 6) == round.(refdict[:MRTlast], digits = 6)

    # MRTinf
    @test round.(dsnca[:, :MRTinf], digits = 5) == round.(refdict[:MRTinf], digits = 5)

    # MRTinf_pred

    # Cllast

    # Clinf
    @test round.(dsnca[:, :Clinf], sigdigits = 6) == round.(refdict[:Clinf], sigdigits = 6)

    # Vzlast

    # Vzinf
    @test round.(dsnca[:, :Vzinf], sigdigits = 6) == round.(refdict[:Vzinf], sigdigits = 6)

    # Vssinf

end

@testset "  #2 Linear up Log down, Dose 100, Dosetime 0.25, tau 9    " begin
    ds = MetidaNCA.pkimport(pkdata2, :Time, :Concentration, [:Subject, :Formulation]; dosetime = MetidaNCA.DoseTime(dose = 100, time = 0.25, tau = 9))
    sort!(ds, :Subject)
    dsnca = MetidaNCA.nca!(ds, adm = :ev, calcm = :luld)
    # Cmax
    @test dsnca[:, :Cmax] == refdict2[:Cmax]

    # Tmax
    @test dsnca[:, :Tmax] == refdict2[:Tmax]

    # Cdose
    @test round.(dsnca[:, :Cdose], sigdigits = 6) == round.(refdict2[:Cdose], sigdigits = 6)

    # Tlag
    #=
    @test round.(dsnca[:, :Tlag], sigdigits = 6) == round.([0
    0
    0
    0
    0.25
    0
    0
    0
    0
    0], sigdigits = 6)
    =#

    # Clast
    @test dsnca[:, :Clast] == refdict2[:Clast]

    # AUClast
    @test round.(dsnca[:, :AUClast], digits = 6) == round.(refdict2[:AUClast], digits = 6)
    # AUMClast / AUMCtau

    @test round.(dsnca[:, :AUMCtau], digits = 6) ==  round.(refdict2[:AUMCtau], digits = 6)

# AUCall
    @test round.(dsnca[:, :AUCall], digits = 6) == round.(refdict2[:AUCall], digits = 6)
# Rsq
# Adjusted Rsq
    # LZint
    @test round.(dsnca[:, :LZint], digits = 6) == round.(refdict2[:LZint], digits = 6)
# Kel
    @test round.(dsnca[:, :Kel], digits = 6) == round.(refdict2[:Kel], digits = 6)
# HL
    @test round.(dsnca[:, :HL], digits = 6) == round.(refdict2[:HL], digits = 6)
# Clast_pred
    @test round.(dsnca[:, :Clast_pred], digits = 6) == round.(refdict2[:Clast_pred], digits = 6)
# AUCinf
    @test round.(dsnca[:, :AUCinf], digits = 6) == round.(refdict2[:AUCinf], digits = 6)
# AUCinf_pred
# AUMCinf
# AUMCinf_pred
# AUCpct
    @test round.(dsnca[:, :AUCpct], digits = 6) == round.(refdict2[:AUCpct], digits = 6)
# MRTlast

# MRTinf / MRTtauinf
    @test round.(dsnca[:, :MRTtauinf], digits = 6) == round.(refdict2[:MRTtauinf], digits = 6)
# MRTinf_pred
# Cllast
# Clinf / Cltau
    @test round.(dsnca[:, :Cltau], digits = 6) == round.(refdict2[:Cltau], digits = 6)
# Vzlast
# Vzinf / Vztau
    @test round.(dsnca[:, :Vztau], digits = 6) == round.(refdict2[:Vztau], digits = 6)
# Vssinf

    # AUCtau
    @test round.(dsnca[:, :AUCtau], digits = 6) == round.(refdict2[:AUCtau], digits = 6)

end


@testset "  #3 Linear trapezoidal, IV, Dose 120, Dosetime 0.0, tau 12" begin
    ds = MetidaNCA.pkimport(pkdata2, :Time, :Concentration, [:Subject, :Formulation]; dosetime = MetidaNCA.DoseTime(dose = 120, time = 0, tau = 12))
    sort!(ds, :Subject)
    dsnca = MetidaNCA.nca!(ds, adm = :iv, calcm = :lint)

    # Cmax
    @test dsnca[:, :Cmax] == refdict3[:Cmax]

    # Tmax
    @test dsnca[:, :Tmax] == refdict3[:Tmax]

    # Cdose
    @test round.(dsnca[:, :Cdose], sigdigits = 6) == round.(refdict3[:Cdose], sigdigits = 6)

    # Tlag
    #@test round.(dsnca[:, :Tlag], sigdigits = 6) == round.(refdict3[:Tlag], sigdigits = 6)

    # Clast
    @test dsnca[:, :Clast] == refdict3[:Clast]

    # AUClast
    @test round.(dsnca[:, :AUClast], sigdigits = 6) == round.(refdict3[:AUClast], sigdigits = 6)
    # AUMClast / AUMCtau

    @test round.(dsnca[:, :AUMCtau], sigdigits = 6) ==  round.(refdict3[:AUMCtau], sigdigits = 6)

# AUCall
    @test round.(dsnca[:, :AUCall], sigdigits = 6) == round.(refdict3[:AUCall], sigdigits = 6)
# Rsq
# Adjusted Rsq
    # LZint
    @test round.(dsnca[:, :LZint], sigdigits = 6) == round.(refdict3[:LZint], sigdigits = 6)
# Kel
    @test round.(dsnca[:, :Kel], sigdigits = 6) == round.(refdict3[:Kel], sigdigits = 6)
# HL
    @test round.(dsnca[:, :HL], sigdigits = 6) == round.(refdict3[:HL], sigdigits = 6)
# Clast_pred
    @test round.(dsnca[:, :Clast_pred], sigdigits = 6) == round.(refdict3[:Clast_pred], sigdigits = 6)
# AUCinf
    @test round.(dsnca[:, :AUCinf], sigdigits = 6) == round.(refdict3[:AUCinf], sigdigits = 6)
# AUCinf_pred
# AUMCinf
# AUMCinf_pred
# AUCpct
    @test round.(dsnca[:, :AUCpct], sigdigits = 6) == round.(refdict3[:AUCpct], sigdigits = 6)
# MRTlast

# MRTinf / MRTtauinf
    @test round.(dsnca[:, :MRTtauinf], sigdigits = 6) == round.(refdict3[:MRTtauinf], sigdigits = 6)
# MRTinf_pred
# Cllast
# Clinf / Cltau
    @test round.(dsnca[:, :Cltau], sigdigits = 6) == round.(refdict3[:Cltau], sigdigits = 6)
# Vzlast
# Vzinf / Vztau
    @test round.(dsnca[:, :Vztau], sigdigits = 6) == round.(refdict3[:Vztau], sigdigits = 6)
# Vssinf

    # AUCtau
    @test round.(dsnca[:, :AUCtau], sigdigits = 6) == round.(refdict3[:AUCtau], sigdigits = 6)

    # AUClast
    # AUMClast / AUMCtau
    # AUCall
    # Rsq
    # Adjusted Rsq
    # Kel
    # HL
    # LZint
    # Clast_pred
    # AUCinf
    # AUCinf_pred
    # AUMCinf
    # AUMCinf_pred
    # AUCpct
    # MRTlast
    # MRTinf / MRTtauinf
    # MRTinf_pred
    # Cllast
    # Clinf / Cltau
    # Vzlast
    # Vzinf / Vztau
    # Vssinf

    # AUCtau
    # Ctau
    # Cavg
    # Ctaumin
    # Accind
    # Fluc
    # Fluctau
    # Swing
    # Swingtau

end



@testset "  Linear up Log down, Dose 120, Dosetime 0, tau 12         " begin
    ds = MetidaNCA.pkimport(pkdata2, :Time, :Concentration, [:Subject, :Formulation]; dosetime = MetidaNCA.DoseTime(dose = 120, time = 0, tau = 12))
    sort!(ds, :Subject)
    dsnca = MetidaNCA.nca!(ds, adm = :ev, calcm = :luld)
    # Cmax
    @test dsnca[:, :Cmax] == [190.869
    261.177
    105.345
    208.542
    169.334
    154.648
    153.254
    138.327
    167.347
    125.482]

    # Tmax
    @test dsnca[:, :Tmax] == [1
    1
    1.5
    1
    4
    2.5
    2.5
    4
    3
    2]

    # Cdose
    @test round.(dsnca[:, :Cdose], sigdigits = 6) == round.([0.0
    0.0
    0.0
    0.0
    0.0
    0.0
    0.0
    0.0
    0.0
    0.0], sigdigits = 6)

    # Tlag
    #=
    @test round.(dsnca[:, :Tlag], sigdigits = 6) == round.([0
    0
    0
    0
    0.5
    0
    0
    0
    0
    0], sigdigits = 6)
    =#

    # Clast
    @test dsnca[:, :Clast] == [112.846
    85.241
    67.901
    97.625
    110.778
    69.501
    58.051
    74.437
    93.44
    42.191]

    # AUClast
    @test round.(dsnca[:, :AUClast], sigdigits = 6) == round.([9573.8106
    10054.286
    5392.4572
    9297.0963
    9519.1809
    6948.9856
    6988.7726
    7073.0922
    8303.3586
    5486.8389], sigdigits = 6)

    # AUMCtau
    @test round.(dsnca[:, :AUMCtau], sigdigits = 6) == round.([9995.4521
14645.001
6030.6644
10314.429
11475.22
8486.2318
9023.7903
6403.6719
10120.542
8376.0054], sigdigits = 6)

    # AUCall
    @test round.(dsnca[:, :AUCall], sigdigits = 6) == round.([9573.8106
10054.286
5392.4572
9297.0963
9519.1809
6948.9856
6988.7726
7073.0922
8303.3586
5486.8389
], sigdigits = 6)

    # Rsq
    @test round.(dsnca[:, :Rsq], digits = 6) == round.([0.78607696
0.99276359
0.81358898
0.91885869
0.85335995
0.95011904
0.97031231
0.94796904
0.94753789
0.88092269], digits = 6)

    # Adjusted Rsq
    @test round.(dsnca[:, :ARsq], digits = 6) == round.([0.71476928
    0.99035145
    0.77630678
    0.83771737
    0.82891994
    0.92517856
    0.96041642
    0.92195356
    0.92130684
    0.86391165], digits = 6)

    # Kel
    @test round.(dsnca[:, :Kel], sigdigits = 6) == round.([0.0033847439
    0.014106315
    0.0032914304
    0.0076953442
    0.0068133279
    0.0076922807
    0.012458956
    0.0089300798
    0.0056458649
    0.017189737], sigdigits = 6)

    # HL
    @test round.(dsnca[:, :HL], sigdigits = 5) == round.([204.78571
    49.137367
    210.59148
    90.073577
    101.73401
    90.10945
    55.634451
    77.619371
    122.77077
    40.323315], sigdigits = 5)

    # Clast_pred
    @test round.(dsnca[:, :Clast_pred], sigdigits = 6) == round.([117.30578
    82.53669
    66.931057
    100.76793
    105.29832
    71.939942
    61.172702
    75.604277
    93.761762
    38.810857], sigdigits = 6)

    # AUCinf
    @test round.(dsnca[:, :AUCinf], sigdigits = 6) == round.([42913.408
16097.041
26022.09
21983.338
25778.196
15984.147
11648.152
15408.628
24853.524
7941.2686], sigdigits = 6)

    # AUCinf_pred

    # AUMCinf

    # AUMCinf_pred

    # AUCpct
    @test round.(dsnca[:, :AUCpct], sigdigits = 5) == round.([77.690398
37.539537
79.277386
57.708442
63.072742
56.525766
40.001017
54.096548
66.59082
30.907274], sigdigits = 5)

    # MRTlast
    #=
    @test round.(dsnca[:, :MRTlast], digits = 6) == round.([34.801023
    29.538786
    34.472406
    33.69408
    32.964438
    32.58076
    31.267574
    33.826053
    33.386807
    27.556657], digits = 6)
=#
    # MRTtauinf
    @test round.(dsnca[:, :MRTtauinf], digits = 4) == round.([302.522490
75.321653
312.895820
148.293830
172.561490
130.198390
91.786365
163.779880
176.558310
68.172066], digits = 4 )

    # MRTinf_pred

    # Cllast

    # Cltau
    @test round.(dsnca[:, :Cltau], sigdigits = 6) == round.([0.071896833
0.050424060
0.122488280
0.070172356
0.069035042
0.085192949
0.083697775
0.110260280
0.073377834
0.092799594], sigdigits = 6)

    # Vzlast

    # Vztau
    @test round.(dsnca[:, :Vztau], sigdigits = 6) == round.([21.241439
3.574574
37.214301
9.118807
10.132353
11.075122
6.717880
12.347065
12.996739
5.398547], sigdigits = 6)

    # Vssinf
end


@testset "  Log trapezoidal ATM, Dose 120, Dosetime 0, tau 12        " begin
    ds = MetidaNCA.pkimport(pkdata2, :Time, :Concentration, [:Subject, :Formulation]; dosetime = MetidaNCA.DoseTime(dose = 120, time = 0, tau = 12))
    sort!(ds, :Subject)
    dsnca = MetidaNCA.nca!(ds, adm = :ev, calcm = :logt)
    @test dsnca[:, :Cmax] == [190.869
    261.177
    105.345
    208.542
    169.334
    154.648
    153.254
    138.327
    167.347
    125.482]

    # Tmax
    @test dsnca[:, :Tmax] == [1
    1
    1.5
    1
    4
    2.5
    2.5
    4
    3
    2]

    # Cdose
    @test round.(dsnca[:, :Cdose], sigdigits = 6) == round.([0.0
    0.0
    0.0
    0.0
    0.0
    0.0
    0.0
    0.0
    0.0
    0.0], sigdigits = 6)

    # Tlag
    #=
    @test round.(dsnca[:, :Tlag], sigdigits = 6) == round.([0
    0
    0
    0
    0.5
    0
    0
    0
    0
    0], sigdigits = 6)
    =#

    # Clast
    @test dsnca[:, :Clast] == [112.846
    85.241
    67.901
    97.625
    110.778
    69.501
    58.051
    74.437
    93.44
    42.191]

    # AUClast
    @test round.(dsnca[:, :AUClast], sigdigits = 6) == round.([9572.8582
    10054.0370
    5391.5322
    9296.2179
    9518.6531
    6948.5757
    6987.0645
    7064.7816
    8298.9634
    5485.6538], sigdigits = 6)

end


@testset "  Linear up Log down ATM, Dose 120, Dosetime 0, tau 12     " begin
    ds = MetidaNCA.pkimport(pkdata2, :Time, :Concentration, [:Subject, :Formulation]; dosetime = MetidaNCA.DoseTime(dose = 120, time = 0, tau = 12))
    sort!(ds, :Subject)
    dsnca = MetidaNCA.nca!(ds, adm = :ev, calcm = :luldt)
    # Cmax
    @test dsnca[:, :Cmax] == [190.869
    261.177
    105.345
    208.542
    169.334
    154.648
    153.254
    138.327
    167.347
    125.482]

    # Tmax
    @test dsnca[:, :Tmax] == [1
    1
    1.5
    1
    4
    2.5
    2.5
    4
    3
    2]

    # Cdose
    @test round.(dsnca[:, :Cdose], sigdigits = 6) == round.([0.0
    0.0
    0.0
    0.0
    0.0
    0.0
    0.0
    0.0
    0.0
    0.0], sigdigits = 6)

    # Tlag
    #=
    @test round.(dsnca[:, :Tlag], sigdigits = 6) == round.([0
    0
    0
    0
    0.5
    0
    0
    0
    0
    0], sigdigits = 6)
    =#

    # Clast
    @test dsnca[:, :Clast] == [112.846
    85.241
    67.901
    97.625
    110.778
    69.501
    58.051
    74.437
    93.44
    42.191]

end



@testset "  Linear trapezoidal, Dose 100, Dosetime 2.0, tau 10       " begin
    ds = MetidaNCA.pkimport(pkdata2, :Time, :Concentration, [:Subject, :Formulation]; dosetime = MetidaNCA.DoseTime(dose = 100, time = 2.0, tau = 10))
    sort!(ds, :Subject)
    dsnca = MetidaNCA.nca!(ds, adm = :ev, calcm = :lint)

    # Cmax
    @test dsnca[:, :Cmax] == [150.854
    234.091
    100.943
    165.177
    169.334
    154.648
    153.254
    138.327
    167.347
    125.482]

    # Tmax
    @test dsnca[:, :Tmax] == [4
    2
    2
    2
    4
    2.5
    2.5
    4
    3
    2]

    # Cdose
    # Tlag
    # Clast
    # AUClast

    @test round.(dsnca[:, :AUClast], digits = 3) == round.([9283.0588
    9774.9218
    5255.0335
    9051.5985
    9441.4205
    6781.2608
    6858.7153
    6885.9150
    8094.8320
    5495.0990], digits = 3)

    # AUMClast / AUMCtau
    @test round.(dsnca[:, :AUMCtau], digits = 3) == round.([6915.4913
    10105.3000
    4166.2103
    7068.4668
    8032.4660
    5775.6840
    6243.2670
    4442.4953
    6999.9440
    5849.5573], digits = 3)
    # AUCall
    # Rsq
    # Adjusted Rsq
    # Kel
    # HL
    # LZint
    @test round.(dsnca[:, :LZint], sigdigits = 6) == round.([5.0084856
    5.4288976
    4.4406461
    5.166885
    5.1473571
    4.8296758
    5.0107459
    4.9684786
    4.9472594
    4.8963611], sigdigits = 6)
    # Clast_pred
    @test round.(dsnca[:, :Clast_pred], digits = 3) == round.([117.30578
    82.53669
    66.931057
    100.76793
    105.29832
    71.939942
    61.172702
    75.604277
    93.761762
    38.810857], digits = 3)
    # AUCinf
    @test round.(dsnca[:, :AUCinf], digits = 3) == round.([42622.6560
    15817.6760
    25884.6660
    21737.8410
    25700.4350
    15816.4220
    11518.0940
    15221.4510
    24644.9980
    7949.5287], digits = 3)
    # AUCinf_pred
    # AUMCinf
    # AUMCinf_pred
    # AUCpct
    # MRTlast
    # MRTinf / MRTtauinf
    @test round.(dsnca[:, :MRTtauinf], digits = 3) == round.([306.68478
    72.36944
    303.54748
    145.34822
    153.74241
    123.86557
    85.934338
    167.95643
    168.74491
    63.074784], digits = 3)

    # MRTinf_pred
    # Cllast
    # Clinf / Cltau
    # Vzlast
    # Vzinf / Vztau
    # Vssinf

    # AUCtau
    @test round.(dsnca[:, :AUCtau], digits = 3) == round.([1367.7388
    2043.0158
    838.8295
    1444.7985
    1618.6205
    1224.6608
    1265.7013
    880.311
    1417.942
    1167.911], digits = 3)
    # Ctau
    @test round.(dsnca[:, :Ctau], sigdigits = 6) == round.([144.964
    196.035
    76.027
    132.257
    154.066
    113.751
    123.37
    134.133
    135.58
    106.476], sigdigits = 6)
    # Cavg
    @test round.(dsnca[:, :Cavg], digits = 3) == round.([136.77388
    204.30158
    83.88295
    144.47985
    161.86205
    122.46608
    126.57013
    88.0311
    141.7942
    116.7911], digits = 3)
    # Ctaumin
    @test round.(dsnca[:, :Ctaumin], sigdigits = 6) == round.([121.239
    196.026
    69.985
    128.398
    151.452
    102.16
    105.513
    22.724
    105.438
    106.476], sigdigits = 6)
    # Accind
    @test round.(dsnca[:, :Accind], digits = 3) == round.([30.047153
    7.600775
    30.884671
    13.501282
    15.182793
    13.506455
    8.5367345
    11.705549
    18.216783
    6.3317425], digits = 3)
    # Fluc
    @test round.(dsnca[:, :Fluc], sigdigits = 6) == round.([21.652527
    18.631770
    36.906189
    25.456145
    11.047679
    42.859216
    37.719011
    131.320640
    43.661165
    16.273500], sigdigits = 6)

    # Swing
    @test round.(dsnca[:, :Swing], sigdigits = 6) == round.([0.24426958
    0.19418342
    0.44235193
    0.28644527
    0.11807041
    0.51378230
    0.45246557
    5.08726460
    0.58716023
    0.17850032], sigdigits = 6)

end

@testset "  Linear trapezoidal, Dose 100, Dosetime 0.0, tau 100      " begin
    ds = MetidaNCA.pkimport(pkdata2, :Time, :Concentration, [:Subject, :Formulation]; dosetime = MetidaNCA.DoseTime(dose = 100, time = 0.0, tau = 100))
    sort!(ds, :Subject)
    dsnca = MetidaNCA.nca!(ds, adm = :ev, calcm = :luld)

    # Cmax
    @test dsnca[:, :Cmax] == [190.869
    261.177
    105.345
    208.542
    169.334
    154.648
    153.254
    138.327
    167.347
    125.482]

    # Tmax
    @test dsnca[:, :Tmax] == [1
    1
    1.5
    1
    4
    2.5
    2.5
    4
    3
    2]

    # Cdose
    # Tlag
    # Clast
    # AUClast

    @test round.(dsnca[:, :AUClast], digits = 3) == round.([9573.810559
    10054.28648
    5392.457219
    9297.096334
    9519.180874
    6948.985621
    6988.772632
    7073.092214
    8303.358586
    5486.838889], digits = 3)

    # AUMClast / AUMCtau
    @test round.(dsnca[:, :AUMCtau], digits = 3) == round.([599491.4507
    469656.5157
    341841.895
    530699.316
    554228.2569
    381802.405
    343144.4367
    402118.7541
    487552.4608
    233710.5713], digits = 3)
    # AUCall
    # Rsq
    # Adjusted Rsq
    # Kel
    # HL
    # LZint
    @test round.(dsnca[:, :LZint], sigdigits = 6) == round.([5.0084856
    5.4288976
    4.4406461
    5.166885
    5.1473571
    4.8296758
    5.0107459
    4.9684786
    4.9472594
    4.8963611], sigdigits = 6)

    # Clast_pred
    @test round.(dsnca[:, :Clast_pred], digits = 3) == round.([117.30578
    82.53669
    66.931057
    100.76793
    105.29832
    71.939942
    61.172702
    75.604277
    93.761762
    38.810857], digits = 3)
    # AUCinf
    @test round.(dsnca[:, :AUCinf], digits = 3) == round.([42913.4078813004
    16097.0411126277
    26022.0900281352
    21983.3384532182
    25778.1957695968
    15984.1473646863
    11648.1518057779
    15408.6280190766
    24853.5242997128
    7941.2685538530], digits = 3)
    # AUCinf_pred
    # AUMCinf
    # AUMCinf_pred
    # AUCpct
    # MRTlast
    # MRTinf / MRTtauinf
    @test round.(dsnca[:, :MRTtauinf], digits = 3) == round.([286.7298872
    73.32812486
    309.1287691
    131.3893596
    155.1600979
    126.8510359
    79.61736185
    117.5547609
    177.1315649
    60.86551733], digits = 3)

    # MRTinf_pred
    # Cllast
    # Clinf / Cltau
    # Vzlast
    # Vzinf / Vztau
    # Vssinf

    # AUCtau
    @test round.(dsnca[:, :AUCtau], digits = 3) == round.([12646.63632
    11996.6718
    7195.902904
    11794.11692
    12274.83395
    8729.151856
    8395.400098
    8930.999936
    10727.4135
    6389.420453], digits = 3)
    # Ctau
    @test round.(dsnca[:, :Ctau], sigdigits = 6) == round.([106.6989373
    55.60460942
    61.0383917
    81.23535402
    87.01010737
    58.00027625
    43.15724396
    58.8781831
    80.05171762
    23.98401112], sigdigits = 6)
    # Cavg
    @test round.(dsnca[:, :Cavg], digits = 3) == round.([126.4663632
    119.966718
    71.95902904
    117.9411692
    122.7483395
    87.29151856
    83.95400098
    89.30999936
    107.274135
    63.89420453], digits = 3)
    # Ctaumin
    @test round.(dsnca[:, :Ctaumin], sigdigits = 6) == round.([0
    0
    0
    0
    0
    0
    0
    0
    0
    0], sigdigits = 6)
    # Accind
    @test round.(dsnca[:, :Accind], digits = 3) == round.([3.482585727
1.322732351
3.565571971
1.862990767
2.024054788
1.863483514
1.403869629
1.693257481
2.318008606
1.218397838], digits = 3)
    # Fluc
    @test round.(dsnca[:, :Fluc], sigdigits = 6) == round.([150.9247164
217.7078813
146.3958052
176.8186643
137.9521716
177.1626872
182.5452012
154.8841126
155.9993935
196.3902687], sigdigits = 6)

    # Swing

    # Swing tau
    @test round.(dsnca[:, :Swingtau], sigdigits = 6) == round.([0.788855679643705
3.697038657747540
0.725880991766796
1.567133516128610
0.946141719805552
1.666332127921700
2.551060863661350
1.349376164748240
1.090486062900490
4.231902178181850], sigdigits = 6)

end

@testset "  Linear up Log down, Dose 100, Dosetime 0.25, tau 9 IV    " begin
    ds = MetidaNCA.pkimport(pkdata2, :Time, :Concentration, [:Subject, :Formulation]; dosetime = MetidaNCA.DoseTime(dose = 100, time = 0.25, tau = 9))
    sort!(ds, :Subject)
    dsnca = MetidaNCA.nca!(ds, adm = :iv, calcm = :luld)

    @test dsnca[:, :Cdose] == [178.949
    62.222
    49.849
    52.421
    0.0
    57.882
    19.95
    142.34985100539438
    113.362
    13.634]

end

@testset "  Linear trapezoidal, Dose 100, Dosetime 0, no tau AUCall  " begin

    dsnca = MetidaNCA.nca(aucallpk, :Time, :Concentration; dosetime = MetidaNCA.DoseTime(dose = 100, time = 0), adm = :ev, calcm = :lint)
    @test dsnca[:AUClast] ≈ 9585.4218
    aucl =  MetidaNCA.linauc(72, 96, 112.846, 0)
    dsnca[:AUClast] + aucl ≈ dsnca[:AUCall]

end

@testset "  set-get*! tests                                          " begin
    ds = MetidaNCA.pkimport(pkdata2, :Time, :Concentration, [:Subject, :Formulation])
    sort!(ds, :Subject)
    #=
    @testset "  #1 setkelauto!                                            " begin
        ka = MetidaNCA.setkelauto!(ds[1], false)
        @test MetidaNCA.getkelauto(ka) == true

        dsnca = MetidaNCA.nca!(ds, adm = :ev, calcm = :luldt)
    end
    =#
    @testset "  setdosetime!                                             " begin
        dt = MetidaNCA.DoseTime(dose = 110, time = 2.1, tau = 10)
        MetidaNCA.setdosetime!(ds[1], dt)
        dts = MetidaNCA.getdosetime(ds[1])
        @test dts.dose == 110
        @test dts.time == 2.1
        @test dts.tau == 10
        dt2 = MetidaNCA.DoseTime(dose = 100, time = 2.2, tau = 9)
        MetidaNCA.setdosetime!(ds, dt2, 4)
        MetidaNCA.setdosetime!(ds, dt2, [1,2,3])
        MetidaNCA.setdosetime!(ds, dt2, Dict(:Formulation => "R"))
        MetidaNCA.setdosetime!(ds, dt2)

        dsnca = MetidaNCA.nca!(ds, adm = :ev, calcm = :luldt)
    end
    @testset "  #2 setkelauto!                                            " begin
        kr =  MetidaNCA.ElimRange(kelstart = 4, kelend = 12)
        MetidaNCA.setkelrange!(ds, kr; kelauto = true)
        MetidaNCA.setkelauto!(ds, false, 4)
        MetidaNCA.setkelauto!(ds, false, [1,2,3])
        MetidaNCA.setkelauto!(ds, false, Dict(:Formulation => "R"))
        MetidaNCA.setkelauto!(ds, false)
        dsnca = MetidaNCA.nca!(ds, adm = :ev, calcm = :luldt)
        @test MetidaNCA.getkelauto(ds[1]) == false
    end
    ds = MetidaNCA.pkimport(pkdata2, :Time, :Concentration, [:Subject, :Formulation])
    sort!(ds, :Subject)
    @testset "  setkelrange!                                             " begin
        kr =  MetidaNCA.ElimRange(kelstart = 4, kelend = 12, kelexcl = Int[5,6])
        MetidaNCA.setkelrange!(ds[1], kr)
        krs = MetidaNCA.getkelrange(ds[1])
        @test krs.kelstart == 4
        @test krs.kelend == 12
        @test krs.kelexcl == [5,6]
        kr2 =  MetidaNCA.ElimRange(kelstart = 3, kelend = 12, kelexcl = Int[7])
        MetidaNCA.setkelrange!(ds, kr2, 4)
        MetidaNCA.setkelrange!(ds, kr2, [1,2,3])
        MetidaNCA.setkelrange!(ds, kr2, Dict(:Formulation => "R"))
        MetidaNCA.setkelrange!(ds, kr2)

        dsnca = MetidaNCA.nca!(ds, adm = :ev, calcm = :luldt)
    end
end
@testset "  applylimitrule!                                          " begin
    ds = MetidaNCA.pkimport(pkdata2, :Time, :Concentration, [:Subject, :Formulation])
    sort!(ds, :Subject)
    lr = MetidaNCA.LimitRule(;lloq = 0.5, btmax = 0.0, atmax = NaN, nan = NaN, rm = true)
    MetidaNCA.applylimitrule!(ds[1], lr)
    dsnca = MetidaNCA.nca!(ds, adm = :ev, calcm = :luldt)

    lr = MetidaNCA.LimitRule(;lloq = 0.5, btmax = 0.5, atmax = NaN, nan = NaN, rm = true)
    function af(sbj)
        sbj.id[:Subject] == 1
    end
    MetidaNCA.applylimitrule!(af, ds, lr)
    @test ds[1].obs[1] ≈ 0.5
    MetidaNCA.applylimitrule!(ds, lr, 2)
    @test ds[2].obs[1] ≈ 0.5
    MetidaNCA.applylimitrule!(ds, lr, 3:4)
    @test ds[3].obs[1] ≈ 0.5
    @test ds[4].obs[1] ≈ 0.5
    MetidaNCA.applylimitrule!(ds, lr, Dict(:Formulation => "R"))
    @test ds[7].obs[1] ≈ 0.5
    MetidaNCA.applylimitrule!(ds, lr)
    @test ds[6].obs[1] ≈ 0.5

    ds = MetidaNCA.pkimport(missingpk, :Time, :Concentration)

    @test ismissing(ds.obs[13])
    @test isnan(ds.obs[15])
    @test length(ds) == 18
    MetidaNCA.applylimitrule!(ds, lr)
    @test length(ds) == 16

    ds = MetidaNCA.pkimport(missingpk, :Time, :Concentration)
    lr = MetidaNCA.LimitRule(;lloq = 180, btmax = 0.0, atmax = 0.5, nan = 1000, rm = false)
    MetidaNCA.applylimitrule!(ds, lr)
    @test ds.obs ==  [0.0
    0.0
  190.869
    0.5
    0.5
    0.5
    0.5
    0.5
    0.5
    0.5
    0.5
    0.5
 1000.0
    0.5
 1000.0
    0.5
    0.5
    0.5]
end

@testset "  kel                                                      " begin
    ds = MetidaNCA.pkimport(pkdata2, :Time, :Concentration, [:Subject, :Formulation])
    sort!(ds, :Subject)
    kr1 = MetidaNCA.ElimRange(kelstart = 4, kelend = 12, kelexcl = Int[5,6])
    kr2 = MetidaNCA.ElimRange(kelstart = 3, kelend = 12, kelexcl = Int[7])

    MetidaNCA.setkelrange!(ds, kr1, Dict(:Formulation => "T"))
    MetidaNCA.setkelrange!(ds, kr2, Dict(:Formulation => "R"))

    sub1 = MetidaNCA.subset(ds, Dict(:Formulation => "T"))
    @test MetidaNCA.getkelrange(sub1[1]) == kr1
    sub2 = MetidaNCA.subset(ds, Dict(:Formulation => "R"))
    @test MetidaNCA.getkelrange(sub2[1]) == kr2
    dsnca = MetidaNCA.nca!(ds, adm = :ev, calcm = :luldt)
end

@testset "  Output                                                   " begin
    io = IOBuffer();
    ds = MetidaNCA.pkimport(pkdata2, :Time, :Concentration, [:Subject, :Formulation]; dosetime = MetidaNCA.DoseTime(dose = 100, time = 0))
    sort!(ds, :Subject)
    @test_nowarn dsnca = MetidaNCA.nca!(ds, adm = :ev, calcm = :lint, verbose = 1, io = io)
    show(io, ds[1])
    dt = MetidaNCA.DoseTime(dose = 100, time = 0.25, tau = 9)
    MetidaNCA.setdosetime!(ds, dt)
    @test_nowarn dsnca = MetidaNCA.nca!(ds, adm = :ev, calcm = :luld, verbose = 1, io = io)
    show(io, ds[1])
    kr =  MetidaNCA.ElimRange(kelstart = 10, kelend = 16, kelexcl = Int[13,14])
    MetidaNCA.setkelrange!(ds, kr; kelauto = false)
    show(io, ds[1])
    @test_nowarn dsnca = MetidaNCA.nca!(ds, adm = :ev, calcm = :lint, verbose = 2, io = io)
end
