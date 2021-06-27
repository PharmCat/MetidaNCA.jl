#using MetidaNCA
using Test
using DataFrames, CSV

path     = dirname(@__FILE__)
io       = IOBuffer();
pkdata2  = CSV.File(joinpath(path, "csv", "pkdata2.csv")) |> DataFrame

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


@testset "  Linear trapezoidal, Dose 100, Dosetime 0, no tau         " begin

ds = MetidaNCA.pkimport(pkdata2, :Time, :Concentration, [:Subject, :Formulation]; dosetime = MetidaNCA.DoseTime(dose = 100, time = 0))
sort!(ds, :Subject)

sbj = MetidaNCA.nca!(ds[1])

# Linear trapezoidal, Dose 100, Dosetime 0, no tau

dsnca = MetidaNCA.nca!(ds, adm = :ev, calcm = :lint)

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
    @test round.(dsnca[:, :AUClast], sigdigits = 6) == round.([9585.4218
    10112.176
    5396.5498
    9317.8358
    9561.26
    6966.598
    7029.5735
    7110.6745
    8315.0803
    5620.8945], sigdigits = 6)

    # AUMClast
    @test round.(dsnca[:, :AUMClast], sigdigits = 6) == round.([333582.48
    298701.39
    186032.06
    313955.9
    315181.56
    226977.06
    219797.71
    240526.05
    277613.98
    154893.06], sigdigits = 6)

    # AUCall
    @test round.(dsnca[:, :AUCall], sigdigits = 6) == round.([9585.4218
    10112.1760
    5396.5498
    9317.8358
    9561.2600
    6966.5980
    7029.5735
    7110.6745
    8315.0803
    5620.8945], sigdigits = 6)

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
    @test round.(dsnca[:, :AUCinf], sigdigits = 6) == round.([42925.019
    16154.93
    26026.183
    22004.078
    25820.275
    16001.76
    11688.953
    15446.21
    24865.246
    8075.3242], sigdigits = 6)

    # AUCinf_pred

    # AUMCinf

    # AUMCinf_pred

    # AUCpct
    @test round.(dsnca[:, :AUCpct], sigdigits = 5) == round.([77.669383
    37.405019
    79.26492
    57.65405
    62.969953
    56.463551
    39.861391
    53.964925
    66.559429
    30.394194], sigdigits = 5)

    # MRTlast
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

    # MRTinf
    @test round.(dsnca[:, :MRTinf], digits = 5) == round.([293.16224
    71.937917
    305.04073
    130.69968
    149.96684
    128.24114
    79.498252
    114.8571
    176.97811
    58.746446], digits = 5)

    # MRTinf_pred

    # Cllast

    # Clinf
    @test round.(dsnca[:, :Clinf], sigdigits = 6) == round.([0.0023296437
    0.0061900608
    0.0038422846
    0.0045446122
    0.0038729255
    0.0062493127
    0.0085550864
    0.0064740799
    0.0040216775
    0.012383404], sigdigits = 6)

    # Vzlast

    # Vzinf
    @test round.(dsnca[:, :Vzinf], sigdigits = 6) == round.([0.68827768
    0.43881487
    1.1673601
    0.59056646
    0.56843374
    0.8124135
    0.68666158
    0.72497447
    0.71232266
    0.72039519], sigdigits = 6)

    # Vssinf

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


@testset "  setkelauto!                                              " begin
end

@testset "  setdosetime!                                             " begin
end

@testset "  setkelrange!                                             " begin
end
