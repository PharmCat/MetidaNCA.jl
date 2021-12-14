
#   Linear trapezoidal, Dose 100, Dosetime 0, no tau

refdict = Dict(
:Cmax => [
    190.869
    261.177
    105.345
    208.542
    169.334
    154.648
    153.254
    138.327
    167.347
    125.482
],
:Tmax => [
    1
    1
    1.5
    1
    4
    2.5
    2.5
    4
    3
    2
],
:Cdose => [
    0.0
    0.0
    0.0
    0.0
    0.0
    0.0
    0.0
    0.0
    0.0
    0.0
],
:Clast => [
    112.846
    85.241
    67.901
    97.625
    110.778
    69.501
    58.051
    74.437
    93.44
    42.191
],
:AUClast => [
    9585.4218
    10112.176
    5396.5498
    9317.8358
    9561.26
    6966.598
    7029.5735
    7110.6745
    8315.0803
    5620.8945
],
:AUMClast => [
    333582.48
    298701.39
    186032.06
    313955.9
    315181.56
    226977.06
    219797.71
    240526.05
    277613.98
    154893.06
],
:AUCall => [
    9585.4218
    10112.1760
    5396.5498
    9317.8358
    9561.2600
    6966.5980
    7029.5735
    7110.6745
    8315.0803
    5620.8945
],
:Rsq => [
    0.78607696
    0.99276359
    0.81358898
    0.91885869
    0.85335995
    0.95011904
    0.97031231
    0.94796904
    0.94753789
    0.88092269
],
:ARsq => [
    0.71476928
    0.99035145
    0.77630678
    0.83771737
    0.82891994
    0.92517856
    0.96041642
    0.92195356
    0.92130684
    0.86391165
],
:Kel => [
    0.0033847439
    0.014106315
    0.0032914304
    0.0076953442
    0.0068133279
    0.0076922807
    0.012458956
    0.0089300798
    0.0056458649
    0.017189737
],
:HL => [
    204.78571
    49.137367
    210.59148
    90.073577
    101.73401
    90.109450057666
    55.634451
    77.619371
    122.77077
    40.323315
],
:Clast_pred => [
    117.30578
    82.53669
    66.931057
    100.76793
    105.29832
    71.939942
    61.172702
    75.604277
    93.761762
    38.810857
],
:AUCinf => [
    42925.019
    16154.93
    26026.183
    22004.078
    25820.275
    16001.76
    11688.953
    15446.21
    24865.246
    8075.3242
],
:AUCpct => [
    77.669383
    37.405019
    79.26492
    57.65405
    62.969953
    56.463551
    39.861391
    53.964925
    66.559429
    30.394194
],
:MRTlast => [
    34.801023
    29.538786
    34.472406
    33.69408
    32.964438
    32.58076
    31.267574
    33.826053
    33.386807
    27.556657
],
:MRTinf => [
    293.16224
    71.937917
    305.04073
    130.69968
    149.96684
    128.24114
    79.498252
    114.8571
    176.97811
    58.746446
],
:Clinf => [
    0.0023296437
    0.0061900608
    0.0038422846
    0.0045446122
    0.0038729255
    0.0062493127
    0.0085550864
    0.0064740799
    0.0040216775
    0.012383404
],
:Vzinf => [
    0.68827768
    0.43881487
    1.1673601
    0.59056646
    0.56843374
    0.8124135
    0.68666158
    0.72497447
    0.71232266
    0.72039519
],
)


# Linear up Log down, Dose 100, Dosetime 0.25, tau 9

refdict2 = Dict(
:Cmax => [
190.869
261.177
105.345
208.542
169.334
154.648
153.254
138.327
167.347
125.482
],
:Tmax => [
1
1
1.5
1
4
2.5
2.5
4
3
2
],
:Cdose => [
121.239
62.222
49.849
52.421
0
57.882
19.95
22.724
105.438
13.634
],
:Clast => [
112.846
85.241
67.901
97.625
110.778
69.501
58.051
74.437
93.44
42.191
],
:AUClast => [
9566.59680869131
10054.28647805950
5392.45721941379
9297.09633445033
9519.18087436122
6948.98562111745
6988.77263241364
7058.81896352039
8302.36808633358
5486.83888944199
],
:AUMCtau => [
5477.20423544297
8367.57088170951
3455.34643479800
6014.64604481587
6609.78830163090
5064.72384740413
4976.96365993911
2863.00517022791
5386.88322025614
4713.47970846693
],
:AUCall => [
9566.59680869131
10054.28647805950
5392.45721941379
9297.09633445033
9519.18087436122
6948.98562111745
6988.77263241364
7058.81896352039
8302.36808633358
5486.83888944199
],
:Rsq => [
0.786076957
0.992763591
0.81358898
0.918858685
0.853359952
0.95011904
0.970312315
0.94796904
0.947537895
0.88092269

],
:ARsq => [
0.714769276
0.990351454
0.776306776
0.83771737
0.828919944
0.92517856
0.96041642
0.92195356
0.921306842
0.863911645
],
# LZint
:LZint => [
5.00848559255328
5.42889759540296
4.44064607555325
5.16688496904739
5.14735707974283
4.82967584017057
5.01074587961482
4.96847859724365
4.94725938794774
4.89636108788302
],
:Kel => [
0.00338474394000776
0.01410631494324980
0.00329143037249282
0.00769534422298109
0.00681332791154901
0.00769228066663777
0.01245895597676470
0.00893007980967252
0.00564586491870971
0.01718973683041960
],
:HL => [
204.785706938398
49.137367437811
210.591476080649
90.073577019460
101.734011566509
90.109450057666
55.634451382012
77.619371308325
122.770769499451
40.323315440951
],
:Clast_pred => [
117.3057799
82.53668981
66.93105694
100.7679335
105.2983206
71.93994201
61.17270231
75.60427664
93.76176158
38.81085735
],
:AUCinf => [
42906.1941313004
16097.0411126277
26022.0900281352
21983.3384532182
25778.1957695968
15984.1473646863
11648.1518057779
15394.3547690766
24852.5337997128
7941.2685538530
],
:AUCpct => [
77.7034598328254
37.5395365663056
79.2773861992505
57.7084419901233
63.0727419426760
56.5257660444885
40.0010169085628
54.1467046238288
66.5934743183831
30.9072744205350
],
:MRTtauinf => [
299.791671096989
74.654997085457
305.919973652938
143.538421744963
173.022067431888
124.653434795141
92.735873637166
175.461862330056
178.810514188399
69.516339753006
],
:Cltau => [
0.078847213948573
0.054590500813083
0.132511872732088
0.074823364534525
0.076283206573122
0.089747243392665
0.092646906460213
0.130442680913677
0.081991954283052
0.103060243120434
],
:Vztau => [
23.2948829648816
3.8699335037324
40.2596615257358
9.7231991664617
11.1961742577834
11.6671826317919
7.4361693413954
14.6071125559694
14.5224789228203
5.9954520617241
],
:AUCtau => [
1268.27563070553
1831.82052757491
754.64936037981
1336.48093242129
1310.90451610924
1114.24035123260
1079.36685444479
766.62024499617
1219.63186357018
970.30626915116
],
)

################################################################################
# Linear Trapezoidal with Linear Interpolation, Dose 120, Dosetime 0.0, tau 12

refdict3 = Dict(
:Cmax => [190.869
261.177
105.345
208.542
169.334
154.648
153.254
138.327
167.347
125.482],

:Tmax => [1
1
1.5
1
4
2.5
2.5
4
3
2],

:Cdose => [0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0],

:Clast => [112.846
85.241
67.901
97.625
110.778
69.501
58.051
74.437
93.44
42.191],

:AUClast => [9585.4218
10112.176
5396.5498
9317.8358
9561.26
6966.598
7029.5735
7110.6745
8315.0803
5620.8945],

:AUMCtau => [9984.8168
14630.0690
6024.4953
10299.7210
11466.1230
8467.3568
9003.0193
6457.0058
10095.8180
8367.3005],

:AUCall => [9585.4218
10112.1760
5396.5498
9317.8358
9561.2600
6966.5980
7029.5735
7110.6745
8315.0803
5620.8945],

:Rsq => [0.78607696
0.99276359
0.81358898
0.91885869
0.86367664
0.95011904
0.97031231
0.94796904
0.94753789
0.87969895],

:ARsq => [0.71476928
0.99035145
0.77630678
0.83771737
0.84420187
0.92517856
0.96041642
0.92195356
0.92130684
0.86766884],

# LZint
:LZint => [5.0084856
5.4288976
4.4406461
5.166885
5.1496027
4.8296758
5.0107459
4.9684786
4.9472594
4.8651403],

:Kel => [0.0033847439
0.0141063150
0.0032914304
0.0076953442
0.0068579883
0.0076922807
0.0124589560
0.0089300798
0.0056458649
0.0165437520],

:HL => [204.785706938398
49.1373674378108
210.591476080649
90.0735770194602
101.071502239954
90.109450057666
55.6344513820121
77.6193713083247
122.770769499451
41.8978220179993],

:Clast_pred => [117.30578
82.53669
66.931057
100.76793
105.19623
71.939942
61.172702
75.604277
93.761762
39.408841],

:AUCinf => [42925.019
16154.93
26026.183
22004.078
25714.393
16001.76
11688.953
15446.21
24865.246
8171.1624],

:AUCpct => [77.669383
37.405019
79.26492
57.6540502829908
62.817478
56.463551
39.861391
53.964925
66.559429
31.210589],

:MRTtauinf => [302.40303
75.590599
312.72083
148.34069
172.0933
130.19061
91.908297
161.57402
176.30461
70.260736],

#Cltau, CLss
:Cltau => [0.07185191
0.050414459
0.12240579
0.070132959
0.06902661
0.085106504
0.083532913
0.10859036
0.073251565
0.092756742],

:Vztau => [21.228167
3.5738929
37.18924
9.113687
10.06514
11.063884
6.7046479
12.160066
12.974374
5.6067536],

:AUCtau => [1670.1018
2380.2695
980.34575
1711.0358
1738.46
1409.998
1436.5595
1105.0705
1638.1903
1293.7065],
)

################################################################################
#4 Log trapezoidal ATM, Dose 120, Dosetime 0, tau 12

refdict4 = Dict(
:Cmax => [
190.869
261.177
105.345
208.542
169.334
154.648
153.254
138.327
167.347
125.482
],
:Tmax => [
1
1
1.5
1
4
2.5
2.5
4
3
2
],
:Cdose => [
0
0
0
0
0
0
0
0
0
0
],
:Clast => [
112.846
85.241
67.901
97.625
110.778
69.501
58.051
74.437
93.44
42.191
],
:AUClast => [
9572.8582
10054.0367665966
5391.5322
9296.2179
9518.6531
6948.5757
6987.0645
7064.7816
8298.9634
5485.6538
],
:AUMCtau => [
9973.8062
14631.1197073321
6022.9286
10307.954
11473.081
8471.0956
8982.0378
6271.7444
10040.829690586
8361.7894
],
:AUCall => [
9572.8582
10054.0367665966
5391.5322
9296.2179
9518.6531
6948.5757
6987.0645
7064.7816
8298.9634
5485.6538
],
:Rsq => [
0.78607696
0.99276359
0.81358898
0.91885869
0.85335995
0.95011904
0.97031231
0.94796904
0.94753789
0.88092269
],
:ARsq => [
0.71476928
0.99035145
0.77630678
0.83771737
0.82891994
0.92517856
0.96041642
0.92195356
0.92130684
0.86391165
],
# LZint
:LZint => [
5.0084856
5.4288976
4.4406461
5.166885
5.1473571
4.8296758
5.0107459
4.9684786
4.9472594
4.8963611
],
:Kel => [
0.003384744
0.014106315
0.00329143
0.007695344
0.006813328
0.007692281
0.012458956
0.00893008
0.005645865
0.017189737
],
:HL => [
204.78571
49.137367
210.59148
90.073577
101.73401
90.109450057666
55.634451
77.619371
122.77077
40.323315
],
:Clast_pred => [
117.30578
82.53669
66.931057
100.76793
105.29832
71.939942
61.172702
75.604277
93.761762
38.810857
],
:AUCinf => [
42912.456
16096.791
26021.165
21982.4599914207
25777.668
15983.737
11646.444
15400.317
24849.129
7940.0834
],
:AUCpct => [
77.692122
37.540119
79.280204
57.710748
63.074033
56.527216
40.006884
54.12574
66.602599
30.911888
],
:MRTtauinf => [
302.63508
75.323724
313.06798
148.31081
172.5577
130.22554
91.866692
164.91799
176.98523
68.167555
],
:Cltau => [
0.071927102
0.050429351
0.12256044
0.070184147
0.069035447
0.0852177496596485
0.08379761
0.11110872
0.073575577
0.092819834
],
:Vztau => [
21.250382
3.5749486
37.236223
9.1203389
10.132412
11.078346
6.7258934
12.442074
13.031764
5.399724
],
:AUCtau => [
1668.3558
2379.5666
979.10878
1709.7878
1738.2375
1408.1573
1432.0218
1080.0233
1630.976
1292.8271
],
)


################################################################################

urefdict = Dict{Symbol, Float64}(
#:N_Samples => 5,
#:Dose => 100,
:Rsq => 0.90549162,
:ARsq => 0.81098324, #Rsq_adjusted
#:Corr_XY => -0.95157323,
#:No_points_lambda_z => 3,
:Kel => 0.13445441, #Lambda_z
#:Lambda_z_intercept => 0.79280975,
#:Lambda_z_lower => 4,
#:Lambda_z_upper => 15,
:HL => 5.1552579, #HL_Lambda_z
#:Span => 2.1337439,
#:Tlag => 0,
:Tmax => 1.5, #Tmax_Rate
:Maxrate => 4, #Max_Rate
#:Mid_Pt_last => 15,
:Rlast => 0.33333333, #Rate_last
#:Rate_last_pred => 0.2940497,
:AUClast => 17.125, #AURC_last
#:AURC_last_D => 0.17125,
:Vol => 11, #Vol_UR
:AR => 16, #Amount_Recovered
:Prec => 16, #Percent_Recovered
:AUCall => 17.125, #AURC_all
:AUCinf => 19.604155, #AURC_INF_obs
:AUCpct => 12.646069, #AURC_%Extrap_obs
#:AURC_INF_pred => 19.311984,
#:AURC_%Extrap_pred => 11.324493,
)
