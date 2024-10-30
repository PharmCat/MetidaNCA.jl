
@testset "  #6 Pharmacodynamics data; Linear-trapezoidal rule        " begin
    io = IOBuffer();

    pd =  MetidaNCA.pdimport(pddata, :time, :obs, :subj; bl = 1.5, th = 5.0)

    @test_nowarn MetidaNCA.pkplot(pd)
    @test_nowarn MetidaNCA.pkplot(pd[1], drawth = true, drawbl = true)

    pd_res = MetidaNCA.nca!(pd[1])
    pd_rds = MetidaNCA.nca!(pd)

    @test last(pd[1].time) - first(pd[1].time) == pd_res[:TABL] + pd_res[:TBBL] == pd_res[:TATH] + pd_res[:TBTH]

    nca_res = MetidaNCA.nca(pddata, :time, :obs, :subj)[1]
    pd_res  = MetidaNCA.nca(pddata, :time, :obs, :subj, type = :pd, bl = 0.0, th = 0.0)[1]

    @test  pd_res[:AUCABL] == pd_res[:AUCATH] == nca_res[:AUClast]

    @test MetidaNCA.getbl(pd[1]) ≈ 1.5
    @test MetidaNCA.getth(pd[1]) ≈ 5.0
    MetidaNCA.setbl!(pd, 2)
    MetidaNCA.setth!(pd, 6)
    @test MetidaNCA.getbl(pd[1]) ≈ 2.0
    @test MetidaNCA.getth(pd[1]) ≈ 6.0
    MetidaNCA.setbl!(pd, 3, 1)
    MetidaNCA.setth!(pd, 4, 1)
    @test MetidaNCA.getbl(pd[1]) ≈ 3.0
    @test MetidaNCA.getth(pd[1]) ≈ 4.0
    MetidaNCA.setbl!(pd, 2, [1])
    MetidaNCA.setth!(pd, 1, [1])
    @test MetidaNCA.getbl(pd[1]) ≈ 2.0
    @test MetidaNCA.getth(pd[1]) ≈ 1.0
    MetidaNCA.setbl!(pd, 0, Dict(:subj => 1))
    MetidaNCA.setth!(pd, 0, Dict(:subj => 1))
    @test MetidaNCA.getbl(pd[1]) ≈ 0.0
    @test MetidaNCA.getth(pd[1]) ≈ 0.0

    @test_throws ErrorException MetidaNCA.setbl!(pd, NaN)
    @test_throws ErrorException MetidaNCA.setth!(pd, NaN)

    pd =  MetidaNCA.pdimport(pddata, :time, :obs; bl = 3.0, th = 1.5, id = Dict(:subj => 1))
    pd_rds = MetidaNCA.nca!(pd)

    @test  pd_rds[:Tmax] ≈ 5.0 atol=1E-6
    @test  pd_rds[:Rmax] ≈ 8.0 atol=1E-6

    @test  pd_rds[:AUCABL] ≈ 7.3857143 atol=1E-6
    @test  pd_rds[:AUCBBL] ≈ 8.7357143 atol=1E-6
    @test  pd_rds[:AUCATH] ≈ 13.959524 atol=1E-6
    @test  pd_rds[:AUCBTH] ≈ 1.8095238 atol=1E-6
    @test  pd_rds[:AUCBTW] ≈ 6.926190 atol=1E-6
    @test  pd_rds[:TABL] ≈ 3.4809524 atol=1E-6
    @test  pd_rds[:TBBL] ≈ 5.5190476 atol=1E-6
    @test  pd_rds[:TATH] ≈ 5.7619048 atol=1E-6
    @test  pd_rds[:TBTH] ≈ 3.2380952 atol=1E-6
    @test  pd_rds[:AUCNETB] ≈ -1.35 atol=1E-2
    @test  pd_rds[:AUCNETT] ≈ 12.15 atol=1E-2
    @test  pd_rds[:TIMEBTW] ≈ 2.2809524 atol=1E-6

    pd =  MetidaNCA.pdimport(pddata, :time, :obs; bl = 1.5, th = 3.0, id = Dict(:subj => 1))
    pd_rds = MetidaNCA.nca!(pd)

    @test  pd_rds[:AUCATH] ≈ 7.3857143 atol=1E-6
    @test  pd_rds[:AUCBTH] ≈ 8.7357143 atol=1E-6
    @test  pd_rds[:AUCABL] ≈ 13.959524 atol=1E-6
    @test  pd_rds[:AUCBBL] ≈ 1.8095238 atol=1E-6
    @test  pd_rds[:AUCBTW] ≈ 6.5738095 atol=1E-6
    @test  pd_rds[:TATH] ≈ 3.4809524 atol=1E-6
    @test  pd_rds[:TBTH] ≈ 5.5190476 atol=1E-6
    @test  pd_rds[:TABL] ≈ 5.7619048 atol=1E-6
    @test  pd_rds[:TBBL] ≈ 3.2380952 atol=1E-6
    @test  pd_rds[:AUCNETT] ≈ -1.35 atol=1E-2
    @test  pd_rds[:AUCNETB] ≈ 12.15 atol=1E-2
    @test  pd_rds[:TIMEBTW] ≈ 2.2809524 atol=1E-6
    #
    pd_rds = MetidaNCA.nca!(pd; calcm = :luld)
end


@testset "  #7 Pharmacodynamics data; Linear-trapezoidal rule; Tau   " begin

    io = IOBuffer();

    pd =  MetidaNCA.pdimport(pddata, :time, :obs, :subj; bl = 1.5, th = 5.0)

    dt = MetidaNCA.DoseTime(dose = 100, time = 0.0, tau = 9.0)
    MetidaNCA.setdosetime!(pd, dt)
    pd_rds = MetidaNCA.nca!(pd)

    @test  pd_rds[1,:AUCATH] ≈ pd_rds[1,:AUCATHtau] atol=1E-6
    @test  pd_rds[1,:AUCBTH] ≈ pd_rds[1,:AUCBTHtau] atol=1E-6
    @test  pd_rds[1,:AUCABL] ≈ pd_rds[1,:AUCABLtau] atol=1E-6
    @test  pd_rds[1,:AUCBBL] ≈ pd_rds[1,:AUCBBLtau] atol=1E-6
    @test  pd_rds[1,:AUCBTW] ≈ pd_rds[1,:AUCBTWtau] atol=1E-6
    @test  pd_rds[1,:TATH] ≈ pd_rds[1,:TATHtau] atol=1E-6
    @test  pd_rds[1,:TBTH] ≈ pd_rds[1,:TBTHtau] atol=1E-6
    @test  pd_rds[1,:TABL] ≈ pd_rds[1,:TABLtau] atol=1E-6
    @test  pd_rds[1,:TBBL] ≈ pd_rds[1,:TBBLtau] atol=1E-6


    dt = MetidaNCA.DoseTime(dose = 100, time = 1.0, tau = 7.0)
    MetidaNCA.setdosetime!(pd, dt)
    pd_rds = MetidaNCA.nca!(pd)
    @test  pd_rds[1,:AUCATH] ≈ pd_rds[1,:AUCATHtau] atol=1E-6
    @test  pd_rds[1,:AUCBTH] - 8.5 ≈ pd_rds[1,:AUCBTHtau] atol=1E-6
    @test  pd_rds[1,:AUCABL] ≈ pd_rds[1,:AUCABLtau] atol=1E-6
    @test  pd_rds[1,:AUCBBL] - 1.5 ≈ pd_rds[1,:AUCBBLtau] atol=1E-6
    @test  pd_rds[1,:AUCBTW] ≈ pd_rds[1,:AUCBTWtau] atol=1E-6
    @test  pd_rds[1,:TATH] ≈ pd_rds[1,:TATHtau] atol=1E-6
    @test  pd_rds[1,:TBTH] - 2 ≈ pd_rds[1,:TBTHtau] atol=1E-6
    @test  pd_rds[1,:TABL] ≈ pd_rds[1,:TABLtau] atol=1E-6
    @test  pd_rds[1,:TBBL] - 2 ≈ pd_rds[1,:TBBLtau] atol=1E-6


    dt = MetidaNCA.DoseTime(dose = 100, time = 0.5, tau = 8.0)
    MetidaNCA.setdosetime!(pd, dt)
    pd_rds = MetidaNCA.nca!(pd)
    @test  pd_rds[1,:AUCATH] ≈ pd_rds[1,:AUCATHtau] atol=1E-6
    @test  pd_rds[1,:AUCBTH] - 4.375 ≈ pd_rds[1,:AUCBTHtau] atol=1E-6
    @test  pd_rds[1,:AUCABL] ≈ pd_rds[1,:AUCABLtau] atol=1E-6
    @test  pd_rds[1,:AUCBBL] - 0.875 ≈ pd_rds[1,:AUCBBLtau] atol=1E-6
    @test  pd_rds[1,:AUCBTW] ≈ pd_rds[1,:AUCBTWtau] atol=1E-6
    @test  pd_rds[1,:TATH] ≈ pd_rds[1,:TATHtau] atol=1E-6
    @test  pd_rds[1,:TBTH] - 1 ≈ pd_rds[1,:TBTHtau] atol=1E-6
    @test  pd_rds[1,:TABL] ≈ pd_rds[1,:TABLtau] atol=1E-6
    @test  pd_rds[1,:TBBL] - 1 ≈ pd_rds[1,:TBBLtau] atol=1E-6

end
