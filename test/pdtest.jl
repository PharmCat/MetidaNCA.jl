
@testset "  PD                                                       " begin
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

    pd =  MetidaNCA.pdimport(pddata, :time, :obs; bl = 1.5, th = 5.0, id = Dict(:subj => 1))
    pd_rds = MetidaNCA.nca!(pd)

    # Not validated
    @test  pd_rds[:Tmax] ≈ 5.0 atol=1E-6
    @test  pd_rds[:Rmax] ≈ 8.0 atol=1E-6

    @test  pd_rds[:AUCABL] ≈ 13.959523809523809 atol=1E-6
    @test  pd_rds[:AUCBBL] ≈ 1.8095238095238095 atol=1E-6
    @test  pd_rds[:AUCATH] ≈ 2.2261904761904767 atol=1E-6
    @test  pd_rds[:AUCBTH] ≈ 13.542857142857143 atol=1E-6
    @test  pd_rds[:AUCBTW] ≈ 11.733333333333333 atol=1E-6

    #
    pd_rds = MetidaNCA.nca!(pd; calcm = :luld)

    

end
