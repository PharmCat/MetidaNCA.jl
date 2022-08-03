
@testset "  PD                                                       " begin
    pd =  MetidaNCA.pdimport(pddata, :time, :obs, :subj; bl = 1.5, th = 5.0)

    pd_res = MetidaNCA.nca!(pd[1])
    pd_rds = MetidaNCA.nca!(pd)

    @test last(pd[1].time) - first(pd[1].time) == pd_res[:TABL] + pd_res[:TBBL] == pd_res[:TATH] + pd_res[:TBTH]

    nca_res = MetidaNCA.nca(pddata, :time, :obs, :subj)[1]
    pd_res  = MetidaNCA.nca(pddata, :time, :obs, :subj, type = :pd, bl = 0.0, th = 0.0)[1]

    @test  pd_res[:AUCABL] == pd_res[:AUCATH] == nca_res[:AUClast]


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

end
