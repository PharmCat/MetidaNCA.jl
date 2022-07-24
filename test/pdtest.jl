
@testset "  PD                                                       " begin
    pd =  MetidaNCA.pdimport(pddata, :time, :obs, :subj; bl = 1.5, th = 5.0)

    pd_res = MetidaNCA.nca!(pd[1])
    pd_rds = MetidaNCA.nca!(pd)

    @test last(pd[1].time) - first(pd[1].time) == pd_res[:TABL] + pd_res[:TBBL] == pd_res[:TATH] + pd_res[:TBTH]

    nca_res = MetidaNCA.nca(pddata, :time, :obs, :subj)[1]
    pd_res  = MetidaNCA.nca(pddata, :time, :obs, :subj, type = :pd, bl = 0.0, th = 0.0)[1]

    @test  pd_res[:AUCABL] == pd_res[:AUCATH] == nca_res[:AUClast]

end
