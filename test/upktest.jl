
@testset "  Urine PK                                                 " begin
    @test_nowarn upkds = MetidaNCA.upkimport(upkdata, :st, :et, :conc, :vol, :subj)
end
