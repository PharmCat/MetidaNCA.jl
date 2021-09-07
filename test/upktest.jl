
@testset "  Urine PK                                                 " begin
    upkds = MetidaNCA.upkimport(upkdata, :st, :et, :conc, :vol, :subj)
    @test_nowarn MetidaNCA.nca!(upkds)
end
