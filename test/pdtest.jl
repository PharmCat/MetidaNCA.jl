
@testset "  PD                                                       " begin
    @test_nowarn MetidaNCA.pdimport(pddata, :time, :obs, :subj; bl = 0, th = NaN)
end
