
@testset "  Urine PK                                                 " begin
    upkds = MetidaNCA.upkimport(upkdata, :st, :et, :conc, :vol, :subj; dosetime =  MetidaNCA.DoseTime(dose = 100))
    unca  = MetidaNCA.nca!(upkds)
    @test unca[1, :Maxrate] == 4.0
    @test unca[1, :AUCall]  == 17.1250
    @test unca[1, :Tmax]    == 1.5
    @test unca[1, :AR]      == 16.0
    @test unca[1, :Vol]     == 11.0
    @test unca[1, :Prec]    == 16.0

    @test unca[1, :HL]      ≈ 5.15525788157323 atol=1E-6
    @test unca[1, :ARsq]    ≈ 0.8109832410336681 atol=1E-6
    @test unca[1, :Kel]     ≈ 0.13445441459631066 atol=1E-6
    @test unca[1, :AUCinf]  ≈ 19.60415499341648 atol=1E-6

    io = IOBuffer();

    @test_nowarn dsnca = MetidaNCA.nca!(upkds, verbose = 2, io = io)
end
