
@testset "  #5 Urine data; Linear-trapezoidal rule                   " begin
    upkds = MetidaNCA.upkimport(upkdata, :st, :et, :conc, :vol, :subj; dosetime =  MetidaNCA.DoseTime(dose = 100.0, time = 0))
    unca  = MetidaNCA.nca!(upkds)
    for k in keys(urefdict)
        @test  unca[1, k]  ≈ urefdict[k] atol=1E-4
    end

    upkdatac = deepcopy(upkdata)
    upkdatac.st = float.(upkdatac.st)
    upkdatac[1, :st] = NaN
    @test_throws ErrorException MetidaNCA.upkimport(upkdatac, :st, :et, :conc, :vol, :subj; dosetime =  MetidaNCA.DoseTime(dose = 100, time = 0))

    upkdatac = deepcopy(upkdata)
    upkdatac.et = float.(upkdatac.et)
    upkdatac[1, :et] = NaN
    @test_throws ErrorException MetidaNCA.upkimport(upkdatac, :st, :et, :conc, :vol, :subj; dosetime =  MetidaNCA.DoseTime(dose = 100, time = 0))

    upkdatac = deepcopy(upkdata)
    upkdatac.et = float.(upkdatac.et)
    upkdatac[1, :et] = 1.5
    @test_throws ErrorException MetidaNCA.upkimport(upkdatac, :st, :et, :conc, :vol, :subj; dosetime =  MetidaNCA.DoseTime(dose = 100, time = 0))
end
