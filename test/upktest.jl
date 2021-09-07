
@testset "  Urine PK                                                 " begin
    upkds = MetidaNCA.upkimport(upkdata, :st, :et, :conc, :vol, :subj)
    unca  = MetidaNCA.nca!(upkds)
    unca[1, :Maxrate] == 4.0
    unca[1, :AUCall]  == 17.1250
    unca[1, :Tmax]    == 1.5
    unca[1, :AR]      == 16.0
    unca[1, :Vol]     == 11.0
end
