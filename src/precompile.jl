
PrecompileTools.@compile_workload begin
    data = metida_table([0.,1.,2.,3.,4.,2.,1.,0.], [0.,1.,2.,3.,4.,5.,6.,7.], names = (:conc, :time))
    pki  = pkimport(data, :time, :conc; dosetime = MetidaNCA.DoseTime(dose = 100, time = 0, tau = 5.5))
    pkr  =  nca!(pki)
end