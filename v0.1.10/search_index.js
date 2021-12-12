var documenterSearchIndex = {"docs":
[{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/#pkimport","page":"API","title":"pkimport","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"MetidaNCA.pkimport","category":"page"},{"location":"api/#MetidaNCA.pkimport","page":"API","title":"MetidaNCA.pkimport","text":"pkimport(data, time, conc, sort; kelauto = true,  elimrange = ElimRange(), dosetime = DoseTime())\n\nImport PK data from table data.\n\ntime - time column;\nconc - concentration column;\nsort - subject sorting columns.\n\nkeywords:\n\nkelauto - if true auto range settings, if false used kelstart/kelend from elimrange;\nelimrange - set elimination range settings;\ndosetime - set dose and dose time, by default dosetime = 0, dose is NaN.\n\n\n\n\n\npkimport(data, time, conc; kelauto = true,  elimrange = ElimRange(), dosetime = DoseTime())\n\nImport PK data from tabular data data, time - time column, conc - concentration column.\n\n\n\n\n\npkimport(time, conc; kelauto = true,  elimrange = ElimRange(), dosetime = DoseTime())\n\nImport PK data from time vector time and concentration vector conc.\n\n\n\n\n\n","category":"function"},{"location":"api/#nca!","page":"API","title":"nca!","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"MetidaNCA.nca!","category":"page"},{"location":"api/#MetidaNCA.nca!","page":"API","title":"MetidaNCA.nca!","text":"nca!(data::PKSubject{T,O}; adm = :ev, calcm = :lint, intpm = nothing, limitrule = nothing, verbose = 0, warn = true, io::IO = stdout, modify! = nothing) where T where O\n\nadm - administration:\n:ev - extra vascular;\n:iv - intravascular bolus;\ncalcm - AUC/AUMC calculation method:\n:lint - linear trapezoidal;\n:logt - log-trapezoidal after Tmax;\n:luld - linar up log down;\n:luldt - linear up log down after Tmax;\nintpm - interpolation method:\n:lint - linear trapezoidal;\n:logt - log-trapezoidal after Tmax;\n:luld - linar up log down;\n:luldt - linear up log down after Tmax;\nlimitrule use limitrule for data;\nverbose - print to io, 1: partial areas table, 2: 1, and results;\nwarn - show warnings;\nio - output stream;\nmodify! - function to modify output paramaters, call modify!(data, result) if difined.\n\n\n\n\n\nnca!(data::DataSet{T1}; adm = :ev, calcm = :lint, intpm = nothing, verbose = 0, warn = true, io::IO = stdout) where T1 <: PKSubject{T,O}  where T  where O\n\nNon-compartmental (NCA) analysis of pharmacokinetic (PK) data.\n\n\n\n\n\n","category":"function"},{"location":"api/#nca","page":"API","title":"nca","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"MetidaNCA.nca","category":"page"},{"location":"api/#MetidaNCA.nca","page":"API","title":"MetidaNCA.nca","text":"nca(args...; kelauto = true,  elimrange = ElimRange(), dosetime = DoseTime(), kwargs...)\n\nnca(data, time, conc, sort; kelauto = true,  elimrange = ElimRange(), dosetime = DoseTime(), kwargs...)\n\nnca(data, time, conc; kelauto = true,  elimrange = ElimRange(), dosetime = DoseTime(), kwargs...)\n\nnca(time, conc; kelauto = true,  elimrange = ElimRange(), dosetime = DoseTime(), kwargs...)\n\nImport data and perform NCA analysis.\n\nSyntax simillar to pkimport\n\nApplicable kwargs see  nca!.\n\n\n\n\n\n","category":"function"},{"location":"api/#DoseTime","page":"API","title":"DoseTime","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"MetidaNCA.DoseTime","category":"page"},{"location":"api/#MetidaNCA.DoseTime","page":"API","title":"MetidaNCA.DoseTime","text":"DoseTime(dose::D, time::T, tau::TAU) where D <: Number where T <: Number where TAU <: Number\n\nDose settings.\n\ndose - dose;\ntime - dose time;\ntau - tau (τ);\n\nDose time set 0 by default.\n\n\n\n\n\n","category":"type"},{"location":"api/#setdosetime!","page":"API","title":"setdosetime!","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"MetidaNCA.setdosetime!","category":"page"},{"location":"api/#MetidaNCA.setdosetime!","page":"API","title":"MetidaNCA.setdosetime!","text":"setdosetime!(data::T, dosetime::DoseTime) where T <: PKSubject\n\nSet dose time dosetime for subject data.\n\n\n\n\n\nsetdosetime!(data::DataSet{T}, dosetime::DoseTime, ind::Int) where T <: PKSubject\n\nind - index in DataSet.\n\n\n\n\n\nsetdosetime!(data::DataSet{T}, dosetime::DoseTime, inds::Union{Vector{Int}, UnitRange{Int}, Tuple{Vararg{Int}}}) where T <: PKSubject\n\ninds - indexes in DataSet.\n\n\n\n\n\nsetdosetime!(data::DataSet{T}, dosetime::DoseTime) where T <: PKSubject\n\nFor all subjects in DataSet.\n\n\n\n\n\nsetdosetime!(data::DataSet{T}, dosetime::DoseTime, sort::Dict) where T <: PKSubject\n\nSet dose time dosetime for subjects if sort ⊆ subject's id.\n\n\n\n\n\n","category":"function"},{"location":"api/#setkelauto!","page":"API","title":"setkelauto!","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"MetidaNCA.setkelauto!","category":"page"},{"location":"api/#MetidaNCA.setkelauto!","page":"API","title":"MetidaNCA.setkelauto!","text":"setkelauto!(data::T, kelauto::Bool) where T <: PKSubject\n\nSet range for elimination parameters calculation for subject.\n\ndata     - PK subject;\nkelauto  - value.\n\n\n\n\n\nsetkelauto!(data::DataSet{T}, kelauto::Bool, ind::Int) where T <: PKSubject\n\n\n\n\n\nsetkelauto!(data::DataSet{T}, kelauto::Bool, inds::Union{Vector{Int}, UnitRange{Int}, Tuple{Vararg{Int}}}) where T <: PKSubject\n\n\n\n\n\nsetkelauto!(data::DataSet{T}, kelauto::Bool) where T <: PKSubject\n\n\n\n\n\nsetkelauto!(data::DataSet{T}, kelauto::Bool, sort::Dict) where T <: PKSubject\n\n\n\n\n\n","category":"function"},{"location":"api/#ElimRange","page":"API","title":"ElimRange","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"MetidaNCA.ElimRange","category":"page"},{"location":"api/#MetidaNCA.ElimRange","page":"API","title":"MetidaNCA.ElimRange","text":"ElimRange(kelstart::Int, kelend::Int, kelexcl::Vector{Int})::ElimRange\n\nElimination settings for PK subject.\n\nkelstart - start point;\nkelend - end point;\nkelexcl - excluded points.\n\n\n\n\n\n","category":"type"},{"location":"api/#setkelrange!","page":"API","title":"setkelrange!","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"MetidaNCA.setkelrange!","category":"page"},{"location":"api/#MetidaNCA.setkelrange!","page":"API","title":"MetidaNCA.setkelrange!","text":"setkelrange!(data::T, range::ElimRange{:point}; kelauto = false) where T <: PKSubject\n\nSet range for subject data. Set kelauto if possible.\n\n\n\n\n\nsetdosetime!(data::DataSet{T}, dosetime::DoseTime, ind::Int) where T <: PKSubject\n\n\n\n\n\nsetkelrange!(data::DataSet{T}, range::ElimRange{:point}, inds::Union{Vector{Int}, UnitRange{Int}, Tuple{Vararg{Int}}}; kelauto = false)\n\n\n\n\n\nsetkelrange!(data::DataSet{T}, range::ElimRange{:point}; kelauto = false) where T <: PKSubject\n\n\n\n\n\nsetkelrange!(data::DataSet{T}, range::ElimRange{:point}, sort::Dict; kelauto = false) where T <: PKSubject\n\n\n\n\n\n","category":"function"},{"location":"api/#getdosetime","page":"API","title":"getdosetime","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"MetidaNCA.getdosetime","category":"page"},{"location":"api/#MetidaNCA.getdosetime","page":"API","title":"MetidaNCA.getdosetime","text":"getdosetime(data::T) where T <: PKSubject\n\nReturn dosetime.\n\n\n\n\n\n","category":"function"},{"location":"api/#getkelauto","page":"API","title":"getkelauto","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"MetidaNCA.getkelauto","category":"page"},{"location":"api/#MetidaNCA.getkelauto","page":"API","title":"MetidaNCA.getkelauto","text":"getkelauto!(data::T) where T <: PKSubject\n\n\n\n\n\n","category":"function"},{"location":"api/#getkelrange","page":"API","title":"getkelrange","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"MetidaNCA.getkelrange","category":"page"},{"location":"api/#MetidaNCA.getkelrange","page":"API","title":"MetidaNCA.getkelrange","text":"getkelrange(data::T) where T <: PKSubject\n\n\n\n\n\n","category":"function"},{"location":"api/#LimitRule","page":"API","title":"LimitRule","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"MetidaNCA.LimitRule","category":"page"},{"location":"api/#MetidaNCA.LimitRule","page":"API","title":"MetidaNCA.LimitRule","text":"LimitRule(lloq::T, btmax, atmax, nan, rm::Bool) where T <: Real\n\nLimitRule(;lloq = NaN, btmax = NaN, atmax = NaN, nan = NaN, rm::Bool = false)\n\nlloq - LLOQ - low limit of quantification;\nbtmax - value for points before Tmax;\natmat - values for points after Tmax;\nnan - values for replacing NaN;\nrm - if true, removee all NaN points.\n\nRule for PK subject.\n\nSTEP 1 (NaN step): replace all NaN and missing values with nan keyword value (if nan not NaN);\nSTEP 2 (LLOQ step): replace values below lloq with btmax value if this value befor Tmax or with atmax if this value after Tmax (if lloq not NaN);\nSTEP 3 (remove NaN): rm == true, then remove all NaN and missing values.\n\n\n\n\n\n","category":"type"},{"location":"api/#applylimitrule!","page":"API","title":"applylimitrule!","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"MetidaNCA.applylimitrule!","category":"page"},{"location":"api/#MetidaNCA.applylimitrule!","page":"API","title":"MetidaNCA.applylimitrule!","text":"applylimitrule!(data::PKSubject, rule::LimitRule)\n\nApply rule to PK subject .\n\nSTEP 1 (NaN step): replace all NaN and missing values with nan keyword value (if nan not NaN);\nSTEP 2 (LLOQ step): replace values below lloq with btmax value if this value befor Tmax or with atmax if this value after Tmax (if lloq not NaN);\nSTEP 3 (remove NaN): rm == true, then remove all NaN and missing values.\n\n\n\n\n\napplylimitrule!(f::Function, data::DataSet{T}, rule::LimitRule) where T <: PKSubject\n\n\n\n\n\napplylimitrule!(data::DataSet{T}, rule::LimitRule, ind::Int) where T <: PKSubject\n\n\n\n\n\napplylimitrule!(data::DataSet{T}, rule::LimitRule, inds::Union{Vector{Int}, UnitRange{Int}, Tuple{Vararg{Int}}}) where T <: PKSubject\n\n\n\n\n\napplylimitrule!(data::DataSet{T}, rule::LimitRule) where T <: PKSubject\n\n\n\n\n\napplylimitrule!(data::DataSet{T}, rule::LimitRule, sort::Dict) where T <: PKSubject\n\n\n\n\n\napplylimitrule!(time, obs, rule::LimitRule)\n\n\n\n\n\n","category":"function"},{"location":"api/#pkplot","page":"API","title":"pkplot","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"MetidaNCA.pkplot","category":"page"},{"location":"api/#MetidaNCA.pkplot","page":"API","title":"MetidaNCA.pkplot","text":"pkplot(subj; ls = false, elim = false, xticksn = :auto, yticksn = :auto, kwargs...)\n\nPlot for subject\n\nls - concentration in log scale;\nelim - draw elimination curve;\nxticksn - number of ticks on x axis;\nyticksn - number of ticks on y axis/\n\n\n\n\n\npkplot(data::DataSet{T};\ntypesort::Union{Nothing, Symbol, AbstractVector{Symbol}} = nothing,\npagesort::Union{Nothing, Symbol, AbstractVector{Symbol}} = nothing,\nsort::Union{Nothing, Dict{Symbol}} = nothing,\nkwargs...) where T <: AbstractSubject\n\nPK plot for subject set.\n\ntypesort - sort on page by this id key;\npagesort - different pages by this id key;\nsort - use only subjects if sort ⊆ subject id.\n\n\n\n\n\n","category":"function"},{"location":"details/#Details","page":"Details","title":"Details","text":"","category":"section"},{"location":"details/#Step-1","page":"Details","title":"Step 1","text":"","category":"section"},{"location":"details/","page":"Details","title":"Details","text":"Filter all values before dose time and NaN or missing values. If TAU set, calculate start and end timepoints for AUCtau.","category":"page"},{"location":"details/#Step-2","page":"Details","title":"Step 2","text":"","category":"section"},{"location":"details/","page":"Details","title":"Details","text":"Cmax, Tmax calculation.","category":"page"},{"location":"details/#Step-3","page":"Details","title":"Step 3","text":"","category":"section"},{"location":"details/","page":"Details","title":"Details","text":"Elimination parameters calculation.","category":"page"},{"location":"details/#Step-4","page":"Details","title":"Step 4","text":"","category":"section"},{"location":"details/","page":"Details","title":"Details","text":"Shift all time values by dose time.","category":"page"},{"location":"details/#Step-5","page":"Details","title":"Step 5","text":"","category":"section"},{"location":"details/","page":"Details","title":"Details","text":"Calculate dose concentration (Cdose).","category":"page"},{"location":"details/#Step-6","page":"Details","title":"Step 6","text":"","category":"section"},{"location":"details/","page":"Details","title":"Details","text":"Calculate areas.","category":"page"},{"location":"details/#Step-7","page":"Details","title":"Step 7","text":"","category":"section"},{"location":"details/","page":"Details","title":"Details","text":"Calculate steady-state parameters.","category":"page"},{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"ENV[\"GKSwstype\"] = \"nul\"","category":"page"},{"location":"examples/#Import","page":"Examples","title":"Import","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"using MetidaNCA, CSV, DataFrames;\n\npkdata2 = CSV.File(joinpath(dirname(pathof(MetidaNCA)), \"..\", \"test\", \"csv\",  \"pkdata2.csv\")) |> DataFrame\n\nds = MetidaNCA.pkimport(pkdata2, :Time, :Concentration, [:Subject, :Formulation]; dosetime = MetidaNCA.DoseTime(dose = 100, time = 0))\n\nsort!(ds, :Subject)","category":"page"},{"location":"examples/#NCA","page":"Examples","title":"NCA","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"dsnca = MetidaNCA.nca!(ds, adm = :ev, calcm = :lint)\n\ndsnca[:, :AUClast]","category":"page"},{"location":"examples/#Print-output","page":"Examples","title":"Print output","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"dsnca = MetidaNCA.nca!(ds[1], adm = :ev, calcm = :lint, verbose = 2);\n","category":"page"},{"location":"examples/#Plotting","page":"Examples","title":"Plotting","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"using Plots\n\np = MetidaNCA.pkplot(ds; typesort = :Subject, pagesort = nothing, sort = Dict(:Formulation => \"R\"))\n\npng(p, \"plot1.png\")\n\np = MetidaNCA.pkplot(ds; typesort = :Formulation, pagesort = nothing, legend = true)\n\npng(p, \"plot2.png\")\n\np = MetidaNCA.pkplot(ds; elim = true, ls = true)\n\npng(p[1], \"plot3.png\")\n\np = MetidaNCA.pkplot(ds; typesort = :Subject, pagesort = :Formulation)\n\npng(p[1], \"plot4.png\")","category":"page"},{"location":"examples/#Plot-1","page":"Examples","title":"Plot 1","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"(Image: )","category":"page"},{"location":"examples/#Plot-2","page":"Examples","title":"Plot 2","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"(Image: )","category":"page"},{"location":"examples/#Plot-3","page":"Examples","title":"Plot 3","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"(Image: )","category":"page"},{"location":"examples/#Plot-4","page":"Examples","title":"Plot 4","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"(Image: )","category":"page"},{"location":"examples/#Set-dose-time","page":"Examples","title":"Set dose time","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"dt = MetidaNCA.DoseTime(dose = 200, time = 0)\n\nMetidaNCA.setdosetime!(ds, dt, Dict(:Formulation => \"R\"))\n\ndsnca = MetidaNCA.nca!(ds)\n\ndsnca[:, :Dose]","category":"page"},{"location":"examples/#Set-range-for-elimination","page":"Examples","title":"Set range for elimination","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"kr =  MetidaNCA.ElimRange(kelstart = 4, kelend = 12, kelexcl = Int[5,6])\n\nMetidaNCA.setkelrange!(ds, kr, [1,2,3])\n\ndsnca = MetidaNCA.nca!(ds)\n\np = MetidaNCA.pkplot(ds[1]; elim = true)\n\npng(p, \"plot5.png\")\n\nMetidaNCA.getkeldata(ds[1])","category":"page"},{"location":"examples/#Plot-5","page":"Examples","title":"Plot 5","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"(Image: )","category":"page"},{"location":"examples/#Without-import","page":"Examples","title":"Without import","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"\ndsnca = MetidaNCA.nca(pkdata2, :Time, :Concentration, [:Subject, :Formulation]; dosetime = MetidaNCA.DoseTime(dose = 100, time = 0))\n\nsort!(dsnca, :Subject)\n\ndsnca[:, :AUClast]","category":"page"},{"location":"#MetidaNCA","page":"Home","title":"MetidaNCA","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = MetidaNCA","category":"page"},{"location":"","page":"Home","title":"Home","text":"Non-compartment pharmacokinetic analysis (NCA).","category":"page"},{"location":"","page":"Home","title":"Home","text":"*This program comes with absolutely no warranty. No liability is accepted for any loss and risk to public health resulting from use of this software.","category":"page"},{"location":"#NCA","page":"Home","title":"NCA","text":"","category":"section"},{"location":"#Validation","page":"Home","title":"Validation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Validation report: validation_report.pdf Appendix 2: Appendix2.1.pdf, Appendix2.2.pdf, Appendix2.3.pdf","category":"page"},{"location":"#Contents","page":"Home","title":"Contents","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\n        \"index.md\",\n        \"examples.md\",\n        \"parameters.md\",\n        \"api.md\",\n      ]\nDepth = 3","category":"page"},{"location":"#Reference","page":"Home","title":"Reference","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Makoid C, Vuchetich J, Banakar V (1996-1999), Basic Pharmacokinetics;\nGabrielsson and Weiner (1997), Pharmacokinetic and Pharmacodynamic Data Analysis: Concepts and Applications;\nGibaldi and Perrier (1982), Pharmacokinetics;\nWagner (1975), Fundamentals of Clinical Pharmacokinetics.","category":"page"},{"location":"parameters/#Parameters","page":"Parameters","title":"Parameters","text":"","category":"section"},{"location":"parameters/#Basic-parameters","page":"Parameters","title":"Basic parameters","text":"","category":"section"},{"location":"parameters/#Cmax:-Maximum-concentration","page":"Parameters","title":"Cmax: Maximum concentration","text":"","category":"section"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"Maximum concentration from dose time to dose time + tau (if tau > 0). Firs observation used.","category":"page"},{"location":"parameters/#Tmax:-Time-at-maximum-concentration","page":"Parameters","title":"Tmax: Time at maximum concentration","text":"","category":"section"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"Time at maximum concentration from dose time to dose time + tau (if tau > 0). Firs observation used.","category":"page"},{"location":"parameters/#Cdose","page":"Parameters","title":"Cdose","text":"","category":"section"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"Concentration at dose time.","category":"page"},{"location":"parameters/#AUC-/-AUMC","page":"Parameters","title":"AUC / AUMC","text":"","category":"section"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"Area under Curve / Area under the Moment Curve.","category":"page"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"AUC = sum_n=1^N AUC_n","category":"page"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"AUMC = sum_n=1^N AUMC_n","category":"page"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"Where AUCn/AUMCn- partial AUC/AUMC.","category":"page"},{"location":"parameters/#Linear-trapezoidal-rule","page":"Parameters","title":"Linear trapezoidal rule","text":"","category":"section"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"AUCmid_t_1^t_2 = delta t times fracC_1 + C_22","category":"page"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"AUMCmid_t_1^t_2 = delta t times fract_1 times C_1 + t_2 times C_22","category":"page"},{"location":"parameters/#Logarithmic-trapezoidal-rule","page":"Parameters","title":"Logarithmic trapezoidal rule","text":"","category":"section"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"AUCmid_t_1^t_2 =   delta t times frac C_2 - C_1ln(C_2C_1)","category":"page"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"AUMCmid_t_1^t_2 = delta t times fract_2 times C_2 - t_1 times C_1ln(C_2C_1) -  delta t^2 times frac C_2 - C_1ln(C_2C_1)^2","category":"page"},{"location":"parameters/#Interpolation","page":"Parameters","title":"Interpolation","text":"","category":"section"},{"location":"parameters/#Linear-interpolation-rule","page":"Parameters","title":"Linear interpolation rule","text":"","category":"section"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"C_x = C_1 + frac(t_x-t_1)times(C_2 - C_1)t_2 - t_1","category":"page"},{"location":"parameters/#Logarithmic-interpolation-rule","page":"Parameters","title":"Logarithmic interpolation rule","text":"","category":"section"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"C_x = expleft(ln(C_1) + frac(t_x-t_1)times(ln(C_2) - ln(C_1))t_2 - t_1right)","category":"page"},{"location":"parameters/#AUClast-/-AUMClast","page":"Parameters","title":"AUClast / AUMClast","text":"","category":"section"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"Area from dose time to last observed concentration (>0).","category":"page"},{"location":"parameters/#AUCall","page":"Parameters","title":"AUCall","text":"","category":"section"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"All values used to calculate AUC.","category":"page"},{"location":"parameters/#𝝺z-elimination-constant","page":"Parameters","title":"𝝺z - elimination constant","text":"","category":"section"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"Linear regression used for logarithmic transformed concentration data.","category":"page"},{"location":"parameters/#Half-Life;-T1/2","page":"Parameters","title":"Half-Life; T1/2","text":"","category":"section"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"HL = ln(2)  lambda_z","category":"page"},{"location":"parameters/#Rsq","page":"Parameters","title":"Rsq","text":"","category":"section"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"R²","category":"page"},{"location":"parameters/#ARsq","page":"Parameters","title":"ARsq","text":"","category":"section"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"Adjusted R²","category":"page"},{"location":"parameters/#If-Kel-calculated","page":"Parameters","title":"If Kel calculated","text":"","category":"section"},{"location":"parameters/#AUCinf","page":"Parameters","title":"AUCinf","text":"","category":"section"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"AUC_infty = AUC_last + fracC_lastlambda_z","category":"page"},{"location":"parameters/#AUMCinf","page":"Parameters","title":"AUMCinf","text":"","category":"section"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"AUMC_infty =  AUMC_last + fract_lasttimes C_lastlambda_z + fracC_lastlambda_z^2","category":"page"},{"location":"parameters/#AUCpct","page":"Parameters","title":"AUCpct","text":"","category":"section"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"AUCpct = (AUC_infty - AUC_last)  AUC_infty * 100 ","category":"page"},{"location":"parameters/#AUCinf_pred:-AUC-to-infinity-from-predicted-concentration","page":"Parameters","title":"AUCinf_pred: AUC to infinity from predicted concentration","text":"","category":"section"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"AUC_infty pred = AUC_last + fracC_last predlambda_z","category":"page"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"result[:AUCinfpred]     = result[:AUClast] + result[:Clastpred] / result[:Kel]","category":"page"},{"location":"parameters/#If-Dose-used","page":"Parameters","title":"If Dose used","text":"","category":"section"},{"location":"parameters/#Clearance","page":"Parameters","title":"Clearance","text":"","category":"section"},{"location":"parameters/#Cllast","page":"Parameters","title":"Cllast","text":"","category":"section"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"CL_last = Dose  AUC_last","category":"page"},{"location":"parameters/#Clinf","page":"Parameters","title":"Clinf","text":"","category":"section"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"CL_infty = Dose  AUC_infty","category":"page"},{"location":"parameters/#Steady-state-parameters-(If-Tau-used)","page":"Parameters","title":"Steady-state parameters (If Tau used)","text":"","category":"section"},{"location":"parameters/#AUCtau-/-AUMCtau","page":"Parameters","title":"AUCtau / AUMCtau","text":"","category":"section"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"Area from dose time to dose time + tau.","category":"page"},{"location":"parameters/#Ctau","page":"Parameters","title":"Ctau","text":"","category":"section"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"Concentration at τ time.","category":"page"},{"location":"parameters/#Ctaumin","page":"Parameters","title":"Ctaumin","text":"","category":"section"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"Minimum concentration from dose time to τ time.","category":"page"},{"location":"parameters/#Cavg","page":"Parameters","title":"Cavg","text":"","category":"section"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"C_avg = AUC_tau  tau","category":"page"},{"location":"parameters/#Fluc:-Fluctuation","page":"Parameters","title":"Fluc: Fluctuation","text":"","category":"section"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"Fluc = ( C_max - C_tau min )  C_avg * 100 ","category":"page"},{"location":"parameters/#Fluctau:-Fluctuation-Tau","page":"Parameters","title":"Fluctau: Fluctuation Tau","text":"","category":"section"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"Fluc = ( C_max - C_tau )  C_avg * 100 ","category":"page"},{"location":"parameters/#Accind:-Accumulation-index","page":"Parameters","title":"Accind: Accumulation index","text":"","category":"section"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"Accind = frac11 - exp(-lambda_z tau)","category":"page"},{"location":"parameters/#MRTtauinf","page":"Parameters","title":"MRTtauinf","text":"","category":"section"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"MRT_tauinf = (AUMC_tau + tau * (AUC_infty - AUC_tau))  AUC_tau","category":"page"},{"location":"parameters/#Swing","page":"Parameters","title":"Swing","text":"","category":"section"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"Swing = (C_max - C_tau min)  C_tau min","category":"page"},{"location":"parameters/#Swingtau","page":"Parameters","title":"Swingtau","text":"","category":"section"},{"location":"parameters/","page":"Parameters","title":"Parameters","text":"Swing_tau = (C_max - C_tau)  C_tau","category":"page"}]
}
