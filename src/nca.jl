# Pharmacokinetics
# Makoid C, Vuchetich J, Banakar V. 1996-1999. Basic Pharmacokinetics.


isnanormissing(x) = isnan(x) || ismissing(x)

function findnextnnom(v, range)
    for i in range
        if !isnanormissing(v[i]) return i end
    end
end

function firstobs(time::Vector, obs::Vector, dosetime)
    @inbounds for i = 1:length(time)
        if time[i] >= dosetime || !isnanormissing(obs[i]) return i end
    end
end

function ctaumin(time::Vector{T}, obs::Vector, fobs, tautime) where T
    ttauminp = 0
    @inbounds for i = length(obs):-1:fobs
        if  time[i] < tautime && !isnanormissing(obs[i])
            ttauminp = i
            break
        end
    end
    min = obs[fobs]
    if ttauminp > fobs
        @inbounds for i = fobs+1:ttauminp
            if !isnanormissing(obs[i]) && obs[i] < min min = obs[i] end
        end
    end
    min
end

function firstabovezero(obs::Vector{T}, range) where T
    @inbounds for i in range
        if obs[i] > zero(T) return obs[i] end
    end
end

function firstnnom(obs::Vector{T}, range) where T
    @inbounds for i in range
        if !isnanormissing(obs[i]) return obs[i] end
    end
end

function ctmax(time::AbstractVector, obs::AbstractVector{T}) where T
    cmax  = obs[1]
    tmaxn = 1
    @inbounds for i = 2:length(obs)
        if obs[i] > cmax
            cmax  = obs[i]
            tmaxn = i
        end
    end
    return cmax, time[tmaxn], tmaxn
end

function logcpredict(t₁, t₂, tx, c₁, c₂)
    return exp(log(c₁) + (tx-t₁)/(t₂-t₁)*(log(c₂) - log(c₁)))
end

#Linear trapezoidal auc
function linauc(t₁, t₂, c₁, c₂)
    return (t₂-t₁)*(c₁+c₂)/2
end
#Linear trapezoidal aumc
function linaumc(t₁, t₂, c₁, c₂)
    return (t₂-t₁)*(t₁*c₁+t₂*c₂)/2
end
#Log trapezoidal auc
function logauc(t₁, t₂, c₁, c₂)
    return  (t₂-t₁)*(c₂-c₁)/log(c₂/c₁)
end
#Log trapezoidal aumc
function logaumc(t₁, t₂, c₁, c₂)
    return (t₂-t₁) * (t₂*c₂-t₁*c₁) / log(c₂/c₁) - (t₂-t₁)^2 * (c₂-c₁) / log(c₂/c₁)^2
end
#Intrapolation
#linear prediction bx from ax, a1 < ax < a2
function linpredict(a₁, a₂, ax, b₁, b₂)
    return (ax - a₁) / (a₂ - a₁)*(b₂ - b₁) + b₁
end


function cpredict(t₁, t₂, tx, c₁, c₂, calcm)
    if calcm == :lint || c₂ >= c₁
        return linpredict(t₁, t₂, tx, c₁, c₂)
    else
        return logcpredict(t₁, t₂, tx, c₁, c₂)
    end
end

function slope(x, y)
    if length(x) != length(y) throw(ArgumentError("Unequal vector length!")) end
    n   = length(x)
    if n < 2 throw(ArgumentError("n < 2!")) end
    Σxy = 0.0
    Σx  = 0.0
    Σy  = 0.0
    Σx2 = 0.0
    Σy2 = 0.0
    @inbounds for i = 1:n
        Σxy += x[i] * y[i]
        Σx  += x[i]
        Σy  += y[i]
        Σx2 += x[i]^2
        Σy2 += y[i]^2
    end
    a   = (n * Σxy - Σx * Σy)/(n * Σx2 - Σx^2)
    b   = (Σy * Σx2 - Σx * Σxy)/(n * Σx2 - Σx^2)
    r2  = (n * Σxy - Σx * Σy)^2/((n * Σx2 - Σx^2)*(n * Σy2 - Σy^2))
    if n > 2
        ar  = 1 - (1 - r2)*(n - 1)/(n - 2)
    else
        ar = NaN
    end
    return a, b, r2, ar
end #end slope

#---------------------------------------------------------------------------
#=
function aucpart(t₁, t₂, c₁, c₂, calcm)
    if calcm == :lint
        auc   =  linauc(t₁, t₂, c₁, c₂)
        aumc  = linaumc(t₁, t₂, c₁, c₂)
    else
        if c₁ > 0 && c₂ > 0
            auc   =  logauc(t₁, t₂, c₁, c₂)
            aumc  = logaumc(t₁, t₂, c₁, c₂)
        else
            auc   =  linauc(t₁, t₂, c₁, c₂)
            aumc  = linaumc(t₁, t₂, c₁, c₂)
        end
    end
    return auc, aumc
end
=#
function aucpart(t₁, t₂, c₁::T, c₂::T, calcm, aftertmax) where T
    if calcm == :lint || c₁ <= zero(T) && c₂ <= zero(T)
        auc   =  linauc(t₁, t₂, c₁, c₂)
        aumc  = linaumc(t₁, t₂, c₁, c₂)
    elseif calcm == :logt && aftertmax && c₁ > zero(T) && c₂ > zero(T)
        auc   =  logauc(t₁, t₂, c₁, c₂)
        aumc  = logaumc(t₁, t₂, c₁, c₂)
    elseif calcm == :luld &&  c₁ > c₂ > zero(T)
        auc   =  logauc(t₁, t₂, c₁, c₂)
        aumc  = logaumc(t₁, t₂, c₁, c₂)
    elseif calcm == :luldt && aftertmax && c₁ > c₂ > zero(T)
        auc   =  logauc(t₁, t₂, c₁, c₂)
        aumc  = logaumc(t₁, t₂, c₁, c₂)
    #elseif calcm == :log && c₁ > zero(T) && c₂ > zero(T)
        #auc   =  logauc(t₁, t₂, c₁, c₂)
        #aumc  = logaumc(t₁, t₂, c₁, c₂)
    else
        auc   =  linauc(t₁, t₂, c₁, c₂)
        aumc  = linaumc(t₁, t₂, c₁, c₂)
    end
    return auc, aumc
end

function interpolate(t₁, t₂, tx, c₁::T, c₂::T, intpm, aftertmax) where T
    if intpm == :lint || c₁ <= zero(T) || c₂ <= zero(T)
        c = linpredict(t₁, t₂, tx, c₁, c₂)
    elseif intpm == :logt && aftertmax && c₁ > zero(T) && c₂ > zero(T)
        c = logcpredict(t₁, t₂, tx, c₁, c₂)
    elseif intpm == :luld && c₁ > c₂ > zero(T)
        c = logcpredict(t₁, t₂, tx, c₁, c₂)
    elseif intpm == :luldt && aftertmax && c₁ > c₂ > zero(T)
        c = logcpredict(t₁, t₂, tx, c₁, c₂)
    #elseif intpm == :log && c₁ > zero(T) && c₂ > zero(T)
        #c = logcpredict(t₁, t₂, tx, c₁, c₂)
    else
        c = linpredict(t₁, t₂, tx, c₁, c₂)
    end
    return c
end

function nca!(data::PKSubject{T,O}; adm = :ev, calcm = :lint, intpm = nothing, verbose = false, warn = true, io::IO = stdout) where T where O

    result   = Dict{Symbol, Union{eltype(data.time), eltype(data.obs)}}()

    if verbose
        println(io, "Non-compartmental Pharmacokinetic Analysis")
        printsortval(io, data.sort)
        println(io, "Settings:")
        println(io, "Method: $(calcm)")
    end
    if isnothing(intpm) intpm = calcm end

    time             = data.time
    obs              = data.obs
    auctype  = promote_type(eltype(time), eltype(obs))
    fobs             = firstobs(data.time, data.obs, data.dosetime.time)


    if length(obs) - fobs < 2
        return NCAResult(data, calcm, result, data.id)
    end


    # If TAU set, calculates start and edt timepoints for AUCtau
    if  data.dosetime.tau > zero(typeof(data.dosetime.tau))
        #taulast, taulastp, result[:Ctaumin] = taulastmin(time, obs, fobs, lobs, tautime)
        tautime = data.dosetime.time + data.dosetime.tau
        result[:Ctaumin] = ctaumin(time, obs, fobs, tautime)
    end

    # Dose concentration
    # Dosetime is first point
    cdoseins = false
    if  data.dosetime.time == time[fobs]
        result[:Cdose] = obs[fobs]
        doseaucpart =  doseaumcpart = 0.0
    # Dosetime before first point
    elseif data.dosetime.time < time[fobs]
        if adm == :iv
            if !isnanormissig(obs[fobs]) && !isnanormissig(obs[fobs + 1]) && obs[fobs] > obs[fobs + 1] > zero(O)
                result[:Cdose] = logcpredict(time[fobs], time[fobs + 1], data.dosetime.time, obs[fobs], obs[fobs+1])
            else
                result[:Cdose] = firstnnom(obs, fobs:lobs)
            end
            doseaucpart, doseaumcpart  = aucpart(data.dosetime.time, time[fobs], result[:Cdose], obs[fobs], calcm, true)
        else
            if  data.dosetime.tau > zero(typeof(data.dosetime.tau))
                result[:Cdose] = result[:Ctaumin]
            else
                result[:Cdose] = zero(O)
            end
            doseaucpart, doseaumcpart  = aucpart(data.dosetime.time, time[fobs], result[:Cdose], obs[fobs], calcm, false)
        end
        cdodeins = true
    else
        error("Dose time after first observation!") # Remove
    end

    aucinds = filter!(x-> x ∉ findall(isnanormissing, obs), collect(fobs:length(obs)))
    #println(fobs)
    time_auc = view(time, aucinds)
    obs_auc  = view(obs, aucinds)
    #println(length(time_auc))
    result[:Obsnum] = obsnum = length(obs_auc)
    #CMAX TMAX FOR TAU RANGE
    result[:Cmax], result[:Tmax], tmaxn = ctmax(time_auc, obs_auc)

    #Areas
    aucpartl  = Array{auctype, 1}(undef, obsnum - 1)
    aumcpartl = Array{auctype, 1}(undef, obsnum - 1)
    #Calculate all UAC part based on data
    for i = 1:(obsnum - 1)
        aucpartl[i], aumcpartl[i] = aucpart(time_auc[i], time_auc[i + 1], obs_auc[i], obs_auc[i + 1], calcm, i >= tmaxn)
    end

    # C last and T last
    tlastn = 0
    for i = obsnum:-1:1
        if obs_auc[i] > zero(O)
            result[:Tlast]   = time_auc[i]
            result[:Clast]   = obs_auc[i]
            tlastn           = i
            break
        end
    end

    #-----------------------------------------------------------------------
    #-----------------------------------------------------------------------

    #Elimination
    keldata                = KelData()
    if data.kelauto
        if (adm != :iv && obsnum - tmaxn > 2) || (adm == :iv && obsnum - tmaxn > 1)
            if adm == :iv
                stimep = tmaxn
            else
                stimep = tmaxn + 1
            end
            timep = collect(stimep:obsnum)
            if length(data.kelrange.kelexcl) > 0
                excltime = view(time, data.kelrange.kelexcl)
                filter!(x-> x ∉ findall(x -> x in excltime, time_auc), timep)
            end
            filter!(x-> x ∉ findall(x -> x <= 0, obs_auc), timep)
            if length(timep) > 2
                logconc    = log.(obs_auc)
                for i = length(timep)-2:-1:1
                    timepv = view(timep, i:length(timep))
                    sl = slope(view(time_auc, timepv), view(logconc, timepv))
                    if sl[1] < 0
                        push!(keldata, time_auc[timep[i]], time_auc[timep[end]], sl[1], sl[2], sl[3], sl[4])
                    end
                end
            end
        end
    else
        stimep = findfirst(x -> x >= time[data.kelrange.kelstart], time_auc)
        etimep = findlast(x -> x <= time[data.kelrange.kelend], time_auc)
        timep = collect(stimep:etimep)
        if length(data.kelrange.kelexcl) > 0
            for i in data.kelrange.kelexcl
                excltime = view(time, data.kelrange.kelexcl)
                filter!(x-> x ∉ findall(x -> x in excltime, time_auc), timep)
            end
        end
        if length(timep) > 1
            sl = slope(view(time_auc, timep), log.(view(conc_auc, timep)))
            push!(keldata, time_auc[stimep], time_auc[etimep], sl[1], sl[2], sl[3], sl[4])
        end
    end

    auclast  = doseaucpart
    aumclast = doseaumcpart
    for i = 1:tlastn-1
        auclast  += aucpartl[i]
        aumclast += aumcpartl[i]
    end
    aucall  = auclast
    if tlastn < obsnum
        for i = tlastn:obsnum-1
            aucall  += aucpartl[i]
        end
    end
    #-----------------------------------------------------------------------
    #-----------------------------------------------------------------------
    result[:AUClast]   = auclast
    result[:AUMClast]  = aumclast

    result[:AUCall]    = aucall
    #---------------------------------------------------------------------------
    #---------------------------------------------------------------------------
    result[:MRTlast]    = result[:AUMClast] / result[:AUClast]
    #---------------------------------------------------------------------------
    if data.dosetime.dose > 0
        result[:Cllast]           = data.dosetime.dose / result[:AUClast]
    end
    #-----------------------------------------------------------------------
    #-----------------------------------------------------------------------

    if  length(keldata) > 0
        data.keldata             = keldata
        result[:ARsq], rsqn      = findmax(keldata.ar)
        result[:Rsq]             = keldata.r[rsqn]
        result[:Kel]             = abs(keldata.a[rsqn])
        result[:LZ]              = keldata.a[rsqn]
        result[:LZint]           = keldata.b[rsqn]
        result[:Rsqn]            = rsqn
        result[:Clast_pred]      = exp(result[:LZint] + result[:LZ]*result[:Tlast])
        result[:HL]              = LOG2 / result[:Kel]
        result[:AUCinf]          = result[:AUClast] + result[:Clast] / result[:Kel]
        result[:AUCinf_pred]     = result[:AUClast] + result[:Clast_pred] / result[:Kel]
        result[:AUCpct]          = (result[:AUCinf] - result[:AUClast]) / result[:AUCinf] * 100.0
        result[:AUMCinf]         = result[:AUMClast] + result[:Tlast] * result[:Clast] / result[:Kel] + result[:Clast] / result[:Kel] ^ 2
        result[:MRTinf]          = result[:AUMCinf] / result[:AUCinf]
        if data.dosetime.dose > 0
            result[:Vzlast]          = data.dosetime.dose / result[:AUClast] / result[:Kel]
            result[:Vzinf]           = data.dosetime.dose / result[:AUCinf] / result[:Kel]
            result[:Clinf]           = data.dosetime.dose / result[:AUCinf]
            result[:Vssinf]          = result[:Clinf] * result[:MRTinf]
        end
    end
        #-----------------------------------------------------------------------
        #Steady-state

    if data.dosetime.tau > 0
        eaucpartl  = eaumcpartl = 0.0
        taulastp = findlast(x -> x < tautime, time_auc)
        if tautime < time_auc[end]
            result[:Ctau] = interpolate(time_auc[taulastp], time_auc[taulastp + 1], tautime, obs_auc[taulastp], obs_auc[taulastp + 1], intpm, true)
            eaucpartl, eaumcpartl = aucpart(time_auc[taulastp], tautime, obs_auc[taulastp], result[:Ctau], calcm, true)
                #remoove part after tau
        elseif tautime > time_obs[end] && result[:LZint] !== NaN
                #extrapolation
            result[:Ctau] = exp(result[:LZint] + result[:LZ] * tautime)
            eaucpartl, eaumcpartl = aucpart(time_auc[end], tautime, obs_auc[end], result[:Ctau], calcm, true)
        else
            result[:Ctau] = obs[end]
        end

        auctau   = eaucpartl  + doseaucpart
        aumctau  = eaumcpartl + doseaumcpart
        for i = 1:taulastp-1
            auctau  += aucpartl[i]
            aumctau += aumcpartl[i]
        end
        result[:AUCtau]   = auctau
        result[:AUMCtau]  = aumctau

        result[:Cavg]     = result[:AUCtau]/data.dosetime.tau
        if result[:Ctaumin] != 0
            result[:Swing]    = (result[:Cmax] - result[:Ctaumin])/result[:Ctaumin]
        end
        if result[:Ctau] != 0
            result[:Swingtau] = (result[:Cmax] - result[:Ctau])/result[:Ctau]
        end
        result[:Fluc]     = (result[:Cmax] - result[:Ctaumin])/result[:Cavg]
        result[:Fluctau]  = (result[:Cmax] - result[:Ctau])/result[:Cavg]
        #If Kel calculated
        if result[:Kel] !== NaN
            result[:Accind]   = 1 / (1 - (exp(-result[:Kel] * data.dosetime.tau)))
        end
        result[:MRTtauinf]       = (result[:AUMCtau] + data.dosetime.tau * (result[:AUCinf] - result[:AUCtau])) / result[:AUCtau]
        result[:Vztau]           = data.dosetime.dose / result[:AUCtau] / result[:Kel]
        result[:Cltau]           = data.dosetime.dose / result[:AUCtau]
        result[:Vsstau]          = result[:Cltau] * result[:MRTtauinf]
    end
    ############################################################################
    #=
    if verbose
        aucpartlsum  = similar(aucpartl)
        aumcpartlsum = similar(aumcpartl)
        for i = 1:length(aucpartl)
            aucpartlsum[i]  = sum(view(aucpartl, 1:i))
            aumcpartlsum[i] = sum(view(aumcpartl, 1:i))
        end
        astx    = Vector{String}(undef, length(time))
        astx[1] = ""
        for i = 1:length(pmask)
            if pmask[i] astx[i+1] = "Yes" else astx[i+1] = "No" end
        end
        mx = hcat(time, obs, round.(vcat([0], aucpartl), digits = 3),  round.(vcat([0], aucpartlsum), digits = 3), round.(vcat([0], aumcpartl), digits = 3),  round.(vcat([0], aumcpartlsum), digits = 3), astx)
        mx = vcat(permutedims(["Time", "Concentrtion", "AUC", "AUC (cumulate)", "AUMC", "AUMC (cumulate)", "Include"]), mx)
        printmatrix(io, mx)
        println(io, "")
        println(io, "Cdose: $(result[:Cdose]), Dose time: $(data.dosetime.time)")
        if data.dosetime.time > time[1]
            println("Dose AUC  part: $(doseaucpart)")
            println("Dose AUMC part: $(doseaumcpart)")
        end
        println(io, "")
        if tautime < time[end] && tautime > 0
            println(io, "Tau + dosetime is less then end time. Interpolation used.")
            println(io, "Interpolation between: $(time[ncae]) - $( time[ncae + 1]), method: $(intp)")
            println(io, "Ctau: $(result[:Ctau])")
            println(io, "AUC  final part: $(eaucpartl)")
            println(io, "AUMC final part: $(eaumcpartl)")
        end
    end
    =#
    #-----------------------------------------------------------------------
    return NCAResult(data, calcm, result, data.id)
end


function nca!(data::DataSet{T1}; adm = :ev, calcm = :lint, intpm = nothing, verbose = false, warn = true, io::IO = stdout) where T1 <: PKSubject{T,O}  where T  where O
    result = Vector{NCAResult{T1}}(undef, length(data))
    for i = 1:length(data)
        result[i] = nca!(data[i]; adm = adm, calcm = calcm, intpm = intpm, verbose = verbose, warn = warn, io = io)
    end
    DataSet(result)
end
