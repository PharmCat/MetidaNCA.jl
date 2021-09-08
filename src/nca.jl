# Pharmacokinetics
# Makoid C, Vuchetich J, Banakar V. 1996-1999. Basic Pharmacokinetics.

function firstobs(time::Vector{T}, obs::Vector, dosetime) where T <: Number
    @inbounds for i = 1:length(time)
        if time[i] >= dosetime && !isnanormissing(obs[i]) return i end
    end
    error("Observations not found")
end
function firstobs(time::Vector{T}, obs, vol, dosetime) where T <: Tuple
    @inbounds for i = 1:length(time)
        if time[i][1] >= dosetime && !isnanormissing(obs[i]) && !isnanormissing(vol[i]) return i end
    end
    error("Observations not found")
end

function ctaumin(time::AbstractVector, obs::AbstractVector, taulastp::Int)
    min = first(obs)
    if length(obs) == 1 return min end
    @inbounds for i = 2:taulastp
        if  obs[i] < min  min = obs[i] end
    end
    min
end

function ctmax(time::AbstractVector, obs::AbstractVector{T}, taulastp) where T
    cmax  = obs[1]
    tmaxn = 1
    if length(obs) == 1 return cmax, first(time), tmaxn end
    taulastp > length(obs) && error("taulastp > length(obs")
    @inbounds for i = 2:taulastp
        if obs[i] > cmax
            cmax  = obs[i]
            tmaxn = i
        end
    end
    return cmax, time[tmaxn], tmaxn
end
function ctmax(time::AbstractVector, obs::AbstractVector{T}) where T
    f = findfirst(!isnanormissing, obs)
    cmax  = obs[f]
    tmaxn = f
    if length(obs) - f == 0 return cmax, time[f], tmaxn end
    @inbounds for i = f+1:length(obs)
        if !isnanormissing(obs[i]) && obs[i] > cmax
            cmax  = obs[i]
            tmaxn = i
        end
    end
    return cmax, time[tmaxn], tmaxn
end
function ctmax(data::PKSubject)
    fobs = firstobs(data.time, data.obs, data.dosetime.time)
    if  data.dosetime.tau > 0
        taulastp = findlast(x -> x <= data.dosetime.time + data.dosetime.tau, data.time)
    else
        taulastp = length(data.obs)
    end
    cmax  = data.obs[fobs]
    tmaxn = fobs
    if length(data.obs) - fobs == 0 return cmax, data.time[fobs], tmaxn end
    @inbounds for i = fobs + 1:length(data.obs)
        if !isnanormissing(data.obs[i]) && data.obs[i] > cmax
            cmax  = data.obs[i]
            tmaxn = i
        end
    end
    return cmax, data.time[tmaxn], tmaxn
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

function slope(x, y)
    if length(x) != length(y) throw(ArgumentError("Unequal vector length!")) end
    n   = length(x)
    if n < 2 throw(ArgumentError("n < 2!")) end
    Σxy::Float64 = zero(Float64)
    Σx::Float64  = zero(Float64)
    Σy::Float64  = zero(Float64)
    Σx2::Float64 = zero(Float64)
    Σy2::Float64 = zero(Float64)
    @inbounds for i = 1:n
        xi = x[i]
        yi = y[i]
        Σxy += xi  * yi
        Σx  += xi
        Σy  += yi
        Σx2 += xi^2
        Σy2 += yi^2
    end
    a   = (n * Σxy - Σx * Σy)/(n * Σx2 - Σx^2)
    b   = (Σy * Σx2 - Σx * Σxy)/(n * Σx2 - Σx^2)
    r2  = (n * Σxy - Σx * Σy)^2/((n * Σx2 - Σx^2)*(n * Σy2 - Σy^2))

    n > 2 ? ar  = 1 - (1 - r2)*(n - 1)/(n - 2) : ar = NaN

    return a, b, r2, ar
end
function logslope(x, y)
    if length(x) != length(y) throw(ArgumentError("Unequal vector length!")) end
    n   = length(x)
    if n < 2 throw(ArgumentError("n < 2!")) end
    Σxy::Float64 = zero(Float64)
    Σx::Float64  = zero(Float64)
    Σy::Float64  = zero(Float64)
    Σx2::Float64 = zero(Float64)
    Σy2::Float64 = zero(Float64)
    @inbounds for i = 1:n
        xi = x[i]
        yi = log(y[i])
        Σxy += xi * yi
        Σx  += xi
        Σy  += yi
        Σx2 += xi^2
        Σy2 += yi^2
    end
    a   = (n * Σxy - Σx * Σy)/(n * Σx2 - Σx^2)
    b   = (Σy * Σx2 - Σx * Σxy)/(n * Σx2 - Σx^2)
    r2  = (n * Σxy - Σx * Σy)^2/((n * Σx2 - Σx^2)*(n * Σy2 - Σy^2))

    n > 2 ? ar  = 1 - (1 - r2)*(n - 1)/(n - 2) : ar = NaN

    return a, b, r2, ar
end #end slope

#---------------------------------------------------------------------------
function aucpart(t₁, t₂, c₁::T, c₂::T, calcm, aftertmax) where T
    if calcm == :lint || c₁ <= zero(T) && c₂ <= zero(T)
        auc   =  linauc(t₁, t₂, c₁, c₂)
    elseif calcm == :logt && aftertmax && c₁ > zero(T) && c₂ > zero(T)
        auc   =  logauc(t₁, t₂, c₁, c₂)
    elseif calcm == :luld &&  c₁ > c₂ > zero(T)
        auc   =  logauc(t₁, t₂, c₁, c₂)
    elseif calcm == :luldt && aftertmax && c₁ > c₂ > zero(T)
        auc   =  logauc(t₁, t₂, c₁, c₂)
    #elseif calcm == :log && c₁ > zero(T) && c₂ > zero(T)
        #auc   =  logauc(t₁, t₂, c₁, c₂)
    else
        auc   =  linauc(t₁, t₂, c₁, c₂)
    end
    return auc
end
function aumcpart(t₁, t₂, c₁::T, c₂::T, calcm, aftertmax) where T
    if calcm == :lint || c₁ <= zero(T) && c₂ <= zero(T)
        aumc  = linaumc(t₁, t₂, c₁, c₂)
    elseif calcm == :logt && aftertmax && c₁ > zero(T) && c₂ > zero(T)
        aumc  = logaumc(t₁, t₂, c₁, c₂)
    elseif calcm == :luld &&  c₁ > c₂ > zero(T)
        aumc  = logaumc(t₁, t₂, c₁, c₂)
    elseif calcm == :luldt && aftertmax && c₁ > c₂ > zero(T)
        aumc  = logaumc(t₁, t₂, c₁, c₂)
    #elseif calcm == :log && c₁ > zero(T) && c₂ > zero(T)
        #aumc  = logaumc(t₁, t₂, c₁, c₂)
    else
        aumc  = linaumc(t₁, t₂, c₁, c₂)
    end
    return aumc
end
#---------------------------------------------------------------------------
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

################################################################################
# STEPS
################################################################################
# 1
function step_1_filterpksubj!(time, obs, dosetime)
    fobs     = firstobs(time, obs, dosetime)
    if fobs > 1
        inds = collect(1:fobs-1)
    else
        inds = Vector{Int}(undef, 0)
    end
    @inbounds for i = fobs:length(obs)
        if isnanormissing(obs[i])
            push!(inds, i)
        end
    end
    if length(inds) > 0
        deleteat!(time, inds)
        deleteat!(obs, inds)
    end
    time, obs
end
#=
function step_1_filterpksubj(time, obs, dosetime)
    fobs     = firstobs(time, obs, dosetime)
    ni = 0
    @inbounds for i = fobs:length(obs)
        isnanormissing(obs[i]) || begin ni += 1 end
    end
    inds = Vector{Int}(undef, ni)
    ni = 1
    @inbounds for i = fobs:length(obs)
        isnanormissing(obs[i]) || begin
        inds[ni] = i
        ni += 1
        end
    end
    time_cp = time[inds]
    obs_cp  = obs[inds]
    time_cp, obs_cp
end
=#
# 3
function step_3_elim!(result, data::PKSubject{T,O}, adm, tmaxn, time_cp, obs_cp, time) where T where O
    obsnum = length(time_cp)
    excltime = time[data.kelrange.kelexcl]
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
                exclinds = findall(x -> x in excltime, time_cp)
                filter!(x-> x ∉ exclinds, timep)
            end
            zcinds = findall(x -> x <= 0, obs_cp)
            filter!(x-> x ∉ zcinds, timep)
            if length(timep) > 2
                logconc    = log.(obs_cp)
                for i = length(timep)-2:-1:1
                    timepv = view(timep, i:length(timep))
                    sl = slope(view(time_cp, timepv), view(logconc, timepv))
                    if sl[1] < 0
                        push!(keldata, time_cp[timep[i]], time_cp[timep[end]], sl[1], sl[2], sl[3], sl[4])
                    end
                end
            end
        end
    else
        stimep = findfirst(x -> x >= time[data.kelrange.kelstart], time_cp)
        etimep = findlast(x -> x <= time[data.kelrange.kelend], time_cp)
        timep = collect(stimep:etimep)
        if length(data.kelrange.kelexcl) > 0
            @inbounds for i in data.kelrange.kelexcl
                filter!(x-> x ∉ findall(x -> x in excltime, time_cp), timep)
            end
        end
        zcinds = findall(x -> x <= 0, obs_cp)
        filter!(x-> x ∉ zcinds, timep)
        if length(timep) > 1
            sl = logslope(view(time_cp, timep), view(obs_cp, timep))
            push!(keldata, time_cp[stimep], time_cp[etimep], sl[1], sl[2], sl[3], sl[4])
        end
    end
    keldata, excltime
end
# 6
function step_6_areas(time_cp, obs_cp, calcm, tmaxn, tlastn)
    obsnum = length(time_cp)
    aucpartl  = Array{Float64, 1}(undef, obsnum - 1)
    aumcpartl = Array{Float64, 1}(undef, obsnum - 1)
    #Calculate all AUC/AUMC part based on data
    for i = 1:(obsnum - 1)
        aucpartl[i]  = aucpart(time_cp[i], time_cp[i + 1], obs_cp[i], obs_cp[i + 1], calcm, i >= tmaxn)
        aumcpartl[i] = aumcpart(time_cp[i], time_cp[i + 1], obs_cp[i], obs_cp[i + 1], calcm, i >= tmaxn)
    end

    #-----------------------------------------------------------------------
    #-----------------------------------------------------------------------
    auclast  = 0
    aumclast = 0
    @inbounds for i = 1:tlastn-1
        auclast  += aucpartl[i]
        aumclast += aumcpartl[i]
    end
    aucall  = auclast
    if tlastn < length(time_cp)
        @inbounds for i = tlastn:obsnum-1
            aucall  += aucpartl[i]
        end
    end
    aucpartl, aumcpartl, auclast, aumclast, aucall
end
"""
    nca(args...; kelauto = true,  elimrange = ElimRange(), dosetime = DoseTime(), kwargs...)

    nca(data, time, conc, sort; kelauto = true,  elimrange = ElimRange(), dosetime = DoseTime(), kwargs...)

    nca(data, time, conc; kelauto = true,  elimrange = ElimRange(), dosetime = DoseTime(), kwargs...)

    nca(time, conc; kelauto = true,  elimrange = ElimRange(), dosetime = DoseTime(), kwargs...)

Import data and perform NCA analysis.

Syntax simillar to [`pkimport`](@ref)

Applicable `kwargs` see  [`nca!`](@ref).

"""
function nca(args...; type = :bps, kelauto = true,  elimrange = ElimRange(), dosetime = DoseTime(), kwargs...)
    pki    = pkimport(args...; kelauto = kelauto,  elimrange = elimrange, dosetime = dosetime)
    #kwargs = Dict{Symbol, Any}(kwargs)
    nca!(pki; kwargs...)
end
"""
    nca!(data::PKSubject{T,O}; adm = :ev, calcm = :lint, intpm = nothing, limitrule = nothing, verbose = 0, warn = true, io::IO = stdout, modify! = nothing) where T where O

* `adm` - administration:
    - `:ev` - extra vascular;
    - `:iv` - intravascular bolus;
* `calcm` - AUC/AUMC calculation method:
    - `:lint` - linear trapezoidal;
    - `:logt` - log-trapezoidal after Tmax;
    - `:luld` - linar up log down;
    - `:luldt` - linear up log down after Tmax;
* `intpm` - interpolation method:
    - `:lint` - linear trapezoidal;
    - `:logt` - log-trapezoidal after Tmax;
    - `:luld` - linar up log down;
    - `:luldt` - linear up log down after Tmax;
* `limitrule` use limitrule for data;
* `verbose` - print to `io`, 1: partial areas table, 2: 1, and results;
* `warn` - show warnings;
* `io` - output stream;
* `modify!` - function to modify output paramaters, call `modify!(data, result)` if difined.

"""
function nca!(data::PKSubject{T,O}; adm = :ev, calcm = :lint, intpm = nothing, limitrule::LimitRule = LimitRule(), verbose = 0, warn = true, io::IO = stdout, modify! = identity) where T where O

    result   = Dict{Symbol, Float64}()

    options =  Dict(:type => :bps, :adm => adm, :calcm => calcm, :intpm => intpm, :limitrule => limitrule, :verbose => verbose, :warn => warn, :modify! => modify!)

    if verbose > 0
        println(io, "  Non-compartmental Pharmacokinetic Analysis")
        if length(data.id) > 0
            print(io, "    Subject: ")
            for (k,v) in data.id
                print(io, "$(k) => $(v); ")
            end
            println(io, "")
        end
        println(io, "    Settings:")
        println(io, "    Method: $(calcm); Dose: $(data.dosetime.dose); Dose time: $(data.dosetime.time)")
        if data.dosetime.tau > 0
            println(io, "    Tau: $(data.dosetime.tau)")
        end
    end
    if isnothing(intpm) intpm = calcm end

    auctype  = promote_type(eltype(data.time), eltype(data.obs))

    if isapplicable(limitrule)
        time, obs = applylimitrule!(deepcopy(data.time), deepcopy(data.obs), limitrule)
    else
        time             = deepcopy(data.time)
        obs              = deepcopy(data.obs)
    end
################################################################################
    # STEP 1 FILTER ALL BEFORE DOSETIME AND ALL NAN OR MISSING VALUES
    time_cp, obs_cp = step_1_filterpksubj!(time, obs, data.dosetime.time)
    if length(obs_cp) < 2
        return NCAResult(data, options, result)
    end
################################################################################
    # STEP 2 - CMAX TMAX FOR TAU RANGE Clast Tlast
    result[:Obsnum] = obsnum = length(obs_cp)
    # If TAU set, calculates start and end timepoints for AUCtau
    if  data.dosetime.tau > zero(typeof(data.dosetime.tau))
        taulastp = findlast(x -> x <= data.dosetime.time + data.dosetime.tau, time_cp)
        result[:Ctaumin] = ctaumin(time_cp, obs_cp, taulastp)
    else
        taulastp = length(obs_cp)
    end
    result[:Cmax], result[:Tmax], tmaxn = ctmax(time_cp, obs_cp, taulastp)
    # C last and T last
    tlastn = 0
    @inbounds for i = length(obs_cp):-1:1
        if obs_cp[i] > zero(O)
            result[:Tlast]   = time_cp[i]
            result[:Clast]   = obs_cp[i]
            tlastn           = i
            break
        end
    end
################################################################################
    # STEP 3
    # Elimination
    keldata, excltime = step_3_elim!(result, data, adm, tmaxn, time_cp, obs_cp, data.time)
################################################################################
    # STEP 4
    if  data.dosetime.time > 0
        time_cp .-= data.dosetime.time
    end
################################################################################
    # STEP 5
    # Dose concentration
    # Dosetime is first point
    #local doseaucpart = zero(Float64)
    #local doseaumcpart = zero(Float64)
    cdoseins = zero(Int)
    if  iszero(first(time_cp))
        result[:Cdose] = first(obs_cp)
    # Dosetime before first point
    else
        if adm == :iv
            if  first(obs_cp) > obs_cp[2] > zero(O)
                result[:Cdose] = logcpredict(first(time_cp), time_cp[2], 0, first(obs_cp), obs_cp[2])
            else
                result[:Cdose] = first(obs_cp)
            end
        else
            if  data.dosetime.tau > zero(typeof(data.dosetime.tau))
                result[:Cdose] = result[:Ctaumin]
            else
                result[:Cdose] = zero(O)
            end
        end
        cdoseins = 1
        pushfirst!(time_cp, zero(T))
        pushfirst!(obs_cp, result[:Cdose])
        # if time-zero point added some points shoud be shifted
        taulastp += 1
        tmaxn += 1
        tlastn += 1
    end
################################################################################
    # STEP 6
    # Areas
    aucpartl, aumcpartl, auclast, aumclast, aucall = step_6_areas(time_cp, obs_cp, calcm, tmaxn, tlastn)
    result[:AUClast]   = auclast
    result[:AUMClast]  = aumclast
    result[:AUCall]    = aucall
################################################################################
    # STEP 7
    # Other parameters
    #---------------------------------------------------------------------------
    result[:MRTlast]    = result[:AUMClast] / result[:AUClast]
    #---------------------------------------------------------------------------
    if data.dosetime.dose > 0
        result[:Cllast]           = data.dosetime.dose / result[:AUClast]
        result[:Dose]             = data.dosetime.dose
    end
    #-----------------------------------------------------------------------
    #-----------------------------------------------------------------------
    tlagn = findfirst(!iszero, obs_cp)
    if tlagn > 1 result[:Tlag] = time_cp[tlagn-1] else result[:Tlag] = zero(Float64) end

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
        result[:AUCpct]          = (result[:AUCinf] - result[:AUClast]) / result[:AUCinf] * 100
        result[:AUMCinf]         = result[:AUMClast] + result[:Tlast] * result[:Clast] / result[:Kel] + result[:Clast] / result[:Kel] ^ 2
        result[:MRTinf]          = result[:AUMCinf] / result[:AUCinf]
        if data.dosetime.dose > 0
            result[:Vzlast]          = data.dosetime.dose / result[:AUClast] / result[:Kel]
            result[:Vzinf]           = data.dosetime.dose / result[:AUCinf] / result[:Kel]
            result[:Clinf]           = data.dosetime.dose / result[:AUCinf]
            result[:Vssinf]          = result[:Clinf] * result[:MRTinf]
        end
    else
        result[:Kel] = NaN
    end
################################################################################
    # STEP 8
    # Steady-state parameters
    if data.dosetime.tau > 0
        eaucpartl  = eaumcpartl = 0.0
        if time_cp[taulastp] < data.dosetime.tau < time_cp[end]
            result[:Ctau] = interpolate(time_cp[taulastp], time_cp[taulastp + 1], data.dosetime.tau, obs_cp[taulastp], obs_cp[taulastp + 1], intpm, true)
            eaucpartl = aucpart(time_cp[taulastp], data.dosetime.tau, obs_cp[taulastp], result[:Ctau], calcm, true)
            eaumcpartl = aumcpart(time_cp[taulastp], data.dosetime.tau, obs_cp[taulastp], result[:Ctau], calcm, true)
                #remoove part after tau
        elseif data.dosetime.tau > time_cp[end] && result[:Kel] !== NaN
                #extrapolation
            result[:Ctau] = exp(result[:LZint] + result[:LZ] * (data.dosetime.tau + data.dosetime.time))
            eaucpartl = aucpart(time_cp[end], data.dosetime.tau, obs_cp[end], result[:Ctau], calcm, true)
            eaumcpartl = aumcpart(time_cp[end], data.dosetime.tau, obs_cp[end], result[:Ctau], calcm, true)
        else
            result[:Ctau] = obs_cp[taulastp]
        end

        auctau   = eaucpartl
        aumctau  = eaumcpartl
        @inbounds for i = 1:taulastp-1
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
        result[:Fluc]     = (result[:Cmax] - result[:Ctaumin])/result[:Cavg] * 100
        result[:Fluctau]  = (result[:Cmax] - result[:Ctau])/result[:Cavg] * 100
        #If Kel calculated
        if result[:Kel] !== NaN
            result[:Accind]   = 1 / (1 - (exp(-result[:Kel] * data.dosetime.tau)))
        end
        result[:MRTtauinf]       = (result[:AUMCtau] + data.dosetime.tau * (result[:AUCinf] - result[:AUCtau])) / result[:AUCtau]
        result[:Vztau]           = data.dosetime.dose / result[:AUCtau] / result[:Kel]
        result[:Cltau]           = data.dosetime.dose / result[:AUCtau]
        result[:Vsstau]          = result[:Cltau] * result[:MRTtauinf]
    end
################################################################################
    # Verbose output
    if verbose > 0
        aucpartlsum  = similar(aucpartl)
        aumcpartlsum = similar(aumcpartl)
        @inbounds for i = 1:length(aucpartl)
            aucpartlsum[i]  = sum(view(aucpartl, 1:i))
            aumcpartlsum[i] = sum(view(aumcpartl, 1:i))
        end
        #aucpartlsum  .+= doseaucpart
        #aumcpartlsum .+= doseaumcpart
        if  data.dosetime.time > 0
            time_cp .+= data.dosetime.time
        end
        hnames = [:Time, :Concentrtion, :AUC, :AUC_cum, :AUMC, :AUMC_cum, :Info]
        mx = metida_table(collect(time_cp), collect(obs_cp), pushfirst!(aucpartl, 0.0),  pushfirst!(aucpartlsum, 0.0), pushfirst!(aumcpartl, 0.0),  pushfirst!(aumcpartlsum, 0.0), fill("", length(obs_cp));
        names = hnames)
        if cdoseins > 0
            #println(io, "    Dose interpolated part AUC $(doseaucpart); AUMC $(doseaumcpart)")
            #mx[1,3] = mx[1,4] = doseaucpart
            #mx[1,5] = mx[1,6] = doseaumcpart
            #pushfirst!(mx, [data.dosetime.time, result[:Cdose], 0.0, 0.0, 0.0, 0.0,"D*"])
            #ins = 1
            mx[1, 7] = "D*"
        else
            mx[1, 7] = "D"
        end
        if !isnan(result[:Kel])
            @inbounds for i = 1:length(time_cp)
                if time_cp[i] >= keldata.s[rsqn] && time_cp[i] <= keldata.e[rsqn]
                    if length(data.kelrange.kelexcl) > 0
                        if time_cp[i] in excltime
                            mx[i, 7] = "Excl"
                        else
                            mx[i, 7] = "E"
                        end
                    else
                        mx[i, 7] = "E"
                    end
                end
            end
        end
        hnames = (["Time" "Concentrtion" "AUC"  "AUC" "AUMC" "AUMC" "Info"],
                  ["" "" "" "(cumulate)" "" "(cumulate)" ""])
        PrettyTables.pretty_table(io, mx; tf = PrettyTables.tf_compact, header = hnames)
        println(io, "")
        println(io, "    Cdose: $(result[:Cdose]), Dose time: $(data.dosetime.time)")
        println(io, "    Kel start: $(keldata.s[rsqn]); end: $(keldata.e[rsqn])")
        println(io, "")
        if data.dosetime.tau < time_cp[end] && data.dosetime.tau > 0
            println(io, "    Tau + dosetime is less then end time. Interpolation used.")
            println(io, "    Ctau: $(result[:Ctau])")
            println(io, "    AUC  final part: $(eaucpartl)")
            println(io, "    AUMC final part: $(eaumcpartl)")
            println(io, "")
        end
        if verbose > 1
            println(io, "    Results:")
            PrettyTables.pretty_table(io, result; tf = PrettyTables.tf_compact)
        end
    end
################################################################################

    ncares = NCAResult(data, options, result)
    modify!(ncares)

    #-----------------------------------------------------------------------
    return ncares
end

"""
    nca!(data::DataSet{T1}; adm = :ev, calcm = :lint, intpm = nothing, verbose = 0, warn = true, io::IO = stdout) where T1 <: PKSubject{T,O}  where T  where O

Non-compartmental (NCA) analysis of pharmacokinetic (PK) data.
"""
function nca!(data::DataSet{Subj}; adm = :ev, calcm = :lint, intpm = nothing, limitrule::LimitRule = LimitRule(), verbose = 0, warn = true, io::IO = stdout, modify! = identity) where Subj <: AbstractSubject
    result = Vector{NCAResult{Subj}}(undef, length(data))
    for i = 1:length(data)
        result[i] = nca!(data[i]; adm = adm, calcm = calcm, intpm = intpm, limitrule = limitrule, verbose = verbose, warn = warn, io = io, modify! = modify!)
    end
    DataSet(result)
end


function maxconc(subj::T) where T <: PKSubject
    maximum(subj.obs)
end
function minconc(subj::T) where T <: PKSubject
    minimum(subj.obs)
end

function exrate(time, conc, vol)
    length(time) == length(conc) == length(vol) || error("")
    er = Vector{Float64}(undef, length(time))
    @inbounds for i = 1:length(time)
        er[i] = conc[i]*vol[i]/(time[i][2] - time[i][1])
    end
    er
end

function step_1_filterupksubj(time, obs, vol, dosetime)
    fobs     = firstobs(time, obs, vol, dosetime)
    ni = 0
    @inbounds for i = fobs:length(obs)
        if !isnanormissing(obs[i]) && !isnanormissing(vol[i])
            ni += 1
        end
    end
    inds = Vector{Int}(undef, ni)
    ni = 1
    @inbounds for i = fobs:length(obs)
        if !isnanormissing(obs[i]) && !isnanormissing(vol[i])
            inds[ni] = i
            ni += 1
        end
    end
    time_cp = time[inds]
    obs_cp  = obs[inds]
    vol_cp  = vol[inds]
    time_cp, obs_cp, vol_cp
end

function nca!(data::UPKSubject{T, O, VOL, V}; adm = :ev, calcm = :lint, intpm = nothing, limitrule::LimitRule = LimitRule(), verbose = 0, warn = true, io::IO = stdout, modify! = identity) where T where O where VOL where V

    result   = Dict{Symbol, Float64}()

    options =  Dict(:type => :urine, :adm => adm, :calcm => calcm, :intpm => intpm, :limitrule => limitrule, :verbose => verbose, :warn => warn, :modify! => modify!)

    if verbose > 0
        println(io, "  Non-compartmental Pharmacokinetic Analysis")
        println(io, "  Matrix: urine")
        if length(data.id) > 0
            print(io, "    Subject: ")
            for (k,v) in data.id
                print(io, "$(k) => $(v); ")
            end
            println(io, "")
        end
        println(io, "    Settings:")
        println(io, "    Method: $(calcm); Dose: $(data.dosetime.dose); Dose time: $(data.dosetime.time)")
    end
    if isnothing(intpm) intpm = calcm end

    #=
    if isapplicable(limitrule)
        mtime, obs = applylimitrule!(mtime, deepcopy(data.obs), limitrule)
    else
        mtime  = mtime
        obs    = data.obs
    end
    =#

    time, obs, vol = step_1_filterupksubj(data.time, data.obs, data.vol, data.dosetime.time)

    mtime  = map(x-> (x[1]+x[2])/2, time)

    exr  = exrate(time, obs, vol)

    result[:AR]      = data.obs' * data.vol
    result[:Vol]  = sum(vol)

    if time[1][1] > data.dosetime.time
        pushfirst!(mtime, 0)
    else
        pushfirst!(mtime, time[1][1])
    end
    pushfirst!(obs, 0)
    pushfirst!(vol, 0)
    pushfirst!(exr, 0)

    result[:Maxrate], result[:Tmax], tmaxn = ctmax(mtime, exr)

    if data.dosetime.dose > 0
        100*Amount_Recovered/Dose
        result[:prec]  = result[:AR]/data.dosetime.dose * 100
    end
    obsnum = length(exr)

    #result[:Kel]
    #result[:HL]
    aucpartl  = Array{Float64, 1}(undef, obsnum - 1)
    #Calculate all AUC part based on data
    for i = 1:(obsnum - 1)
        aucpartl[i] = aucpart(mtime[i], mtime[i + 1], exr[i], exr[i + 1], calcm, i >= tmaxn)
    end

    result[:AUCall]  = sum(aucpartl)

    ncares = NCAResult(data, options, result)
    modify!(ncares)
    #-----------------------------------------------------------------------
    return ncares
end
