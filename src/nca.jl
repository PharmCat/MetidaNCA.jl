# Pharmacokinetics
# Makoid C, Vuchetich J, Banakar V. 1996-1999. Basic Pharmacokinetics.

function validobsn(time::Vector, obs::Vector)
    if length(time) != length(obs) 
        error("Vector length `time` not equal `obs`") 
    end
    n = 0
    @inbounds for i = 1:length(time)
        if !isnanormissing(time[i]) && !isnanormissing(obs[i]) n+=1 end
    end
    n
end

function firstobs(time::Vector{T}, obs::Vector, dosetime) where T <: Number
    @inbounds for i = 1:length(time)
        if time[i] >= dosetime && !isnanormissing(obs[i]) return i end #unitful dosetime?
    end
    error("Observations not found")
end
function firstobs(time::Vector{<:Tuple}, obs, vol, dosetime)
    @inbounds for i = 1:length(time)
        if time[i][1] >= dosetime && !isnanormissing(obs[i]) && !isnanormissing(vol[i]) return i end
    end
    error("Observations not found")
end
#=
function ctauminmax(time::AbstractVector, obs::AbstractVector, taulastp::Int)
    fi = 0
    min = obs[1]
    max = obs[1]
    for i = 1:taulastp
        if !isnanormissing(obs[i])
            fi = i
            min = obs[i]
            max = obs[i]
            break
        end
    end
    if length(obs) == 1 return min, max end
    @inbounds for i = fi:taulastp
        if !isnanormissing(obs[i]) && obs[i] < min  min = obs[i] end
        if !isnanormissing(obs[i]) && obs[i] > max  max = obs[i] end
    end
    min, max
end
=#
function ctmax(time::AbstractVector, obs::AbstractVector, taulastp::Int = length(obs)) 
    taulastp > length(obs) && error("taulastp > length(obs")
    f = findfirst(!isnanormissing, view(obs, 1:taulastp))
    if isnothing(f) return NaN, NaN, 0, NaN end
    cmax  = obs[f]
    cmin  = obs[f]
    tmaxn = f
    if taulastp == f return cmax, first(time), tmaxn, cmin end
    @inbounds for i = f+1:taulastp
        if !isnanormissing(obs[i])
            if obs[i] > cmax
                cmax  = obs[i]
                tmaxn = i
            elseif obs[i] < cmin
                cmin = obs[i]
            end
        end
    end
    return cmax, time[tmaxn], tmaxn, cmin
end
   #=
function ctmax(time::AbstractVector, obs::AbstractVector)
 
    f = findfirst(!isnanormissing, obs)
    cmax  = obs[f]
    tmaxn = f
    if length(obs) - f == 0 return cmax, time[f], tmaxn end
    @inbounds for i = f + 1:length(obs)
        if !isnanormissing(obs[i]) && obs[i] > cmax
            cmax  = obs[i]
            tmaxn = i
        end
    end
    return cmax, time[tmaxn], tmaxn
   
end
 =#
function ctmax(data::PKSubject, dt::Int = 0)
    if dt == 0 dt = length(data.dosetime) end
    datadt = data.dosetime[dt]
    fobs = firstobs(gettime(data), getobs(data), datadt.time)
    if  datadt.tau > 0
        taulastp = findlast(x -> x <= datadt.time + datadt.tau, gettime(data))
    else
        taulastp = length(data)
    end
    cmax  = getobs(data)[fobs]
    tmaxn = fobs
    if length(data) - fobs == 0 return cmax, gettime(data)[fobs], tmaxn end
    @inbounds for i = fobs + 1:length(data)
        if !isnanormissing(getobs(data)[i]) && getobs(data)[i] > cmax
            cmax  = getobs(data)[i]
            tmaxn = i
        end
    end
    return cmax, data.time[tmaxn], tmaxn
end

function logtpredict(c₁, c₂, cx, t₁, t₂)
    return log(cx/c₁)/log(c₂/c₁)*(t₂-t₁)+t₁
end

function logcpredict(t₁, t₂, tx, c₁::T, c₂::T) where T
    return exp(log(c₁/oneunit(c₁)) + (tx-t₁)/(t₂-t₁)*(log(c₂/oneunit(c₂)) - log(c₁/oneunit(c₁))))*oneunit(c₁)
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

function slope(x::AbstractVector{X}, y::AbstractVector{Y}; funk::Function = identity) where X where Y
    XYT = promote_type(X, Y)
    if length(x) != length(y) throw(ArgumentError("Unequal vector length!")) end
    n   = length(x)
    if n < 2 throw(ArgumentError("n < 2!")) end
    Σxy = zero(XYT)
    Σx  = zero(X)
    Σy  = zero(Y)
    Σx2 = zero(X)
    Σy2 = zero(Y)
    @inbounds for i = 1:n
        xi = x[i]
        yi = funk(y[i])
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

    return a, b, r2, ar, n
end
#=
function logslope(x, y)
    if length(x) != length(y) throw(ArgumentError("Unequal vector length!")) end
    n   = length(x)
    if n < 2 throw(ArgumentError("n < 2!")) end
    Σxy = zero(Float64)
    Σx  = zero(Float64)
    Σy  = zero(Float64)
    Σx2 = zero(Float64)
    Σy2 = zero(Float64)
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

    return a, b, r2, ar, n
end #end slope
=#
#---------------------------------------------------------------------------
function aucpart(t₁, t₂, c₁, c₂, calcm, aftertmax)
    if calcm == :lint || c₁ <= zero(c₁) && c₂ <= zero(c₂)
        auc   =  linauc(t₁, t₂, c₁, c₂)
    elseif calcm == :logt && aftertmax && c₁ > zero(c₁) && c₂ > zero(c₂)
        auc   =  logauc(t₁, t₂, c₁, c₂)
    elseif calcm == :luld &&  c₁ > c₂ > zero(c₂)
        auc   =  logauc(t₁, t₂, c₁, c₂)
    elseif calcm == :luldt && aftertmax && c₁ > c₂ > zero(c₂)
        auc   =  logauc(t₁, t₂, c₁, c₂)
    #elseif calcm == :log && c₁ > zero(T) && c₂ > zero(T)
        #auc   =  logauc(t₁, t₂, c₁, c₂)
    else
        auc   =  linauc(t₁, t₂, c₁, c₂)
    end
    return auc
end
function aumcpart(t₁, t₂, c₁, c₂, calcm, aftertmax)
    if calcm == :lint || c₁ <= zero(c₁) && c₂ <= zero(c₂)
        aumc  = linaumc(t₁, t₂, c₁, c₂)
    elseif calcm == :logt && aftertmax && c₁ > zero(c₁) && c₂ > zero(c₂)
        aumc  = logaumc(t₁, t₂, c₁, c₂)
    elseif calcm == :luld &&  c₁ > c₂ > zero(c₂)
        aumc  = logaumc(t₁, t₂, c₁, c₂)
    elseif calcm == :luldt && aftertmax && c₁ > c₂ > zero(c₂)
        aumc  = logaumc(t₁, t₂, c₁, c₂)
    #elseif calcm == :log && c₁ > zero(T) && c₂ > zero(T)
        #aumc  = logaumc(t₁, t₂, c₁, c₂)
    else
        aumc  = linaumc(t₁, t₂, c₁, c₂)
    end
    return aumc
end
#---------------------------------------------------------------------------
function interpolate(t₁, t₂, tx, c₁, c₂, intpm, aftertmax)
    if intpm == :lint || c₁ <= zero(c₁) || c₂ <= zero(c₂)
        c = linpredict(t₁, t₂, tx, c₁, c₂)
    elseif intpm == :logt && aftertmax && c₁ > zero(c₁) && c₂ > zero(c₂)
        c = logcpredict(t₁, t₂, tx, c₁, c₂)
    elseif intpm == :luld && c₁ > c₂ > zero(c₂)
        c = logcpredict(t₁, t₂, tx, c₁, c₂)
    elseif intpm == :luldt && aftertmax && c₁ > c₂ > zero(c₂)
        c = logcpredict(t₁, t₂, tx, c₁, c₂)
    #elseif intpm == :log && c₁ > zero(T) && c₂ > zero(T)
        #c = logcpredict(t₁, t₂, tx, c₁, c₂)
    else
        c = linpredict(t₁, t₂, tx, c₁, c₂)
    end
    return c
end

# Time interpolation
function tinterpolate(c₁, c₂, cx, t₁, t₂, intpm, aftertmax)
    if intpm == :lint || c₁ <= zero(c₁) || c₂ <= zero(c₂)
        t = linpredict(c₁, c₂, cx, t₁, t₂)
    elseif intpm == :logt && aftertmax && c₁ > zero(c₁) && c₂ > zero(c₂)
        t = logtpredict(c₁, c₂, cx, t₁, t₂)
    elseif intpm == :luld && c₁ > c₂ > zero(c₂)
        t = logtpredict(c₁, c₂, cx, t₁, t₂)
    elseif intpm == :luldt && aftertmax && c₁ > c₂ > zero(c₂)
        t = logtpredict(c₁, c₂, cx, t₁, t₂)
    #elseif intpm == :log && c₁ > zero(T) && c₂ > zero(T)
        #t = logtpredict(c₁, c₂, cx, t₁, t₂)
    else
        t = linpredict(c₁, c₂, cx, t₁, t₂)
    end
    return t
end

################################################################################
# STEPS
################################################################################
# 1
# Make observation vector and time vector, no points befor Dosetime and nopoints after last nonmissing concentration
function step_1_filterpksubj(_time::AbstractVector{T1}, _obs::AbstractVector{T2}, dosetime) where T1 where T2
    fobs     = firstobs(_time, _obs, dosetime)
    li       = findlast(!isnanormissing, _obs)
    n        = li - fobs + 1
    obstype  = typeof(zero(T2))
    nanobst  = NaN * zero(T2)

    time     = Vector{T1}(undef, n)
    obs      = Vector{obstype}(undef, n)
    inds     = Vector{Int}(undef, 0)

    ii       = 0

    for i = fobs:li
        ii += 1
        time[ii] = _time[i]
        if isnanormissing(_obs[i])
            obs[ii] = nanobst
            push!(inds, ii)
        else
            obs[ii] = _obs[i]
        end
    end
    time, obs, inds # return time, obs, and NaN or Missing value list
end
#2

function step_2_interpolate!(time, obs::AbstractVector{T}, inds, tmaxn, intpm) where T
    if length(inds) > 0
        vals = Vector{T}(undef, length(inds))
        for i = 1:length(inds)
            aftertmax = inds[i] > tmaxn
            i₁ = findlast(!isnanormissing, obs[1:inds[i] - 1])
            i₂ = findfirst(!isnanormissing, obs[inds[i] + 1:end]) + inds[i]
            vals[i] = interpolate(time[i₁], time[i₂], time[inds[i]], obs[i₁], obs[i₂], intpm, aftertmax)
        end
        for i = 1:length(inds)
            obs[inds[i]] = vals[i]
        end
    end
end


function getexcltimes(time::AbstractVector{T}, excl::ElimRange{:time}) where T
    if length(excl.kelexcl) > 0
        inds = findall(x -> x in excl.kelexcl, time)
        if isnothing(inds)
            return T[], Int[] 
        else
            return time[inds], inds
        end
    else
        return T[], Int[]
    end
end
function getexcltimes(time, excl::ElimRange{:point})
    return time[excl.kelexcl], excl.kelexcl
end

function getstimep(::Any, time_cp, excl::ElimRange{:time})
    return findfirst(x -> x >= excl.kelstart, time_cp)
end
function getstimep(time, time_cp, excl::ElimRange{:point}) 
    return findfirst(x -> x >= time[excl.kelstart], time_cp)
end
function getetimep(::Any, time_cp,  excl::ElimRange{:time})
    return findlast(x -> x <= excl.kelend, time_cp)
end
function getetimep(time, time_cp, excl::ElimRange{:point}) 
    return findlast(x -> x <= time[excl.kelend], time_cp)
end
# 3
# Elimination, TlastN, Tlast
function step_3_elim!(result, data, adm, tmaxn, time_cp::AbstractVector{T}, obs_cp::AbstractVector{O}, keldata, kelauto) where T where O
    resize!(keldata)

    time   = gettime(data)
    obsnum = length(time_cp)
    # data.kelrange.kelexcl - indexes; excltime - time values
    excltime, _ = getexcltimes(time, data.kelrange)
    # Unitful values
    r_time_cp = reinterpret(typeof(one(T)), time_cp)
    r_obs_cp  = reinterpret(typeof(one(O)), obs_cp)

    tlastn   = findlast(x-> x > zero(x), r_obs_cp)
    tlast    = r_time_cp[tlastn]

    if kelauto
        if (adm != :iv && obsnum - tmaxn > 2) || (adm == :iv && obsnum - tmaxn > 1)
            if adm == :iv
                stimep = tmaxn
            else
                stimep = tmaxn + 1
            end
            timep = collect(stimep:obsnum) # time points (and conc) - indexes for time vector from start to end
            # Esclude all indexes in data.kelrange.kelexcl
            if length(data.kelrange.kelexcl) > 0
                exclinds = findall(x -> x in excltime, time_cp)
                filter!(x-> x ∉ exclinds, timep)
            end
            # find all concentrations <= 0 - indexes
            zcinds = findall(x -> x <= zero(O), obs_cp)
            # exclude concentrations <= 0 from time vector
            filter!(x-> x ∉ zcinds, timep)
            if length(timep) > 2
                logconc    = log.(r_obs_cp)
                for i = length(timep)-2:-1:1
                    timepv = view(timep, i:length(timep))
                    sl = slope(view(r_time_cp, timepv), view(logconc, timepv))
                    if sl[1] < 0
                        push!(keldata, time_cp[timep[i]], time_cp[timep[end]], sl[1], sl[2], sl[3], sl[4], sl[5])
                    end
                end
            end
        end
    else
        stimep = getstimep(time, time_cp, data.kelrange) 
        if isnothing(stimep)
            @warn "Start-time not found - automatic kel calclation used."
            step_3_elim!(result, data, adm, tmaxn, time_cp, obs_cp, keldata, kelauto)
        end
        etimep = getetimep(time, time_cp, data.kelrange) 
        if isnothing(stimep)
            @warn "End-time not found - automatic kel calclation used."
            step_3_elim!(result, data, adm, tmaxn, time_cp, obs_cp, keldata, kelauto) 
        end
        timep = collect(stimep:etimep)
        if length(data.kelrange.kelexcl) > 0
            @inbounds for i in data.kelrange.kelexcl
                filter!(x-> x ∉ findall(x -> x in excltime, time_cp), timep)
            end
        end
        zcinds = findall(x -> x <= 0, obs_cp)
        filter!(x-> x ∉ zcinds, timep)
        if length(timep) > 1
            sl = slope(view(r_time_cp, timep), view(r_obs_cp, timep), funk = log)
            push!(keldata, time_cp[stimep], time_cp[etimep], sl[1], sl[2], sl[3], sl[4], sl[5])
        end
    end
    # keldata - New KelData, excltime excluded time values, t last N and tlast
    keldata, excltime, tlastn, tlast
end
# 6
function step_6_areas(time_cp::AbstractVector{T}, obs_cp::AbstractVector{O}, calcm, tmaxn, tlastn) where T where O
    obsnum     = length(time_cp)
    auctype    = typeof(zero(eltype(time_cp))*zero(eltype(obs_cp)))
    aumctype   = typeof(zero(eltype(time_cp))^2*zero(eltype(obs_cp)))
    aucpartl   = Array{auctype, 1}(undef, obsnum - 1)
    aumcpartl  = Array{aumctype, 1}(undef, obsnum - 1)
    #Calculate all AUC/AUMC part based on data
    for i = 1:(obsnum - 1)
        aucpartl[i]  = aucpart(time_cp[i], time_cp[i + 1], obs_cp[i], obs_cp[i + 1], calcm, i >= tmaxn)
        aumcpartl[i] = aumcpart(time_cp[i], time_cp[i + 1], obs_cp[i], obs_cp[i + 1], calcm, i >= tmaxn)
    end

    #-----------------------------------------------------------------------
    #-----------------------------------------------------------------------
    auclast  = zero(auctype)
    aumclast = zero(aumctype)
    @inbounds for i = 1:tlastn-1
        auclast  += aucpartl[i]
        aumclast += aumcpartl[i]
    end
    if auclast == zero(auctype) auclast = NaN * zero(auctype) end
    if aumclast == zero(aumctype) aumclast = NaN * zero(aumctype) end
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

    nca(data, time, obs, sort; kelauto = true,  elimrange = ElimRange(), dosetime = DoseTime(), kwargs...)

    nca(data, time, obs; kelauto = true,  elimrange = ElimRange(), dosetime = DoseTime(), kwargs...)

    nca(time, obs; kelauto = true,  elimrange = ElimRange(), dosetime = DoseTime(), kwargs...)

Import data and perform NCA analysis.

Syntax simillar to [`pkimport`](@ref)

Applicable `kwargs` see  [`nca!`](@ref).

See also: [`ElimRange`](@ref), [`DoseTime`](@ref), [`LimitRule`](@ref).
"""
function nca(args...; type::Symbol = :bps, bl = 0, th = 0, kelauto = true,  elimrange = ElimRange(), dosetime = DoseTime(), limitrule::Union{Nothing, LimitRule} = nothing, io = stdout, wm = false, kwargs...)
    if !(type in (:bps, :ur, :pd)) error("Unknown type") end
    if type == :bps
        pki    = pkimport(args...; kelauto = kelauto,  elimrange = elimrange, dosetime = dosetime, limitrule = limitrule, kwargs...)
    elseif type == :ur
        pki    = upkimport(args...; kelauto = kelauto,  elimrange = elimrange, dosetime = dosetime, kwargs...)
    elseif type == :pd
        pki    = pdimport(args...; th = th, bl = bl, limitrule = limitrule, kwargs...)
    end
    nca!(pki; io = io, wm = wm, kwargs...)
end

"""
    nca!(data::DataSet{Subj}, obs::Union{Symbol, Nothing} = nothing; adm = :ev, calcm = :lint, intpm = nothing, verbose = 0, warn = true, wm::Bool = false, io::IO = stdout, modify! = identity) where Subj <: AbstractSubject

Non-compartmental (NCA) analysis of PK/PD data.

**Arguments**

* data - DataSet to analyze;
* obs - observation to use as a PK metric (by default used first observations vector).

 **Keywords**

* `wm` - write output to DataSet metadata (`:ncalog` index).

Other keywords is the same as in `nca!` method for `PKSubject`.
"""
function nca!(data::DataSet{Subj}, obs::Union{Symbol, Nothing} = nothing;  wm::Bool = false, io = stdout, kwargs...) where Subj <: AbstractSubject

    writetostdout = false
    if wm
        if io == stdout
            io = IOBuffer()
            writetostdout = true
        end
    end

    ds = map(x -> nca!(x, obs; io = io, kwargs...), data)
    
    if wm
        io2 = deepcopy(io)
        ds.metadata[:ncalog] = String(take!(io2))
        if writetostdout 
            write(stdout, String(take!(io))) 
        end
    end
    ds
end

"""
    nca!(data::PKSubject{T,O}, obs::Union{Symbol, Nothing} = nothing; adm = :ev, calcm = :lint,  usedose::Int = 0, intpm = nothing,  partials = nothing, prtext = :err, verbose = 0, warn = true, io::IO = stdout, modify! = nothing) where T where O

**Arguments**

* data - PK Subject to analyze;
* obs - observation to use as a PK metric (by default used first observations vector).

**Keywords**

* `adm` - administration:
    - `:ev` - extra vascular;
    - `:iv` - intravascular bolus;
* `calcm` - AUC/AUMC calculation method:
    - `:lint` - linear trapezoidal;
    - `:luld` - linear up log down;
    - `:luldt` - linear up log down after Tmax;
    - `:logt` - log-trapezoidal after Tmax (Not Recommended);
* `usedose` - number of dose used for calculation; if `usedose = 0` (by default) - last dose used.
* `intpm` - interpolation method:
    - `:lint` - linear trapezoidal;
    - `:luld` - linear up log down;
    - `:luldt` - linear up log down after Tmax;
    - `:logt` - log-trapezoidal after Tmax;
    - `:calcm` - same as `calcm`;
* `partials` - calculate partial AUC vor vector of time intervals (`:err` (default) - throw error if end time > last oservation time; `:last` - no extrapolation; `:extr` - if `Kel` calculated used extrapolation or `NaN` if no `Kel`);
* `prtext` - extrapolation rule for partials AUC;
* `verbose` - print to `io`, 1: partial areas table, 2: 1, and results;
* `warn` - show warnings;
* `io` - output stream;
* `modify!` - function to modify output paramaters, call `modify!(ncaresult)` if difined.

Results:

* Cmax
* Tmax
* Cdose
* Tlag
* Clast
* AUClast
* AUMClast
* AUCall
* Rsq
* ARsq
* Kel
* HL
* LZint
* NpLZ
* Clast_pred
* AUCinf
* AUCinf_pred
* AUMCinf
* AUMCinf_pred
* AUCpct
* MRTlast
* MRTinf
* MRTinf_pred
* Cllast
* Clinf
* Vzlast
* Vzinf
* Vssinf

Steady-state parameters (tau used):

* AUCtau
* AUMCtau
* Ctau
* Cavg
* Ctaumax
* Ctaumin
* Accind
* Fluc
* Fluctau
* Swing
* Swingtau
* MRTtauinf
* MRTtauinf_pred
* Cltau
* Vztau

`partials` is a vector of vectors, tuples or pairs. Example: `partials = [(1,2), (3,4)]`, `partials = [[1,2], (3,4)]`

"""
function nca!(data::PKSubject{T, OBS}, obs::Union{Symbol, Nothing} = NCARESOBS; 
    adm = :ev, calcm = :lint, usedose::Int = 0,
    intpm = :calcm,  
    partials = nothing, 
    prtext = :err, 
    verbose = 0, 
    warn = true, 
    io::IO = stdout, 
    modify! = identity,
    kwargs...) where T where OBS

    O = eltype(getobs(data, obs))

    ptype  = promote_type(Float64, T, O)

    result   = Dict{Symbol, ptype}()

    if intpm == :calcm intpm = calcm end

    options =  Dict(:obsname => obs, :type => :bps, :adm => adm, :calcm => calcm, :intpm => intpm, :verbose => verbose, :warn => warn, :modify! => modify!)

    if usedose == 0
        usedose = length(data.dosetime)
    else
        if usedose < 0 || usedose > length(data.dosetime) error("Wrong 'usedose' value.") end
    end

    used_dosetime_tau  = data.dosetime[usedose].tau
    used_dosetime_time = data.dosetime[usedose].time
    used_dosetime_dose = data.dosetime[usedose].dose


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
        println(io, "    Method: $(calcm); Dose: $(used_dosetime_dose); Dose time: $(used_dosetime_time)")
        if used_dosetime_tau > 0
            println(io, "    Tau: $(used_dosetime_tau)")
        end
    end
################################################################################
    # STEP 1 FILTER ALL BEFORE DOSETIME AND ALL NAN OR MISSING VALUES
    if validobsn(gettime(data), getobs(data, obs)) == 0
        return NCAResult(data, options, result) 
    end
    time_cp, obs_cp, einds = step_1_filterpksubj(gettime(data), getobs(data, obs), used_dosetime_time)
    if length(obs_cp) < 2
        return NCAResult(data, options, result)
    end
################################################################################
    # STEP 2 - CMAX TMAX FOR TAU RANGE Clast Tlast; interpolate NaN and missings
    #result[:Obsnum] = obsnum = length(obs_cp)
    result[:Obsnum] = validobsn(time_cp, obs_cp)
    # If TAU set, calculates start and end timepoints for AUCtau

    if  used_dosetime_tau > zero(typeof(used_dosetime_tau))
        taulastp = findlast(x -> x <= used_dosetime_time + used_dosetime_tau, time_cp)
        result[:Cmax], result[:Tmax], tmaxn, result[:Ctaumin] = ctmax(time_cp, obs_cp, taulastp)
    else
        taulastp = length(obs_cp)
        result[:Cmax], result[:Tmax], tmaxn, _cmin = ctmax(time_cp, obs_cp)
    end
    

    step_2_interpolate!(time_cp, obs_cp, einds, tmaxn, intpm)

################################################################################
    # STEP 3
    # Elimination, add interpolated inds to elimination exclusion
    # Get last concentration
    if length(einds) > 0
        for ei in einds
            if ei ∉ data.kelrange.kelexcl push!(data.kelrange.kelexcl, ei) end
        end
        sort!(data.kelrange.kelexcl)
        if data.kelrange.kelstart in data.kelrange.kelexcl || data.kelrange.kelend in data.kelrange.kelexcl
            data.kelauto = true
        end
    end


    keldata, excltime, tlastn, tlast = step_3_elim!(result, data, adm, tmaxn, time_cp, obs_cp, KelData(T[], T[], Float64[], Float64[], Float64[], Float64[], Int[]), data.kelauto)
    # C last and T last
    result[:Tlast]   = time_cp[tlastn]
    result[:Clast]   = obs_cp[tlastn]


################################################################################
    # STEP 4
    if  used_dosetime_tau > zero(used_dosetime_tau)
        time_cp .-= used_dosetime_time
    end
################################################################################
    # STEP 5
    # Dose concentration
    # Dosetime is first point
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
            if  used_dosetime_tau > zero(typeof(used_dosetime_tau))
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
    if used_dosetime_dose > zero(used_dosetime_dose)
        result[:Cllast]           = used_dosetime_dose / result[:AUClast]
        result[:Dose]             = used_dosetime_dose
    end
    #-----------------------------------------------------------------------
    #-----------------------------------------------------------------------
    tlagn = findfirst(!iszero, obs_cp)

    if isnothing(tlagn)
        result[:Tlag] = NaN*zero(T)
    elseif tlagn > 1
        result[:Tlag] = time_cp[tlagn-1]
    else
        result[:Tlag] = zero(T)
    end

    if  length(keldata) > 0
        #data.keldata             = keldata
        result[:ARsq], rsqn      = findmax(keldata.ar)
        result[:Rsq]             = keldata.r[rsqn]
        #kel                      = abs(keldata.a[rsqn])
        result[:Kel]             = abs(keldata.a[rsqn]) / oneunit(T)
        result[:LZ]              = keldata.a[rsqn]
        result[:LZint]           = keldata.b[rsqn]
        result[:Rsqn]            = rsqn
        result[:NpLZ]            = keldata.n[rsqn]
        result[:Clast_pred]      = exp(result[:LZint] + result[:LZ] * tlast) * oneunit(O)
        result[:HL]              = LOG2 / result[:Kel]
        result[:AUCinf]          = result[:AUClast] + result[:Clast] / result[:Kel]
        result[:AUCinf_pred]     = result[:AUClast] + result[:Clast_pred] / result[:Kel]
        result[:AUCpct]          = (result[:AUCinf] - result[:AUClast]) / result[:AUCinf] * 100
        result[:AUMCinf]         = result[:AUMClast] + result[:Tlast] * result[:Clast] / result[:Kel] + result[:Clast] / result[:Kel] ^ 2
        result[:AUMCinf_pred]    = result[:AUMClast] + result[:Tlast] * result[:Clast_pred] / result[:Kel] + result[:Clast_pred] / result[:Kel] ^ 2
        result[:MRTinf]          = result[:AUMCinf] / result[:AUCinf]
        result[:MRTinf_pred]     = result[:AUMCinf_pred] / result[:AUCinf_pred]
        if used_dosetime_dose > zero(used_dosetime_dose)
            result[:Vzlast]          = used_dosetime_dose / result[:AUClast] / result[:Kel]
            result[:Vzinf]           = used_dosetime_dose / result[:AUCinf] / result[:Kel]
            result[:Clinf]           = used_dosetime_dose / result[:AUCinf]
            result[:Vssinf]          = result[:Clinf] * result[:MRTinf]
        end
    else
        result[:Kel] = NaN / oneunit(T)
    end
################################################################################
    # STEP 8
    # Steady-state parameters
    if used_dosetime_tau > zero(used_dosetime_tau)
        eaucpartl  = zero(T)*zero(O)
        eaumcpartl = zero(T)^2*zero(O)
        if time_cp[taulastp] < used_dosetime_tau < time_cp[end]
            result[:Ctau] = interpolate(time_cp[taulastp], time_cp[taulastp + 1], used_dosetime_tau, obs_cp[taulastp], obs_cp[taulastp + 1], intpm, true)
            eaucpartl = aucpart(time_cp[taulastp], used_dosetime_tau, obs_cp[taulastp], result[:Ctau], calcm, true)
            eaumcpartl = aumcpart(time_cp[taulastp], used_dosetime_tau, obs_cp[taulastp], result[:Ctau], calcm, true)
                #remoove part after tau
        elseif used_dosetime_tau > time_cp[end] && result[:Kel] !== NaN
                #extrapolation
            result[:Ctau] = exp(result[:LZint] + result[:LZ] * (used_dosetime_tau + used_dosetime_time))
            eaucpartl = aucpart(time_cp[end], used_dosetime_tau, obs_cp[end], result[:Ctau], calcm, true)
            eaumcpartl = aumcpart(time_cp[end], used_dosetime_tau, obs_cp[end], result[:Ctau], calcm, true)
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

        result[:Cavg]     = result[:AUCtau]/used_dosetime_tau
        if result[:Ctaumin] != 0
            result[:Swing]    = (result[:Cmax] - result[:Ctaumin])/result[:Ctaumin]
        end
        if result[:Ctau] != 0
            result[:Swingtau] = (result[:Cmax] - result[:Ctau])/result[:Ctau]
        end
        result[:Fluc]     = (result[:Cmax] - result[:Ctaumin])/result[:Cavg] * 100
        result[:Fluctau]  = (result[:Cmax] - result[:Ctau])/result[:Cavg] * 100
        #If Kel calculated
        result[:Cltau]           = used_dosetime_dose / result[:AUCtau]
        if !isnan(result[:Kel])
            result[:Accind]   = 1 / (1 - (exp(-result[:Kel] * used_dosetime_tau)))
            result[:MRTtauinf]       = (result[:AUMCtau] + used_dosetime_tau * (result[:AUCinf] - result[:AUCtau])) / result[:AUCtau]
            result[:MRTtauinf_pred]  = (result[:AUMCtau] + used_dosetime_tau * (result[:AUCinf_pred] - result[:AUCtau])) / result[:AUCtau]
            result[:Vztau]           = used_dosetime_dose / result[:AUCtau] / result[:Kel]
            result[:Vsstau]          = result[:Cltau] * result[:MRTtauinf]
        end


    end
    #partials
    if !isnothing(partials)
        for prt in partials
            stime = prt[1]
            etime = prt[2]
            if stime <  used_dosetime_time error("Start time can't be less than dose time!") end
            if stime <  used_dosetime_time error("End time can't be less than dose time!") end
            if etime <= stime error("End time can't be less or equal start time!") end
            suffix = "_"*string(stime)*"_"*string(etime)
            stime = stime - used_dosetime_time
            etime = etime - used_dosetime_time
            if etime > last(time_cp) && prtext == :err error("End time can't be greater than last time point ($(last(time_cp)))! Use keyword `prtext=:last` or `prtext=:extr`...") end 
            #first point
            firstp = findfirst(x -> x >= stime, time_cp)
            #last point
            lastp  = findlast(x -> x <= etime, time_cp)
            firstpart = zero(T)*zero(O)
            lastpart  = zero(T)*zero(O)
            if stime < time_cp[firstp]
                firstpartc = interpolate(time_cp[firstp - 1], time_cp[firstp], stime, obs_cp[firstp - 1], obs_cp[firstp], intpm,  stime > result[:Tmax])
                firstpart += aucpart(stime, time_cp[firstp], firstpartc, obs_cp[firstp], calcm, stime > result[:Tmax])
                #println("firstpartc = $firstpartc , firstpart = $firstpart")
            end
            if etime > time_cp[lastp] && etime < last(time_cp) # if last time > etime -> interpolation
                lastpartc  = interpolate(time_cp[lastp], time_cp[lastp + 1], etime, obs_cp[lastp], obs_cp[lastp + 1], intpm,  time_cp[lastp] > result[:Tmax])
                lastpart +=  aucpart(time_cp[lastp], etime, obs_cp[lastp], lastpartc, calcm, time_cp[lastp] > result[:Tmax])
                #println("lastpartc = $lastpartc , lastpart = $lastpart")
            elseif etime >= time_cp[lastp] && prtext == :last
                lastpartc = zero(O)
            elseif etime > time_cp[lastp] && prtext == :extr && !isnan(result[:Kel])
                lastpartc = exp(result[:LZint] + result[:LZ] * etime)
                lastpart += aucpart(time_cp[lastp], etime, obs_cp[end], lastpartc, calcm, time_cp[lastp] > result[:Tmax])
            else
                lastpartc = NaN
                lastpart += lastpartc
            end


            aucpartial = zero(T)*zero(O)
            if firstp != lastp
                aucpartn = lastp - firstp
                for i = 1:aucpartn
                    aucpartial += aucpart(time_cp[firstp + i - 1], time_cp[firstp + i], obs_cp[firstp + i - 1], obs_cp[firstp + i], calcm, time_cp[firstp + i - 1] >= result[:Tmax])
                    #println("first = $(firstp + i - 1) , scns = $(firstp + i) , aucpartial = $aucpartial")
                end
            end
            aucpartial += firstpart + lastpart
            result[Symbol("AUC"*suffix)] = aucpartial
            if verbose > 2
                println(io, "    Partial $(prt[1]) - $(prt[2]); from point $firstp to $lastp")
                if stime < time_cp[firstp]
                    println(io, "      Interpolation values: first = $firstpartc, last = $lastpartc")
                end
                if etime > time_cp[lastp]
                    println(io, "      Interpolation parts: first = $firstpart, last = $lastpart")
                end
            end
        end
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
        if  used_dosetime_time > 0
            time_cp .+= used_dosetime_time
        end
        hnames = [:Time, :Concentrtion, :AUC, :AUC_cum, :AUMC, :AUMC_cum, :Info]
        mx = metida_table(time_cp, obs_cp, pushfirst!(aucpartl, 0.0),  pushfirst!(aucpartlsum, 0.0), pushfirst!(aumcpartl, 0.0),  pushfirst!(aumcpartlsum, 0.0), fill("", length(obs_cp));
        names = hnames)
        if cdoseins > 0
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
                if i in einds
                    mx[i, 7] *= "@"
                end
            end
        end
        hnames = (["Time" "Conc." "AUC"  "AUC" "AUMC" "AUMC" "Info"],
                  ["" "" "" "(cum.)" "" "(cum.)" ""])
        PrettyTables.pretty_table(io, mx; tf = PrettyTables.tf_compact, header = hnames, formatters = PrettyTables.ft_printf("%3.4g"))
        println(io, "")
        println(io, "    Cdose: $(result[:Cdose]), Dose time: $(used_dosetime_time)")
        if isnan(result[:Kel])
            println(io, "    Elimination not calculated")
        else
            println(io, "    Kel start: $(keldata.s[rsqn]); end: $(keldata.e[rsqn])")
        end
        if length(einds) > 0
            println(io, "    @ - Interpolated points ($(length(einds)))")
        end
        println(io, "")
        if used_dosetime_tau < time_cp[end] && used_dosetime_tau > 0
            println(io, "    Tau + dosetime is less then end time. Interpolation used.")
            println(io, "    Ctau: $(result[:Ctau])")
            println(io, "    AUC  final part: $(eaucpartl)")
            println(io, "    AUMC final part: $(eaumcpartl)")
            println(io, "")
        end

        if verbose > 1
            println(io, "    Results:")
            PrettyTables.pretty_table(io, result; tf = PrettyTables.tf_compact, header = ["Parameter", "Value"], formatters = PrettyTables.ft_printf("%4.6g"))
        end
    end
################################################################################

    ncares = NCAResult(data, options, obs, keldata, result)
    modify!(ncares)

    #-----------------------------------------------------------------------
    return ncares
end

function maxconc(subj::T; obsname::Union{Symbol, Nothing} = nothing) where T <: AbstractSubject
    maximum(Iterators.filter(x-> !(ismissing(x) || isnan(x)), getobs(subj, obsname)))
end
function minconc(subj::T, pos = false; obsname::Union{Symbol, Nothing} = nothing) where T <: AbstractSubject
    if pos
        return minimum(Iterators.filter(x-> x > zero(x), getobs(subj, obsname)))
    else
        return minimum(Iterators.filter(x-> !(ismissing(x) || isnan(x)), getobs(subj, obsname)))
    end
end

function exrate(time::AbstractVector{Tuple{S, E}}, conc::AbstractVector{C}, vol::AbstractVector{V})  where S where E where C where V
    T = promote_type(S, E)
    length(time) == length(conc) == length(vol) || error("")
    er = Vector{typeof(oneunit(C)*oneunit(V)/oneunit(T))}(undef, length(time))
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

"""
    nca!(data::UPKSubject{T, O, VOL, V}; adm = :ev, calcm = :lint, intpm = nothing, verbose = 0, warn = true, io::IO = stdout, modify! = identity) where T where O where VOL where V

Non-compartmental (NCA) analysis of pharmacokinetic for urine data.

Results:

* AUCall
* AUClast
* Rlast
* Maxrate
* Tmax
* AR
* Vol
* Prec
* ARsq
* Rsq
* Kel
* LZ
* LZint
* Rsqn
* HL
* AUCinf
"""
function nca!(data::UPKSubject{Tuple{S, E}, O, VOL, V}, obs::Union{Symbol, Nothing} = nothing; adm = :ev, calcm = :lint, intpm = nothing, verbose = 0, warn = true, io::IO = stdout, modify! = identity, kwargs...) where S where E where O where VOL where V

    ptype  = promote_type(Float64, S, E, O, VOL)
    ttype  = promote_type(S, E)
    result   = Dict{Symbol, ptype}()
    
    options =  Dict(:type => :urine, :adm => adm, :calcm => calcm, :intpm => intpm, :verbose => verbose, :warn => warn, :modify! => modify!)

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


    time, obsvals, vol = step_1_filterupksubj(gettime(data), getobs(data, obs), data.vol, data.dosetime.time)

    mtime  = map(x-> (x[1]+x[2])/2, time)

    exr  = exrate(time, obsvals, vol)

    result[:AR]   = getobs(data, obs)' * data.vol
    result[:Vol]  = sum(vol)

    if time[1][1] > data.dosetime.time
        pushfirst!(mtime, 0)
    else
        pushfirst!(mtime, time[1][1])
    end
    pushfirst!(obsvals, zero(O))
    pushfirst!(vol, zero(VOL))
    pushfirst!(exr, zero(eltype(exr)))

    result[:Maxrate], result[:Tmax], tmaxn, result[:Minrate] = ctmax(mtime, exr)

    if data.dosetime.dose > zero(data.dosetime.dose)
        result[:Prec]  = result[:AR]/data.dosetime.dose * 100
    end
    obsnum = length(exr)

    lastobs = length(exr)
    for i = length(exr):-1:1
        if exr[i] > zero(eltype(exr))
            break
        else
            lastobs = i
        end
    end

    # STEP 3
    # Elimination
    keldata, excltime = step_3_elim!(result, data, adm, tmaxn, mtime, exr, KelData(ttype[], ttype[], Float64[], Float64[], Float64[], Float64[], Int[]), data.kelauto)

    #result[:Kel]
    #result[:HL]
    aucpartl  = Vector{typeof(zero(eltype(exr))*zero(ttype))}(undef, obsnum - 1)
    #Calculate all AUC part based on data
    for i = 1:(obsnum - 1)
        aucpartl[i] = aucpart(mtime[i], mtime[i + 1], exr[i], exr[i + 1], calcm, i >= tmaxn)
    end

    result[:AUCall]  = sum(aucpartl)
    result[:AUClast] = sum(aucpartl[1:lastobs-1])
    result[:Rlast]  = exr[end]

    if  length(keldata) > 0
        result[:ARsq], rsqn      = findmax(keldata.ar)
        result[:Rsq]             = keldata.r[rsqn]
        result[:Kel]             = abs(keldata.a[rsqn]) / oneunit(ttype)
        result[:LZ]              = keldata.a[rsqn]
        result[:LZint]           = keldata.b[rsqn]
        result[:Rsqn]            = rsqn
        result[:NpLZ]            = keldata.n[rsqn]
        result[:HL]              = LOG2 / result[:Kel]

        result[:AUCinf]          = result[:AUClast] + result[:Rlast] / result[:Kel]
        result[:AUCpct]          = (result[:AUCinf] - result[:AUClast]) / result[:AUCinf] * 100
    end

    ncares = NCAResult(data, options, obs, keldata, result)
    modify!(ncares)
    #-----------------------------------------------------------------------
    return ncares
end


function auctspl(c1, c2, t1, t2, sl, calcm)
    slu = sl*oneunit(c1)
    if c1 >= slu && c2 >= slu
        auca = aucpart(t1, t2, c1, c2, calcm, true) - (t2 - t1)*slu
        aucb = zero(auca)
        ta   = t2 - t1
        tb   = zero(ta)
    elseif c1 <= slu && c2 <= slu
        aucb = (t2 - t1)*slu - aucpart(t1, t2, c1, c2, calcm, true)
        auca = zero(aucb)
        tb   = t2 - t1
        ta   = zero(tb)
    else
        tint = tinterpolate(c1, c2, slu, t1, t2, calcm, true)
        l = aucpart(t1, tint, c1, slu, calcm, true)
        r = aucpart(tint, t2, slu, c2, calcm, true)
        if c1 > slu && c2 < slu
            auca = l - (tint - t1)*slu
            aucb = (t2 - tint)*slu - r
            ta   = tint - t1
            tb   = t2 - tint
        elseif c1 < slu && c2 > slu
            auca = r - (t2 - tint)*slu
            aucb = (tint - t1)*slu - l
            ta   = t2 - tint
            tb   = tint - t1
        else
            error("!!")
        end
    end
    auca, aucb, ta, tb
end

function auctblth(c1, c2, t1, t2, bl, th, calcm)

    aucabl, aucbbl, tabl, tbbl = auctspl(c1, c2, t1, t2, bl, calcm)

    aucath, aucbth, tath, tbth = auctspl(c1, c2, t1, t2, th, calcm)

    if th > bl
        aucbtw = aucabl - aucath
    else
        aucbtw = aucbbl - aucbth
    end

    aucabl, aucbbl, tabl, tbbl, aucath, aucbth, tath, tbth, aucbtw
end

"""
    nca!(data::PDSubject{T,O}; calcm = :lint, intpm = nothing, verbose = 0, warn = true, io::IO = stdout, modify! = identity, kwargs...) where T where O

Non-compartmental (NCA) analysis of pharmacodynamic data.

Results:

* Rmax - max responce;
* Tmax - time for maximum responce;
* AUCABL - AUC above baseline;
* AUCBBL - AUC below baseline;
* AUCATH - AUC above threshold;
* AUCBTH - AUC below threshold;
* AUCNETB - AUCABL - AUCBBL;
* AUCNETT - AUCATH - AUCBTH;
* TABL - time above baseline;
* TBBL - time below baseline;
* TATH - time above threshold;
* TBTH - time below threshold;
* AUCBTW - AUC between baseline and threshold;

"""
function nca!(data::PDSubject{T,O}, obs::Union{Symbol, Nothing} = nothing; calcm = :lint, intpm = nothing, verbose = 0, warn = true, io::IO = stdout, modify! = identity, kwargs...) where T where O
    
    ptype  = promote_type(Float64, T, O)

    result   = Dict{Symbol, ptype}()

    if isnothing(intpm) intpm = calcm end

    options =  Dict(:type => :pd, :calcm => calcm, :intpm => intpm, :verbose => verbose, :warn => warn, :modify! => modify!)

    if verbose > 0
        println(io, "  Pharmacodynamic Analysis")
        if length(data.id) > 0
            print(io, "    Subject: ")
            for (k,v) in data.id
                print(io, "$(k) => $(v); ")
            end
            println(io, "")
        end
        println(io, "    Settings:")
        println(io, "    Method: $(calcm);")
    end

    auctype  = promote_type(eltype(gettime(data)), eltype(getobs(data)))

################################################################################
    # STEP 1 FILTER ALL BEFORE DOSETIME AND ALL NAN OR MISSING VALUES
    if validobsn(gettime(data), getobs(data)) == 0 return NCAResult(data, options, result) end
    time_cp, obs_cp, einds = step_1_filterpksubj(gettime(data), getobs(data), first(gettime(data)))
    if length(obs_cp) < 2
        return NCAResult(data, options, result)
    end
################################################################################
    result[:Obsnum] = obsnum = length(obs_cp)

    result[:BL] = data.bl
 
    result[:TH] = data.th


    result[:Rmax], result[:Tmax], tmaxn, result[:Rmin] = ctmax(time_cp, obs_cp, length(obs_cp))


    
    # ALL NAN AND MISSING VALUES LINEAR INTERPOLATED
    step_2_interpolate!(time_cp, obs_cp, einds, 1, :lint)

    aucpartabl  = Array{ptype, 1}(undef, obsnum - 1)
    aucpartbbl  = Array{ptype, 1}(undef, obsnum - 1)
    aucpartath  = Array{ptype, 1}(undef, obsnum - 1)
    aucpartbth  = Array{ptype, 1}(undef, obsnum - 1)
    tpartabl    = Array{ptype, 1}(undef, obsnum - 1)
    tpartbbl    = Array{ptype, 1}(undef, obsnum - 1)
    tpartath    = Array{ptype, 1}(undef, obsnum - 1)
    tpartbth    = Array{ptype, 1}(undef, obsnum - 1)
    aucpartbtw  = Array{ptype, 1}(undef, obsnum - 1)

    for i = 1:(obsnum - 1)
        aucpartabl[i], aucpartbbl[i], tpartabl[i], tpartbbl[i], aucpartath[i], aucpartbth[i], tpartath[i], tpartbth[i], aucpartbtw[i] = auctblth( obs_cp[i], obs_cp[i + 1], time_cp[i], time_cp[i + 1], data.bl, data.th, calcm)
    end

    result[:AUCABL] = sum(aucpartabl)
    result[:AUCBBL] = sum(aucpartbbl)
    result[:AUCATH] = sum(aucpartath)
    result[:AUCBTH] = sum(aucpartbth)
    result[:TABL]   = sum(tpartabl)
    result[:TBBL]   = sum(tpartbbl)
    result[:TATH]   = sum(tpartath)
    result[:TBTH]   = sum(tpartbth)
    result[:AUCBTW] = sum(aucpartbtw)

    result[:AUCNETB] = result[:AUCABL] - result[:AUCBBL]
    result[:AUCNETT] = result[:AUCATH] - result[:AUCBTH]

    if data.th > data.bl
        result[:TIMEBTW] = result[:TBTH] - result[:TBBL]
    else
        result[:TIMEBTW] = result[:TBBL] - result[:TBTH]
    end

   
    # Tau parameters
    if data.dosetime.tau > zero(data.dosetime.tau)
        taufirstp = findfirst(x -> x >= data.dosetime.time, time_cp)
        taulastp  = findlast(x -> x <= data.dosetime.time + data.dosetime.tau, time_cp)
        taupn = taulastp - taufirstp

        aucpartabl_tau  = Array{ptype, 1}(undef, taupn)
        aucpartbbl_tau  = Array{ptype, 1}(undef, taupn)
        aucpartath_tau  = Array{ptype, 1}(undef, taupn)
        aucpartbth_tau  = Array{ptype, 1}(undef, taupn)
        tpartabl_tau    = Array{ptype, 1}(undef, taupn)
        tpartbbl_tau    = Array{ptype, 1}(undef, taupn)
        tpartath_tau    = Array{ptype, 1}(undef, taupn)
        tpartbth_tau    = Array{ptype, 1}(undef, taupn)
        aucpartbtw_tau  = Array{ptype, 1}(undef, taupn)

        for i = taufirstp:taulastp-1
            j = i - taufirstp + 1
            aucpartabl_tau[j], aucpartbbl_tau[j], tpartabl_tau[j], tpartbbl_tau[j], aucpartath_tau[j], aucpartbth_tau[j], tpartath_tau[j], tpartbth_tau[j], aucpartbtw_tau[j] = auctblth(obs_cp[i], obs_cp[i + 1], time_cp[i], time_cp[i + 1], data.bl, data.th, calcm)
        end
        # before first responce point
        if data.dosetime.time < time_cp[taufirstp]
            if taufirstp > 1
                # Interpolate responce between observations
                intpmr  = interpolate(time_cp[taufirstp-1], time_cp[taufirstp], data.dosetime.time, obs_cp[taufirstp-1], obs_cp[taufirstp], calcm, true)
                aucpartabl_, aucpartbbl_, tpartabl_, tpartbbl_, aucpartath_, aucpartbth_, tpartath_, tpartbth_, aucpartbtw_ = auctblth(intpmr, obs_cp[taufirstp], data.dosetime.time, time_cp[taufirstp], data.bl, data.th, calcm)
            else # else set first point to zero
                aucpartabl_, aucpartbbl_, tpartabl_, tpartbbl_, aucpartath_, aucpartbth_, tpartath_, tpartbth_, aucpartbtw_ = auctblth(0, obs_cp[taufirstp], data.dosetime.time, time_cp[taufirstp], data.bl, data.th, calcm)
            end
            pushfirst!(aucpartabl_tau, aucpartabl_)
            pushfirst!(aucpartbbl_tau, aucpartbbl_)
            pushfirst!(aucpartath_tau, aucpartath_)
            pushfirst!(aucpartbth_tau, aucpartbth_)
            pushfirst!(tpartabl_tau, tpartabl_)
            pushfirst!(tpartbbl_tau, tpartbbl_)
            pushfirst!(tpartath_tau, tpartath_)
            pushfirst!(tpartbth_tau, tpartbth_)
            pushfirst!(aucpartbtw_tau, aucpartbtw_)
        end
        # after last responce point
        if data.dosetime.time + data.dosetime.tau > time_cp[taulastp]
            if  taulastp < length(time_cp)
                # Interpolate responce between observations
                intpmr  = interpolate(time_cp[taulastp], time_cp[taulastp+1], data.dosetime.time + data.dosetime.tau, obs_cp[taulastp], obs_cp[taulastp+1], calcm, true)

                aucpartabl_, aucpartbbl_, tpartabl_, tpartbbl_, aucpartath_, aucpartbth_, tpartath_, tpartbth_, aucpartbtw_ = auctblth(obs_cp[taulastp], intpmr, time_cp[taulastp], data.dosetime.time + data.dosetime.tau, data.bl, data.th, calcm)
            else # else set first point to zero
                aucpartabl_, aucpartbbl_, tpartabl_, tpartbbl_, aucpartath_, aucpartbth_, tpartath_, tpartbth_, aucpartbtw_ =  auctblth(obs_cp[taulastp], 0, time_cp[taulastp], data.dosetime.time + data.dosetime.tau, data.bl, data.th, calcm)
            end
            push!(aucpartabl_tau, aucpartabl_)
            push!(aucpartbbl_tau, aucpartbbl_)
            push!(aucpartath_tau, aucpartath_)
            push!(aucpartbth_tau, aucpartbth_)
            push!(tpartabl_tau, tpartabl_)
            push!(tpartbbl_tau, tpartbbl_)
            push!(tpartath_tau, tpartath_)
            push!(tpartbth_tau, tpartbth_)
            push!(aucpartbtw_tau, aucpartbtw_)
        end

        result[:AUCABLtau] = sum(aucpartabl_tau)
        result[:AUCBBLtau] = sum(aucpartbbl_tau)
        result[:AUCATHtau] = sum(aucpartath_tau)
        result[:AUCBTHtau] = sum(aucpartbth_tau)
        result[:TABLtau]   = sum(tpartabl_tau)
        result[:TBBLtau]   = sum(tpartbbl_tau)
        result[:TATHtau]   = sum(tpartath_tau)
        result[:TBTHtau]   = sum(tpartbth_tau)
        result[:AUCBTWtau] = sum(aucpartbtw_tau)
    end

    # Verbose output
    if verbose > 0
        hnames = [:Time, :Observation, :AUCABL, :AUCBBL, :AUCATH, :AUCBTH]
        mx = metida_table(collect(time_cp),
        collect(obs_cp),
        pushfirst!(aucpartabl, zero(first(aucpartabl))),
        pushfirst!(aucpartbbl, zero(first(aucpartbbl))),
        pushfirst!(aucpartath, zero(first(aucpartath))),
        pushfirst!(aucpartbth, zero(first(aucpartbth)));
        names = hnames)

        hnames = (["Time" "Obs." "AUCABL"  "AUCBBL" "AUCATH" "AUCBTH"],
                  ["" "" "" "" "" "" ""])
        PrettyTables.pretty_table(io, mx; tf = PrettyTables.tf_compact, header = hnames, formatters = PrettyTables.ft_printf("%3.4g"))
        println(io, "")
        if verbose > 1
            println(io, "    Results:")
            PrettyTables.pretty_table(io, result; tf = PrettyTables.tf_compact, header = ["Parameter", "Value"], formatters = PrettyTables.ft_printf("%4.6g"))
        end
    end

    ncares = NCAResult(data, options, result)
    modify!(ncares)
    #-----------------------------------------------------------------------
    return ncares
end
