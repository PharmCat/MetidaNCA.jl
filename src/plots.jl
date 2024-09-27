struct NoPageSort end

const PKPLOTSTYLE = (
(:solid, :blue, :circle, :blue),
(:solid, :red, :utriangle, :red),
(:solid, :green, :diamond, :green),
(:solid, :magenta, :pentagon, :magenta),
(:solid, :purple, :heptagon, :purple),
(:solid, :indigo, :octagon, :indigo),
(:solid, :gold, :star, :gold),
(:solid, :yellow, :rect, :yellow),
(:solid, :gray, :xcross, :gray),
(:solid, :cyan, :cross, :cyan),
(:dot, :blue, :utriangle, :blue),
(:dot, :red, :circle, :red),
(:dot, :green, :rect, :green),
(:dot, :yellow, :diamond, :yellow),
(:dot, :gray, :cross, :gray),
(:dot, :cyan, :xcross, :cyan),
(:dot, :gold, :pentagon, :gold),
(:dot, :magenta, :star, :magenta),
(:dot, :purple, :octagon, :purple),
(:dot, :indigo, :heptagon, :indigo),
)

function plotstyle(n)
    if isnothing(n) n = 1 end
    if n >= 1 && n <= 20
        return PKPLOTSTYLE[n]
    elseif n > 1900
        return (:auto, :auto, :auto, :auto)
    end
    n -= 20
    linestyle   = (:dash, :dashdot, :dashdotdot)
    linecolor   = (:yellow, :gray, :cyan, :gold, :magenta, :purple, :indigo)
    markershape = (:star4, :star6, :star7, :star8, :rect, :star5, :diamond, :hexagon, :cross, :xcross, :utriangle, :dtriangle, :rtriangle, :ltriangle, :pentagon, :heptagon, :octagon, :vline, :hline)
    markercolor = (:gray, :gold, :magenta, :purple, :indigo)
    a = n % 3 + 1
    b = n % 7 + 1
    c = n % 19 + 1
    d = n % 5 + 1
    (linestyle[a], linecolor[b], markershape[c], markercolor[d])
end

@userplot PKPlot
@userplot PKElimpPlot
@userplot PKElimpDrop
@userplot PdHLine

function luceil(x)
    fl = Int(floor(log10(x)))
    if fl < 0 fl = 0 end
    ceil(x/10^fl)*10^fl
end

@recipe function f(subj::PKPlot; lcd = :auto, tcd = :auto)
    x, y = subj.args

    if isa(lcd, Real)
        lc = luceil(maximum(x->isnan(x) ? -Inf : x, y)/lcd)
        yt = 0:lc:lc*lcd
    elseif isa(lcd, StepRange)
        yt = lcd
    elseif lcd == :all
        yt = y
    end

    if isa(tcd, Real)
        tc = luceil(maximum(x)/tcd)
        xt = 0:tc:tc*tcd
    elseif isa(tcd, StepRange)
        xt = tcd
    elseif tcd == :all
        xt = x
    end

    widen             --> true
    seriestype        --> :line
    xguide            --> "Time"
    link              --> :both
    legend            --> true
    grid              --> true
    gridstyle         --> :auto
    #ticks       := [nothing :auto nothing]
    #xlims             --> (minimum(x), maximum(x)*1.1)
    #ylims             --> (0, maximum(y)*1.1)

    if !isa(lcd, Symbol) || lcd != :auto
        yticks            --> yt
    end
    if !isa(tcd, Symbol) || tcd != :auto
        xticks            --> xt
    end

    seriescolor       --> :blue
    markershape       --> :circle
    markersize        --> 3
    markercolor       --> :match
    markerstrokealpha --> 0
    (x, y)
end

@recipe function f(subj::PKElimpPlot)
    x, y = subj.args
    seriestype        --> :line
    legend            --> false
    markersize        --> 0
    markerstrokealpha --> 0
    (x, y)
end

@recipe function f(subj::PKElimpDrop)
    x, y = subj.args
    seriestype        --> :scatter
    legend            --> false
    markersize        --> 4
    markercolor       --> :red
    markershape       --> :xcross
    (x, y)
end

@recipe function f(subj::PdHLine)
    x, y = subj.args
    seriestype        --> :straightline
    (x, [y, y])
end

# Text label from ID
function plotlabel(d, ld = nothing)
    title = ""
    if isnothing(d) return title end
    if length(d) > 0
        for (k, v) in d
            kv = k
            if !isnothing(ld) && haskey(ld, kv)
                kv = ld[kv]
            end
            title *= "$(kv) = $(v); "
        end
    end
    return title
end


"""
    pkplot(subj; ls = false, elim = false, xticksn = :auto, yticksn = :auto, kwargs...)

Plot for subject

* `ls` - concentration in log scale;
* `elim` - draw elimination curve;
* `xticksn` - number of ticks on x axis;
* `yticksn` - number of ticks on y axis;

*Other keywords:*

* `plotstyle` - predefined plot style from PKPLOTSTYLE;
* `drawbl` (`false`) - draw baseline, only for PDSubject;
* `drawth` (`false`) - draw threshold, only for PDSubject;

"""
function pkplot(subj::AbstractSubject; ls = false, elim = false, xticksn = :auto, yticksn = :auto, kwargs...)
    time = subj.time
    obs  = subj.obs
    kwargs = Dict{Symbol, Any}(kwargs)
    k = keys(kwargs)

    if !(:plotstyle in k)
        kwargs[:linestyle], kwargs[:linecolor], kwargs[:markershape],  kwargs[:markercolor]  = PKPLOTSTYLE[1]
    else
        kwargs[:linestyle], kwargs[:linecolor], kwargs[:markershape],  kwargs[:markercolor]  = kwargs[:plotstyle]
    end
    if !(:drawbl in k)
        kwargs[:drawbl] = false
    end
    if !(:drawth in k)
        kwargs[:drawth] = false
    end
    if !(:title in k)
        kwargs[:title] = plotlabel(subj.id)
    end
    if !(:legend in k)
        kwargs[:legend] = true
    end
    if !(:ylabel in k)
        kwargs[:ylabel] = "Concentration"
    end
    if !(:xlims in k)
        kwargs[:xlims] = (minimum(subj.time), maximum(subj.time)*1.1)
    end
    if :yscale in k
        if kwargs[:yscale] in [:ln, :log, :log2, :log10]
            ls = false
            if !(:minorticks in k) kwargs[:minorticks] = true end

            inds = findall(x-> x > 0, subj.obs)
            time = subj.time[inds]
            obs  = subj.obs[inds]
            if !(:yticks in k)
                if kwargs[:yscale] == :log10
                    b = 10
                elseif kwargs[:yscale] == :log2
                    b = 2
                elseif kwargs[:yscale] == :ln || kwargs[:yscale] == :log
                    b = ℯ
                end
                
                t = collect(floor(log(b, minimum(obs))):ceil(log(b, maximum(obs))))
                pushfirst!(t, first(t) - 1)
                kwargs[:yticks] = b .^ t
                
            end
            if !(:ylims in k)
                kwargs[:ylims] = (minimum(obs)*0.5, maximum(obs)*2.)
            else
                if kwargs[:ylims][1] <= 0
                    kwargs[:ylims] = (minimum(obs)/b, kwargs[:ylims][2])
                end
            end
        end
    else
        if !(:ylims in k)
            kwargs[:ylims] = (minimum(subj.obs), maximum(subj.obs)*1.15)
        end
    end

    if ls == true
        inds = findall(x-> x > 0, subj.obs)
        time = subj.time[inds]
        obs = log.(subj.obs[inds])
        if (:ylims in k)
            kwargs[:ylims] = (0, log(kwargs[:ylims][2]))
        end
    end

    p = pkplot(time, obs;  lcd = yticksn, tcd = xticksn, kwargs...)
    if elim
        if length(subj.keldata) > 0
            arsq, rsqn = findmax(subj.keldata.ar)
            lz        = subj.keldata.a[rsqn]
            lzint     = subj.keldata.b[rsqn]
            ts        = subj.keldata.s[rsqn]
            te        = subj.keldata.e[rsqn]
            if ls true
                x = [ts, te]
                y = [lzint + lz * x[1], lzint + lz * x[2]]
            else
                x = collect(ts:(te-ts)/100:te)
                y = exp.(lzint .+ lz .* x)
            end
            pkelimpplot!(p, x, y; title = kwargs[:title]*"\n($(round(lzint, sigdigits = 4)) + $(round(lz, sigdigits = 4)) * Time; aR² = $(round(arsq, sigdigits = 4))) ")
            if length(subj.kelrange.kelexcl) > 0
                pkelimpdrop!(p, time[subj.kelrange.kelexcl], obs[subj.kelrange.kelexcl])
            end
        end
    end
    if isa(subj, PDSubject)
        if kwargs[:drawth] == true
            pdhline!(p, [minimum(subj.time), maximum(subj.time)], getth(subj), lc = :blue, label = "TH")
        end
        if kwargs[:drawbl] == true
            pdhline!(p, [minimum(subj.time), maximum(subj.time)], getbl(subj), lc = :red, label = "BL")
        end
    end
    return p
end

function pkplot!(subj; ls = false, elim = false, xticksn = :auto, yticksn = :auto, kwargs...)
    time = subj.time
    obs = subj.obs
    kwargs = Dict{Symbol, Any}(kwargs)
    k = keys(kwargs)
    if !(:plotstyle in k)
        kwargs[:linestyle], kwargs[:linecolor], kwargs[:markershape],  kwargs[:markercolor]  = PKPLOTSTYLE[1]
    else
        kwargs[:linestyle], kwargs[:linecolor], kwargs[:markershape],  kwargs[:markercolor]  = kwargs[:plotstyle]
    end

    if !(:legend in k)
        kwargs[:legend] = true
    end
    if !(:xlims in k)
        kwargs[:xlims] = (minimum(subj.time), maximum(subj.time)*1.1)
    end
    if :yscale in k
        if kwargs[:yscale] in [:ln, :log, :log2, :log10]
            ls = false
            if !(:minorticks in k) kwargs[:minorticks] = true end
            inds = findall(x-> x > 0, subj.obs)
            time = subj.time[inds]
            obs  = subj.obs[inds]
            if !(:yticks in k)
                if kwargs[:yscale] == :log10
                    b = 10
                elseif kwargs[:yscale] == :log2
                    b = 2
                elseif kwargs[:yscale] == :ln || kwargs[:yscale] == :log
                    b = ℯ
                end
                t = collect(floor(log(b, minimum(obs))):ceil(log(b, maximum(obs))))
                pushfirst!(t, first(t) - 1)
                kwargs[:yticks] = b .^ t
            end
            if !(:ylims in k)
                kwargs[:ylims] = (minimum(obs)*0.5, maximum(obs)*2.)
            else
                if kwargs[:ylims][1] <= 0
                    kwargs[:ylims] = (minimum(obs)/b, kwargs[:ylims][2])
                end
            end
        end
    else
        if !(:ylims in k)
            kwargs[:ylims] = (minimum(subj.obs), maximum(subj.obs)*1.15)
        end
    end

    if ls == true
        inds = findall(x-> x > 0, subj.obs)
        time = subj.time[inds]
        obs = log.(subj.obs[inds])
        if (:ylims in k)
            kwargs[:ylims] = (0, log(kwargs[:ylims][2]))
        end
    end

    p = pkplot!(time, obs;  lcd = yticksn, tcd = xticksn, kwargs...)
    return p
end

function pageplot(data, id, ulist; kwargs...)
    kwargs = Dict{Symbol, Any}(kwargs)
    k = keys(kwargs)
    if !(:title in k)
        kwargs[:title] = plotlabel(id, kwargs[:ldict])
    end

    fst       = true
    p         = nothing
    labvec    = Vector{Int}(undef, 0)
    # Make subdata by ID
    isnothing(id) ? subdata   = data : subdata   = subset(data, id)
    # Y lims

    if !(:ylims in k) && length(subdata) > 1
        ysc   = :yscale in k
        ylmin = findmin(x->minconc(x, ysc), getdata(subdata))[1]
        ylmax = findmax(x->maxconc(x), getdata(subdata))[1]*1.15
        if ysc 
            ylmax *= 5
        end
        kwargs[:ylims] = (ylmin, ylmax)
    end
    # Plotting subdata
    if length(subdata) > 1 kwargs[:elim] = false end
    for subj in subdata
            if !isnothing(ulist)
                num = findfirst(x-> x ⊆ subj.id, ulist)
                if !isnothing(num)
                    style = plotstyle(num)
                    if num ∈ labvec
                        kwargs[:label] = nothing
                    else
                        kwargs[:label] = plotlabel(ulist[num], kwargs[:ldict])
                        push!(labvec, num)
                    end
                else
                    style = plotstyle(1)
                end
            else
                style = plotstyle(1)
            end
            if fst
                p = pkplot(subj; plotstyle = style, kwargs...)
                fst = false
            else
                pkplot!(subj; plotstyle = style, kwargs...)
            end
    end
    p
end

"""
    pkplot(data::DataSet{T};
    typesort::Union{Nothing, Symbol, AbstractVector{Symbol}} = nothing,
    pagesort::Union{Nothing, Symbol, AbstractVector{Symbol}, NoPageSort} = nothing,
    filter::Union{Nothing, Dict{Symbol}} = nothing,
    uylims::Bool = false,
    ldict = nothing,
    savepath::Union{Nothing, AbstractString} = nothing,
    namepref::Union{Nothing, AbstractString} = nothing,
    onlyplots = false,
    kwargs...) where T <: AbstractSubject

PK plot for subject set.

* `typesort` - sort on page by this id key;
* `pagesort` - different pages by this id key;
* `filter` - use only subjects if filter ⊆ subject id;
* `uylims` - same ylims for all dataset;
* `ldict` - Dict with labels for replace;
* `savepath` - path for plot saving;
* `namepref` - name prefix for saving files.
* `onlyplots` - if `true` return only vetor of plots;

Use `pagesort = MetidaNCA.NoPageSort()` to prevent page plotting (return single plot).

Return vector of pairs: `Page ID` => `Plot`.
"""
function pkplot(data::DataSet{T};
    typesort::Union{Nothing, Symbol, AbstractVector{Symbol}} = nothing,
    pagesort::Union{Nothing, Symbol, AbstractVector{Symbol}, NoPageSort} = nothing,
    filter::Union{Nothing, Dict{Symbol}} = nothing,
    uylims::Bool = false,
    ldict = nothing,
    savepath::Union{Nothing, AbstractString} = nothing,
    namepref::Union{Nothing, AbstractString} = nothing,
    onlyplots = false,
    kwargs...) where T <: AbstractSubject

    kwargs = Dict{Symbol, Any}(kwargs)
    k = keys(kwargs)
    if !(:ls in k)
        kwargs[:ls] = false
    end
    if !(:elim in k)
        kwargs[:elim] = false
    end
    if !(:drawbl in k)
        kwargs[:drawbl] = false
    end
    if !(:drawth in k)
        kwargs[:drawth] = false
    end
    if uylims && !(:ylims in k)
        kwargs[:ylims] = (findmin(x -> minconc(x), getdata(data))[1], findmax(x -> maxconc(x), getdata(data))[1]*1.15)
    end
    if !isnothing(filter) data = subset(data, filter) end

    if !isnothing(typesort)
        if isa(typesort, Symbol) typesort = [typesort] end
        typelist = uniqueidlist(data, typesort)
    else
        typelist = nothing
        if !(:legend in k)
            kwargs[:legend] = false
        end
    end
    p = []
    if isnothing(typesort) && isnothing(pagesort)
        printtitle = false
        if !(:title in k)
            printtitle = true
        end
        for subj in data
            if printtitle
                kwargs[:title] = plotlabel(subj.id, ldict)
            end
            if !(:legend in k)
                kwargs[:legend] = false
            end
            push!(p, subj.id => pkplot(subj; kwargs...))
        end
    elseif !isnothing(typesort) && isnothing(pagesort)
        printtitle = false
        if !(:title in k)
            printtitle = true
        end
        for subj in data
            if printtitle
                kwargs[:title] = plotlabel(subj.id, ldict)
            end
            if !(:legend in k)
                kwargs[:legend] = false
            end
            push!(p, subj.id =>  pageplot(data, subj.id, typelist; ldict, kwargs...))
        end
    elseif !isnothing(pagesort) && !isa(pagesort, NoPageSort)  
        if isa(pagesort, Symbol) pagesort = [pagesort] end
        pagelist = uniqueidlist(data, pagesort)
        for id in pagelist
            push!(p, id => pageplot(data, id, typelist; ldict, kwargs...))
        end
    else
        if !(:title in k) && !isnothing(filter)
            kwargs[:title] = plotlabel(filter)
        end
        push!(p, pageplot(data, nothing, typelist; ldict, kwargs...))
    end

    if !isnothing(savepath)
        if @isdefined savefig
            if isfile(savepath)
                error("File found on this path...")
            elseif !isdir(savepath)
                mkpath(savepath)
            end
            if isnothing(namepref) namepref = "plot" end
            for i = 1:length(p) 
                if isa(p[i], Pair)
                    savefig(p[i][2], joinpath(savepath, namepref*"_$(i).png"))
                else
                    savefig(p[i], joinpath(savepath, namepref*"_$(i).png"))
                end
            end
        else
            @warn "savefig not defined, install Plots.jl for plot writing... plots NOT saved..."
        end
    end
    if isa(pagesort, NoPageSort)
        return p[1]
    end

    if onlyplots return  getindex.(p, 2) end
    return p

end


"""
    pkplot(data::DataSet{T}; kwargs...) where T <: NCAResult
"""
function pkplot(data::DataSet{T}; kwargs...) where T <: NCAResult
    ds = map(x-> x.data, data)
    pkplot(ds; kwargs...)
end

"""
    pkplot(data::NCAResult; kwargs...) 
"""
function pkplot(data::NCAResult; kwargs...)
    pkplot(data.data; kwargs...)
end