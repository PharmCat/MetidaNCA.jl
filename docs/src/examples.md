# Examples

```@setup ncaexample
ENV["GKSwstype"] = "nul"
```

## Import

Use [`pkimport`](@ref) to import PK data from table to subject set.

```@example ncaexample
using MetidaNCA, CSV, DataFrames;

pkdata2 = CSV.File(joinpath(dirname(pathof(MetidaNCA)), "..", "test", "csv",  "pkdata2.csv")) |> DataFrame

ds = pkimport(pkdata2, :Time, :Concentration, [:Subject, :Formulation]; dosetime = DoseTime(dose = 100, time = 0))

sort!(ds, :Subject)
```

## NCA

Perform NCA analysis with [`nca!`](@ref). Access to result set is similar to DataFrame or any table.
Find parameter list [here](@ref parameter_list).

```@example ncaexample
dsnca = nca!(ds, adm = :ev, calcm = :lint)

dsnca[:, :AUClast]
```

## Partial AUC

```@example ncaexample
dsnca = nca!(ds, adm = :ev, calcm = :lint, partials = [(1, 7)])

dsnca[:, :AUC_1_7]
```

## Result modification or custom parameter function 

```@example ncaexample

# Define modify! function for new parameter
function newparam(data)
    data.result[:AUChalf] = data.result[:AUClast] / 2
end

dsncanp = nca!(ds, modify!  = newparam)

dsncanp[1][:AUChalf]
```

Function `newparam` applyed to [`NCAResult`](@ref).


## Print output

```@example ncaexample
dsnca = nca!(ds[1], adm = :ev, calcm = :lint, verbose = 2);

```

## Plotting

```@example ncaexample
using Plots

# plot 1

p = pkplot(ds; typesort = :Subject, pagesort = NoPageSort(), filter = Dict(:Formulation => "R"))

png(p, "plot1.png")

# plot 2

p = pkplot(ds; typesort = :Formulation, pagesort = NoPageSort(), legend = true)

png(p, "plot2.png")

# plot 3

p = pkplot(ds; elim = true, ls = true)

png(p[1][2], "plot3.png")

# plot 4
# If pagesort used - return pairs with `Page ID` => `Plot`

p = pkplot(ds; typesort = :Subject, pagesort = :Formulation)

png(p[1][2], "plot4.png")

# plot 5

p = vpcplot(ds)

png(p, "vpcplot.png")

```

#### Plot 1

![](plot1.png)

#### Plot 2

![](plot2.png)

#### Plot 3

![](plot3.png)

#### Plot 4

![](plot4.png)

#### Plot 5

![](vpcplot.png)

## Set dose time

You can set dose time with [`setdosetime!`](@ref) for whole subject set or for
selected subjects.

```@example ncaexample
dt = DoseTime(dose = 200, time = 0)

setdosetime!(ds, dt, Dict(:Formulation => "R"))

dsnca = nca!(ds)

dsnca[:, :Dose]
```

## Set range for elimination

By default no exclusion or range specified. With [`setkelrange!`](@ref) elimination range and exclusion
can be specified for whole subject set or for any selected subjects.

```@example ncaexample
kr =  ElimRange(kelstart = 4, kelend = 12, kelexcl = Int[5,6])

setkelrange!(ds, kr, [1,2,3])

dsnca = nca!(ds)

p = pkplot(ds[1]; elim = true)

png(p, "plot5.png")

getkeldata(ds[1])
```

#### Plot 5

![](plot5.png)


## Without import

You  can use [`nca`](@ref) for NCA analysis directly from tabular data.

```@example ncaexample

dsnca = nca(pkdata2, :Time, :Concentration, [:Subject, :Formulation]; dosetime = DoseTime(dose = 100, time = 0))

sort!(dsnca, :Subject)

dsnca[:, :AUClast]
```

## PD subject

Use [`pdimport`](@ref) to import PD data from table to subject set.

#### Import & NCA

```@example ncaexample

pddata = CSV.File(joinpath(dirname(pathof(MetidaNCA)), "..", "test", "csv",  "pddata.csv")) |> DataFrame

pd =  MetidaNCA.pdimport(pddata, :time, :obs, :subj; bl = 1.5, th = 5.0)

MetidaNCA.nca!(pd[1])
```

#### PD subject plotting

```@example ncaexample

p = MetidaNCA.pkplot(pd[1], drawth = true, drawbl = true)

png(p, "plot6.png")

# Plot DoseTime (can be used for PK plots)

dt = DoseTime(dose = 200, time = 1.5)
setdosetime!(pd, dt)

p = MetidaNCA.pkplot(pd[1], drawdt = true)


png(p, "plot7.png")
```

##### Plot 6

![](plot6.png)

##### Plot 7

![](plot7.png)
