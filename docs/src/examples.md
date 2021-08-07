# Examples

```@setup ncaexample
ENV["GKSwstype"] = "nul"
```

## Import

```@example ncaexample
using MetidaNCA, CSV, DataFrames;

pkdata2 = CSV.File(joinpath(dirname(pathof(MetidaNCA)), "..", "test", "csv",  "pkdata2.csv")) |> DataFrame

ds = MetidaNCA.pkimport(pkdata2, :Time, :Concentration, [:Subject, :Formulation]; dosetime = MetidaNCA.DoseTime(dose = 100, time = 0))

sort!(ds, :Subject)
```

## NCA

```@example ncaexample
dsnca = MetidaNCA.nca!(ds, adm = :ev, calcm = :lint)

dsnca[:, :AUClast]
```

## Plotting

```@example ncaexample
using Plots

p = MetidaNCA.pkplot(ds; typesort = :Subject, pagesort = nothing, sort = Dict(:Formulation => "R"))

png(p, "plot1.png")

p = MetidaNCA.pkplot(ds; typesort = :Formulation, pagesort = nothing, legend = true)

png(p, "plot2.png")

p = MetidaNCA.pkplot(ds; elim = true, ls = true)

png(p[1], "plot3.png")

p = MetidaNCA.pkplot(ds; typesort = :Subject, pagesort = :Formulation)

png(p[1], "plot4.png")
```

#### Plot 1

![](plot1.png)

#### Plot 2

![](plot2.png)

#### Plot 3

![](plot3.png)

#### Plot 4

![](plot4.png)

## Set dose time

```@example ncaexample
dt = MetidaNCA.DoseTime(dose = 200, time = 0)

MetidaNCA.setdosetime!(ds, dt, Dict(:Formulation => "R"))

dsnca = MetidaNCA.nca!(ds)

dsnca[:, :Dose]
```

## Set range for elimination

```@example ncaexample
kr =  MetidaNCA.ElimRange(kelstart = 4, kelend = 12, kelexcl = Int[5,6])

MetidaNCA.setkelrange!(ds, kr, [1,2,3])

dsnca = MetidaNCA.nca!(ds)

p = MetidaNCA.pkplot(ds[1]; elim = true)

png(p, "plot5.png")

getkeldata(ds[1])
```

#### Plot 5

![](plot5.png)


## Without import

```@example ncaexample

dsnca = MetidaNCA.nca(pkdata2, :Time, :Concentration, [:Subject, :Formulation]; dosetime = MetidaNCA.DoseTime(dose = 100, time = 0))

sort!(dsnca, :Subject)

dsnca[:, :AUClast]
```
