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

### Plot 1

![](plot1.png)

### Plot 2

![](plot2.png)

### Plot 3

![](plot3.png)

### Plot 4

![](plot4.png)
