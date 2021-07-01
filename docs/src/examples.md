# Examples

```@example ncaexample
using MetidaNCA, CSV, DataFrames;

pkdata2 = CSV.File(joinpath(dirname(pathof(MetidaNCA)), "..", "test", "csv",  "pkdata2.csv")) |> DataFrame

ds = MetidaNCA.pkimport(pkdata2, :Time, :Concentration, [:Subject, :Formulation]; dosetime = MetidaNCA.DoseTime(dose = 100, time = 0))

sort!(ds, :Subject)

dsnca = MetidaNCA.nca!(ds, adm = :ev, calcm = :lint)

dsnca[:, :AUClast]
```
