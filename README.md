# MetidaNCA

This program comes with absolutely no warranty. No liability is accepted for any loss and risk to public health resulting from use of this software.

| Status | Cover | Build | Docs |
|--------|-------|-------|------|
|[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)|[![codecov](https://codecov.io/gh/PharmCat/MetidaNCA.jl/branch/main/graph/badge.svg?token=A9eyT9g0WZ)](https://codecov.io/gh/PharmCat/MetidaNCA.jl)|![Tier 1](https://github.com/PharmCat/MetidaNCA.jl/workflows/Tier%201/badge.svg) | [![Latest docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://pharmcat.github.io/MetidaNCA.jl/dev/) [![Stable docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://pharmcat.github.io/MetidaNCA.jl/stable/)|


Non-compartment PK analysis (NCA).

Pharmacokinetics, sometimes abbreviated as PK, is a branch of pharmacology dedicated to determine the fate of substances administered to a living organism.

When analyzing pharmacokinetic data, one generally employs either model fitting using nonlinear regression analysis or non-compartmental analysis techniques (NCA). The method one actually employs depends on what is required from the analysis. If the primary requirement is to determine the degree of exposure following administration of a drug (such as AUC), and perhaps the drug's associated pharmacokinetic parameters, such as clearance, elimination half-life, T (max), C (max), etc., then NCA is generally the preferred methodology to use in that it requires fewer assumptions than model-based approaches.[*]

PK urine parameters and PD parameters such as Time Above/Below  Baseline/Threshold can be also calculated.

Also this package include recipes for plotting PK/PD data.

* Gabrielsson J, Weiner D. Non-compartmental analysis. Methods Mol Biol. 2012;929:377-89. doi: 10.1007/978-1-62703-050-2_16. PMID: 23007438.

## Installation

```julia
import Pkg; Pkg.add("MetidaNCA")
```

## Test

```julia
Pkg.test("MetidaNCA")
```

## First step

```julia
using DataFrames, CSV, MetidaNCA

pkdata2  = CSV.File(joinpath(dirname(pathof(MetidaNCA)), "..", "test", "csv", "pkdata2.csv")) |> DataFrame

ds = MetidaNCA.pkimport(pkdata2, :Time, :Concentration, [:Subject, :Formulation]; dosetime = MetidaNCA.DoseTime(dose = 100, time = 0))

sort!(ds, :Subject)

dsnca = MetidaNCA.nca!(ds, adm = :ev, calcm = :lint, verbose = 2)
```
