---
author: "Vladimir Arnautov"
mainfont: "romanuni.ttf"
sansfont: "NotoSans-Regular.ttf"
monofont: "NotoSansMono-Regular.ttf"
title: "MetidaNCA validation report"
date: "2021-07-10"
mathfont: "texgyredejavu-math.otf"
---




# Introduction and package description

Non-compartment anlysis software.

See documentation: [https://pharmcat.github.io/MetidaNCA.jl/dev/](https://pharmcat.github.io/MetidaNCA.jl/dev/)

## Validation purpose

The main validation purpose is confirmation by examination and provision of objective evidence that software
specifications conform to user needs and intended uses, and that the particular requirements
implemented through software can be consistently fulfilled.

## Requirements

  * Julia 1.5.* (or higher) installed
  * Julia packages from dependence list installed (see [Project.toml](https://github.com/PharmCat/MetidaNCA.jl/blob/main/Project.toml))

## Developer software life cycle

  * Development stage
  * Testing procedures development
  * Performing testing procedures on local machine
  * Push to master branch
  * Performing testing procedures with GitHub Actions
  * Make pull request to the official registry of general Julia packages (if nessesary)
  * Make release (if previous completed)

### Versions

  * X.Y.Z - patch release (no breaking changes)
  * X.Y.# - minor release (may include breaking changes)
  * X.#.# - major release (breaking changes, changes in public API)
  * 0.#.# - no stable public API
  * ≥1.#.# - stable public API


## Build support

### Tier 1

  * julia-version: 1.5, 1.6
  * julia-arch: x64
  * os: ubuntu-18.04, macos-10.15, windows-2019

# Installation

## System information

 * Julia version: v"1.6.0-rc1"
 * Current machine: "x86_64-w64-mingw32"

## Installation method

MetidaNCA.jl can be installed by executing the following command in REPL.

```julia
import Pkg; Pkg.add("MetidaNCA")
```



## Version check

The installation process is checking within each testing job via GitHub Actions.
Also GitHub Action [chek](https://github.com/JuliaRegistries/General/blob/master/.github/workflows/automerge.yml)
performed before merging into JuliaRegistries/General repository
(see [Automatic merging of pull requests](https://github.com/JuliaRegistries/General#automatic-merging-of-pull-requests)).




Current package version:
"0.1.1"




# Operation qualification

This part of validation based on testing procedures entails running software products under known conditions with defined inputs and
documented outcomes that can be compared to their predefined expectations. All documented public API included in testing procedures and part of
critical internal methods.

## Coverage

Code coverage report available on [Codecov.io](https://app.codecov.io/gh/PharmCat/MetidaNCA.jl).
Test procedures include all  public API methods check.

* Coverage goal: ≥ 90.0%

## Data

Validation data available in the repository and included in the package. See Appendix 1.

## Testing results

```julia
Pkg.test("MetidaNCA")
```

```
Test Summary:                                               |
   Simple test                                              | No tests
Test Summary:                                               | Pass  Total
  Linear trapezoidal, Dose 100, Dosetime 0, no tau          |   18     18
Test Summary:                                               | Pass  Total
  Linear up Log down, Dose 120, Dosetime 0, tau 12          |   17     17
Test Summary:                                               | Pass  Total
  Log trapezoidal ATM, Dose 120, Dosetime 0, tau 12         |    5      5
Test Summary:                                               | Pass  Total
  Linear up Log down ATM, Dose 120, Dosetime 0, tau 12      |    4      4
Test Summary:                                               | Pass  Total
  Linear up Log down, Dose 100, Dosetime 0.25, tau 9        |    8      8
Test Summary:                                               | Pass  Total
  Linear trapezoidal, Dose 100, Dosetime 2.0, tau 10        |   15     15
Test Summary:                                               | Pass  Total
  Linear trapezoidal, Dose 100, Dosetime 0.0, tau 100       |   15     15
Test Summary:                                               | Pass  Total
  Linear trapezoidal, Dose 120, Dosetime 0.0, tau 12        |    4      4
Test Summary:                                               | Pass  Total
  set-get*! tests                                           |    7      7
Test Summary:                                               |
  applylimitrule!                                           | No tests
Test Summary:                                               |
  kel                                                       | No tests
```





# Performance qualification

Purpose of this testing procedures to demonstrate performance for some critical tasks.

#### Results




# Glossary

  * Installation qualification (IQ) - Establishing confidence that process equipment and ancillary systems are compliant with appropriate codes and approved design intentions, and that manufacturer's recommendations are suitably considered.
  * Operational qualification (OQ) Establishing confidence that process equipment and sub-systems are capable of consistently operating within established limits and tolerances.
  * Product performance qualification (PQ) - Establishing confidence through appropriate testing that the finished product produced by a specified process meets all release requirements for functionality and safety.
  * Repository - GitHub repository: https://github.com/PharmCat/MetidaNCA.jl
  * Master branch - main branch on GitHub ([link](https://github.com/PharmCat/MetidaNCA.jl/tree/main)).
  * Current machine - pc that used for validation report generating.

# Reference

* [General Principles of Software Validation; Final Guidance for Industry and FDA Staff](https://www.fda.gov/media/73141/download)
* [Guidance for Industry Process Validation: General Principles and Practices](https://www.fda.gov/files/drugs/published/Process-Validation--General-Principles-and-Practices.pdf)
* [Glossary of Computer System Software Development Terminology](https://www.fda.gov/inspections-compliance-enforcement-and-criminal-investigations/inspection-guides/glossary-computer-system-software-development-terminology-895)

# Appendix 1


```
.---------.-------------.---------.---------------.
| Subject | Formulation |    Time | Concentration |
|   Int64 |      String | Float64 |       Float64 |
:---------+-------------+---------+---------------:
|       1 |           T |     0.0 |           0.0 |
|       1 |           T |     0.5 |       178.949 |
|       1 |           T |     1.0 |       190.869 |
|       1 |           T |     1.5 |       164.927 |
|       1 |           T |     2.0 |       139.962 |
|       1 |           T |     2.5 |        129.59 |
|       1 |           T |     3.0 |       131.369 |
|       1 |           T |     4.0 |       150.854 |
|       1 |           T |     5.0 |       121.239 |
|       1 |           T |     6.0 |       139.229 |
|       1 |           T |     8.0 |        128.52 |
|       1 |           T |    10.0 |       143.243 |
|       1 |           T |    12.0 |       144.964 |
|       1 |           T |    24.0 |        133.16 |
|       1 |           T |    48.0 |       137.271 |
|       1 |           T |    72.0 |       112.846 |
|       2 |           R |     0.0 |           0.0 |
|       2 |           R |     0.5 |        62.222 |
|       2 |           R |     1.0 |       261.177 |
|       2 |           R |     1.5 |       234.063 |
|       2 |           R |     2.0 |       234.091 |
|       2 |           R |     2.5 |       222.881 |
|       2 |           R |     3.0 |       213.896 |
|       2 |           R |     4.0 |       196.026 |
|       2 |           R |     5.0 |       199.634 |
|       2 |           R |     6.0 |       196.037 |
|       2 |           R |     8.0 |       213.352 |
|       2 |           R |    10.0 |       200.088 |
|       2 |           R |    12.0 |       196.035 |
|       2 |           R |    24.0 |       160.338 |
|       2 |           R |    48.0 |        110.28 |
|       2 |           R |    72.0 |        85.241 |
|       3 |           R |     0.0 |           0.0 |
|       3 |           R |     0.5 |        49.849 |
|       3 |           R |     1.0 |        77.367 |
|       3 |           R |     1.5 |       105.345 |
|       3 |           R |     2.0 |       100.943 |
|       3 |           R |     2.5 |        72.746 |
|       3 |           R |     3.0 |        69.985 |
|       3 |           R |     4.0 |        93.565 |
|       3 |           R |     5.0 |        91.981 |
|       3 |           R |     6.0 |         82.71 |
|       3 |           R |     8.0 |        84.205 |
|       3 |           R |    10.0 |        85.342 |
|       3 |           R |    12.0 |        76.027 |
|       3 |           R |    24.0 |        81.259 |
|       3 |           R |    48.0 |        70.107 |
|       3 |           R |    72.0 |        67.901 |
|       4 |           R |     0.0 |           0.0 |
|       4 |           R |     0.5 |        52.421 |
|       4 |           R |     1.0 |       208.542 |
|       4 |           R |     1.5 |       188.923 |
|       4 |           R |     2.0 |       165.177 |
|       4 |           R |     2.5 |       146.996 |
|       4 |           R |     3.0 |       152.701 |
|       4 |           R |     4.0 |       154.345 |
|       4 |           R |     5.0 |       128.398 |
|       4 |           R |     6.0 |       149.807 |
|       4 |           R |     8.0 |       151.066 |
|       4 |           R |    10.0 |       136.819 |
|       4 |           R |    12.0 |       132.257 |
|       4 |           R |    24.0 |       141.247 |
|       4 |           R |    48.0 |       129.138 |
|       4 |           R |    72.0 |        97.625 |
|       5 |           T |     0.0 |           0.0 |
|       5 |           T |     0.5 |           0.0 |
|       5 |           T |     1.0 |         9.545 |
|       5 |           T |     1.5 |       153.964 |
|       5 |           T |     2.0 |        152.34 |
|       5 |           T |     2.5 |       151.452 |
|       5 |           T |     3.0 |       161.312 |
|       5 |           T |     4.0 |       169.334 |
|       5 |           T |     5.0 |       162.907 |
|       5 |           T |     6.0 |       166.651 |
|       5 |           T |     8.0 |       168.668 |
|       5 |           T |    10.0 |       155.103 |
|       5 |           T |    12.0 |       154.066 |
|       5 |           T |    24.0 |       162.974 |
|       5 |           T |    48.0 |       109.814 |
|       5 |           T |    72.0 |       110.778 |
|       6 |           T |     0.0 |           0.0 |
|       6 |           T |     0.5 |        57.882 |
|       6 |           T |     1.0 |       100.498 |
|       6 |           T |     1.5 |       138.651 |
|       6 |           T |     2.0 |       147.287 |
|       6 |           T |     2.5 |       154.648 |
|       6 |           T |     3.0 |       122.316 |
|       6 |           T |     4.0 |       132.857 |
|       6 |           T |     5.0 |       126.067 |
|       6 |           T |     6.0 |       140.466 |
|       6 |           T |     8.0 |       115.542 |
|       6 |           T |    10.0 |        102.16 |
|       6 |           T |    12.0 |       113.751 |
|       6 |           T |    24.0 |       101.049 |
|       6 |           T |    48.0 |         92.55 |
|       6 |           T |    72.0 |        69.501 |
|       7 |           R |     0.0 |           0.0 |
|       7 |           R |     0.5 |         19.95 |
|       7 |           R |     1.0 |       128.405 |
|       7 |           R |     1.5 |       136.807 |
|       7 |           R |     2.0 |       113.109 |
|       7 |           R |     2.5 |       153.254 |
|       7 |           R |     3.0 |       123.606 |
|       7 |           R |     4.0 |       142.655 |
|       7 |           R |     5.0 |       112.347 |
|       7 |           R |     6.0 |       139.919 |
|       7 |           R |     8.0 |       105.513 |
|       7 |           R |    10.0 |       134.408 |
|       7 |           R |    12.0 |        123.37 |
|       7 |           R |    24.0 |       110.511 |
|       7 |           R |    48.0 |        90.291 |
|       7 |           R |    72.0 |        58.051 |
|       8 |           R |     0.0 |           0.0 |
|       8 |           R |     0.5 |        136.91 |
|       8 |           R |     1.0 |       126.646 |
|       8 |           R |     1.5 |         118.5 |
|       8 |           R |     2.0 |       134.926 |
|       8 |           R |     2.5 |       113.213 |
|       8 |           R |     3.0 |       130.896 |
|       8 |           R |     4.0 |       138.327 |
|       8 |           R |     5.0 |        22.724 |
|       8 |           R |     6.0 |        53.774 |
|       8 |           R |     8.0 |        55.107 |
|       8 |           R |    10.0 |       102.871 |
|       8 |           R |    12.0 |       134.133 |
|       8 |           R |    24.0 |       108.021 |
|       8 |           R |    48.0 |        98.466 |
|       8 |           R |    72.0 |        74.437 |
|       9 |           T |     0.0 |           0.0 |
|       9 |           T |     0.5 |       113.362 |
|       9 |           T |     1.0 |       128.273 |
|       9 |           T |     1.5 |       125.395 |
|       9 |           T |     2.0 |       146.933 |
|       9 |           T |     2.5 |       140.559 |
|       9 |           T |     3.0 |       167.347 |
|       9 |           T |     4.0 |       157.504 |
|       9 |           T |     5.0 |        141.35 |
|       9 |           T |     6.0 |       140.282 |
|       9 |           T |     8.0 |       105.438 |
|       9 |           T |    10.0 |       164.843 |
|       9 |           T |    12.0 |        135.58 |
|       9 |           T |    24.0 |       117.125 |
|       9 |           T |    48.0 |       109.745 |
|       9 |           T |    72.0 |         93.44 |
|      10 |           R |     0.0 |           0.0 |
|      10 |           R |     0.5 |        13.634 |
|      10 |           R |     1.0 |        62.561 |
|      10 |           R |     1.5 |       112.655 |
|      10 |           R |     2.0 |       125.482 |
|      10 |           R |     2.5 |       116.255 |
|      10 |           R |     3.0 |       112.674 |
|      10 |           R |     4.0 |       116.986 |
|      10 |           R |     5.0 |        119.81 |
|      10 |           R |     6.0 |       107.557 |
|      10 |           R |     8.0 |       120.479 |
|      10 |           R |    10.0 |       124.171 |
|      10 |           R |    12.0 |       106.476 |
|      10 |           R |    24.0 |       116.508 |
|      10 |           R |    48.0 |        45.204 |
|      10 |           R |    72.0 |        42.191 |
'---------'-------------'---------'---------------'
```


