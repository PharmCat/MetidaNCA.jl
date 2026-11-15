# Details

## Data import

### Observation values

[`pkimport`](@ref) and [`pdimport`](@ref) convert every observation to one of three kinds:

* a number - a measured value; integer values are converted to float;
* `NaN` - a value below the limit of quantification (BLQ); `NaN` is the internal representation of BLQ;
* `missing` - no value: a missing cell, an empty string, any string that is not a number and is not the BLQ label.

The `blq` keyword sets the label of BLQ values in the source data. It can be `NaN` (default for `pkimport`), a string, or `nothing` (default for `pdimport`).

| Value in the data            | `blq = NaN`   | `blq = "BLQ"` | `blq = nothing` |
|:-----------------------------|:--------------|:--------------|:----------------|
| number                       | number        | number        | number          |
| `NaN`                        | `NaN` (BLQ)   | `missing`     | `missing`       |
| `missing`, `nothing`         | `missing`     | `missing`     | `missing`       |
| string `"BLQ"`               | `missing`     | `NaN` (BLQ)   | `missing`       |
| string `"NaN"`               | `NaN` (BLQ)   | `missing`     | `missing`       |
| numeric string, e.g. `"1.5"` | `1.5`         | `1.5`         | `1.5`           |
| any other string, `""` too   | `missing`     | `missing`     | `missing`       |

`blq = "NaN"` works as `blq = NaN`. If `blq` is not `NaN`, `NaN` values in the data are not BLQ and become `missing`.

The label is compared with the string in the data exactly: case and spaces matter, so `"BLQ "` and `"blq"` do not match `blq = "BLQ"` and become `missing`. Numeric strings are parsed with `tryparse(Float64, s)`, which accepts leading and trailing spaces.

```julia
ds = pkimport(data, :Time, :Concentration, :Subject; blq = "BLQ")
```

With `warn = true` a string column gives the warning `Some concentration values maybe not a number, try to fix.` and one warning `Value ... parsed as missing` for each string that becomes `missing`; an integer column gives `Concentration values transformed to float.`

[`upkimport`](@ref) does not take `blq`: concentration and volume columns of urine data must be numeric.

### Time values

Time values must be numbers and must be unique within a subject. `pkimport` raises an error for a subject with repeated time values and names the subject and the repeated values. Time values are sorted, observations follow their time values. The source table is not changed.

### BLQ and missing values in the analysis

Without a limit rule [`nca!`](@ref) treats BLQ (`NaN`) and `missing` values in the same way:

* values before the first measurable observation are skipped;
* values between measurable observations are interpolated (`intpm` keyword);
* values after the last measurable observation are dropped.

The BLQ label therefore changes the result only through a [`LimitRule`](@ref), see the next section.

## Using LimitRule

```julia
ll = LimitRule(;lloq = 0.1, btmax = 0.0, atmax = NaN, nan = NaN, rm = true)

```

It means that all numeric values at or below `lloq` will be replaced by `btmax` before Tmax and replaced by `atmax` after Tmax; BLQ (`NaN`) and `missing` values will be replaced by `nan` (if `nan` is not `NaN`). If `rm` is `true`, all BLQ (`NaN`) and `missing` values will be deleted.

The steps are applied in this order: `nan`, then `lloq`, then `rm`. The `lloq` step does not change BLQ (`NaN`) values. To apply `btmax` and `atmax` to BLQ values, set `nan` to a value at or below `lloq`: BLQ values are replaced by `nan` in the first step and then by `btmax` or `atmax` in the second. `missing` values are changed in the same way, because the `nan` step does not distinguish BLQ from `missing`.

```julia
# BLQ before Tmax -> 0, BLQ after Tmax -> NaN
ll = LimitRule(; nan = 0.0, lloq = 0.1, btmax = 0.0, atmax = NaN)
```

See also: [`applylimitrule!`](@ref).

## Using DoseTime

```julia

dt = DoseTime(dose = 200.0, time = 0.0)

```

DoseTime can be appliet to each subject or dataset and can be used with [`pkimport`](@ref).  


```julia

ds = pkimport(pkdata2, :Time, :Concentration, [:Subject, :Formulation]; dosetime = dt)

```

DoseTime for staedy-state PK:

```julia
dt = DoseTime(dose = 100.0, time = 0.25, tau = 9.0)
```

See also: [`setdosetime!`](@ref).


## Calculation steps for PK NCA

### Step 1

Filter all values before dose time and `NaN` or `missing` values after last measurable concentration.
If TAU set, calculate start and end timepoints for AUCtau.

!!! note
    By default last `DoseTime` used to get dose time. For using another dose time - set `usedose` keyword as number of this dose.

### Step 2

Cmax, Tmax calculation. Interpolate `NaN` and `missing` values.

!!! note
    If more than one maximum - only first observation used for define Tmax.

### Step 3

Exclude interpolated points from calculation (add to `excltime`). Elimination parameters calculation. Find last concentration > 0 and time for last concentration > 0.

!!! note
    If `kelstart` or `kelend` in `excltime` then `kelauto` set to `true`.


!!! note
    If `kelauto` is `true` than range of observations for elimination will start from Tmax if administration set as `iv`, and from next observation after Tmax in other cases.

Adjusted R square user to choose best fit.

```math
R^2_{adj} = 1 \frac{(1 - R^2) \cdot (n - 1)}{(n - 2)}
```

### Step 4

Shift all time values by dose time.

### Step 5

Calculate dose concentration (Cdose).

!!! note
    If there is no concentration for dosing time:
    * If administration set as `iv` if 1st observation > than 2nd and both > 0 - Dose concentration is log-extrapolated, else set as 1st observation.
    * If administration not `iv`, than if Tau used  Dose concentration set as minimal concentration, in other case set as 0.  

### Step 6

Calculate areas.

!!! note
    If AUClast is 0, than AUClast, AUMClast and AUCall set as `NaN`, so other dependent parameters is `NaN` too.   

### Step 7

Calculate steady-state parameters.

!!! note
    If end of tau interval lies between two observation, than interpolation used to compute Ctau and partial AUCs; `intpm` keyword used to define interpolation method.
    If end of tau interval lies after all observation, than extrapolation used to compute Ctau and partial AUCs. Extrapolation based on using elimination parameters.


## [Unitful details](@id unitful_details)

!!! warning
  **Unitful.jl** MetidaNCA can work with [Unitful.jl](https://painterqubits.github.io/Unitful.jl/stable/). There is no guarantee that all functions will work without errors. All validation procedures with Unitful should be done manually before use.


!!! warning
  **Dose and time settings** If you are using Unitful, check `dositime` settings: `DoseTime(dose = 100u"mg", time = 0u"hr")`. For properly results all values should have units (including time and concentration data in data table).
