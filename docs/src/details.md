# Details

## Using LimitRule

```julia
ll = LimitRule(;lloq = 0.1, btmax = 0.0, atmax = NaN, nan = NaN, rm = true)

```

It means that all values below `lloq` will be replaced by `btmax` before Tmax and replaced by `atmax` after Tmax; `NaN` values will be replaced by `nan`. If `rm` is `true`, all `NaN` values will be deleted. 


## Calculation steps for PK NCA

### Step 1

Filter all values before dose time and `NaN` or `missing` values after last measurable concentration.
If TAU set, calculate start and end timepoints for AUCtau.

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
  **Unitful.jl**
  MetidaNCA can work with [Unitful.jl](https://painterqubits.github.io/Unitful.jl/stable/).
  There is no guarantee that all functions will work without errors.
  All validation procedures with Unitful should be done manually before use.


!!! warning
  **Dose and time settings**
  If you are using Unitful, check `dositime` settings: `DoseTime(dose = 100u"mg", time = 0u"hr")`.
  For properly results all values should have units (including time and concentration data in data table).
