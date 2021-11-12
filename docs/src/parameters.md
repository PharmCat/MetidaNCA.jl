# [Parameter list](@id parameter_list)

## Basic parameters

### :Cmax

Maximum concentration from dose time to dose time + tau (if tau > 0). Firs observation used.

### :Tmax

Time at maximum concentration from dose time to dose time + tau (if tau > 0). Firs observation used.

### :Cdose

By default dose time is 0. If concentration at dose time present in observation list - this concentration will be used.
For extravascular setting (:ev) if τ used (τ > 0) Cdose set as minimum concentration from dose time to τ time [:Ctaumin](#:Ctaumin), else set equal to zero.
For IV (:iv) if 1-st observation > 2-nd observation > 0 then logarithmic extrapolation used, else set equal to 1-st observation.

### AUC / AUMC

Area under Curve / Area under the Moment Curve.

```math
AUC = \sum_{n=1}^N AUC_{n}
```

```math
AUMC = \sum_{n=1}^N AUMC_{n}
```

Where `AUCn`/`AUMCn`- partial AUC/AUMC.

#### Linear trapezoidal rule

```math
AUC\mid_{t_1}^{t_2} = \delta t \times \frac{C_1 + C_2}{2}
```

```math
AUMC\mid_{t_1}^{t_2} = \delta t \times \frac{t_1 \times C_1 + t_2 \times C_2}{2}
```

#### Logarithmic trapezoidal rule

```math
AUC\mid_{t_1}^{t_2} =   \delta t \times \frac{ C_2 - C_1}{ln(C_2/C_1)}
```

```math
AUMC\mid_{t_1}^{t_2} = \delta t \times \frac{t_2 \times C_2 - t_1 \times C_1}{ln(C_2/C_1)} -  \delta t^2 \times \frac{ C_2 - C_1}{ln(C_2/C_1)^2}
```

#### Interpolation

##### Linear interpolation rule

```math
C_x = C_1 + \frac{(t_x-t_1)\times(C_2 - C_1)}{t_2 - t_1}
```

##### Logarithmic interpolation rule

```math
C_x = exp\left(ln(C_1) + \frac{(t_x-t_1)\times(ln(C_2) - ln(C_1))}{t_2 - t_1}\right)
```

#### :AUClast

Area under the curve from dose time to last observed concentration (>0).

#### :AUMClast

Area under the Moment Curve from dose time to last observed concentration (>0).
Dose time is the starting point for this calculation.

#### :AUCall

All values used to calculate AUC.

### :Kel

𝝺z - elimination constant. Linear regression at the terminal phase used for logarithmic transformed concentration data.

### :HL

Half-Life; T1/2

```math
HL = ln(2) / \lambda_z
```

### :Rsq

 Coefficient of determination (R²).

### :ARsq

Adjusted coefficient of determination (R²).

### :MRTlast

Mean residence time (MRT) from the dose time to the time of the last observed concentration.

```math
MRT_{last} = AUMC_{last} / AUC_{last}
```

## If :Kel calculated

### :AUCinf

AUC extrapolated from the last observed concentration to infinity.

```math
AUC_\infty = AUC_{last} + \frac{C_{last}}{\lambda_z}
```

### :AUMCinf

AUMC extrapolated from the last observed concentration to infinity.

```math
AUMC_\infty =  AUMC_{last} + \frac{t_{last}\times C_{last}}{\lambda_z} + \frac{C_{last}}{\lambda_z^2}
```

### :AUCpct

Percentage of AUCinf due to extrapolation from the last observed concentration to infinity.

```math
AUCpct = (AUC_\infty - AUC_{last}) / AUC_\infty * 100 \%
```

#### :AUCinf_pred

AUC extrapolated to infinity from the predicted concentration.

```math
AUC_{\infty pred} = AUC_{last} + \frac{C_{last pred}}{\lambda_z}
```

## If Dose used

### Clearance

#### :Cllast

```math
CL_{last} = Dose / AUC_{last}
```

#### :Clinf

Total body clearance for extravascular administration.

```math
CL_\infty = Dose / AUC_\infty
```

#### :Vzinf

Volume of distribution based on the terminal phase.

##  Steady-state parameters (If τ used)

τ-time = dose_time + τ

### :AUCtau

Area under the curve from dose time to τ-time.

### :AUMCtau

Area under the Moment Curve from the dose time to τ-time.

### :Ctau

Concentration at τ-time.

### :Ctaumin

Minimum concentration from the dose time to τ-time.

### :Cavg

```math
C_{avg} = AUC_\tau / \tau
```

### :Fluc

Fluctuation

```math
Fluc = ( C_{max} - C_{\tau min} ) / C_{avg} * 100 \%
```

### :Fluctau

Fluctuation Tau

```math
Fluc\tau = ( C_{max} - C_{\tau} ) / C_{avg} * 100 \%
```

### :Accind

Accumulation index.

```math
Accind = \frac{1}{1 - exp(-\lambda_z \tau)}
```

### :MRTtauinf

```math
MRT_{\tau\inf} = (AUMC_\tau + \tau * (AUC_\infty - AUC_\tau)) / AUC_\tau
```

### :Swing

```math
Swing = (C_{max} - C_{\tau min}) / C_{\tau min}
```

### :Swingtau

```math
Swing_{\tau} = (C_{max} - C_{\tau}) / C_{\tau}
```
