# Parameters

## Basic parameters

### Cmax

Maximum concentration from dose time to dose time + tau (if tau > 0). Firs observation used.

### Tmax

Time at maximum concentration from dose time to dose time + tau (if tau > 0). Firs observation used.

### Cdose

Concentration at dose time.

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

#### AUClast / AUMClast

Area from dose time to last observed concentration (>0).

#### AUCall / AUMCall

All values used to calculate AUC/AUMC.


### ```math \lambda_z ```- elimination constant

Linear regression used for logarithmic transformed concentration data.

### Half-Life; T1/2

```math
HL = ln(2) / \lambda_z
```

## If Kel calculated

### AUCinf

```math
AUCinf = AUClast + \frac{Clast}{\lambda_z}
```

### AUMCinf

```math
AUMCinf =  AUMClast + \frac{tlast\times Clast}{\lambda_z} + \frac{Clast}{\lambda_z^2}
```

### AUCpct

```math
AUCpct = (AUCinf - AUClast) / AUCinf * 100.0%
```

## If Dose used

### Clearance

#### Cllast

```math
Cllast = Dose / AUClast
```

#### Clinf

```math
Clinf = Dose / AUCinf
```

## If Tau used

### AUCtau / AUMCtau

Area from dose time to dose time + tau.

### Accumulation index

```math
Accind = \frac{1}{1 - exp(-\lambda_z \tau)}
```

### MRTtauinf

```math
MRTtauinf = (AUMCtau + \tau * (AUCinf - AUCtau)) / AUCtau
```
