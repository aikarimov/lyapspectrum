# lyapspectrum

`lyapspectrum` Calculate the Lyapunov spectrum for a particular system

```
[L, LSPAN, LEXP]  = lyapspectrum(ODEFUN,TSPAN,Y0)
[L, LSPAN, LEXP]  = lyapspectrum(ODEFUN,TSPAN,Y0,'NAME',VALUE)
```

### Input:
`ODEFUN` is a function handle `ODEFUN(T,Y)`

`TSPAN = [T0 T1 ... TFINAL]` is a vector of times

`Y0` is a vector of initial conditions

`'NAME',VALUE` gives an additional input options:

`'jacobian',JFCN` is a Jacobian matrix ODEFUN(T,Y)

> Example: `L  = lyapspectrum(ODEFUN,TSPAN,Y0,'jacobian',JFCN)`

`'disp','none'` - sets disply style, 'none' shows no plot (default),

`'disp','2d'` - shows 2d plots time vs values of Lyapunov exponents,

`'disp','3d'` - shows 3d plots with colors corresponding to local values of Lyapunov exponents,

`'disp','all'` - shows all plots

> Example: `L  = lyapspectrum(ODEFUN,TSPAN,Y0,'disp','3d')`

`'df',N` - divides each step by N to calculate the local Lyapunov exponents (default N = 30) 

> Example: `L  = lyapspectrum(ODEFUN,TSPAN,Y0,'df',10)`

`'trans',TTRANS` - skips TTRANS time before calculating the spectrum

> Example: `L  = lyapspectrum(ODEFUN,TSPAN,Y0,'trans',30)`

### Output:
`L` - vector of averaged Lyapunov exponents (base e),

`LSPAN` - matrix of local Lyapunov exponents evolution over times `TSPAN`

`LEXP` - matrix of global Lyapunov exponents evolution over times `TSPAN`
 
Copyright (C) 2023, Karimov A.I.

## Illustrative example

```
[L,~,Lexp] = lyapspectrum(@lorenz2,tspan,y0,'disp','all','view',[0 1 0],'jacobian',@Jlorenz2,'trans',10);
```

shows two pictures.

![This is an image 1](https://github.com/aikarimov/lyapspectrum/blob/main/Lorenz_all_lyap.jpg)
                      
![This is an image 2](https://github.com/aikarimov/lyapspectrum/blob/main/Lorenz_all_series.jpg)

