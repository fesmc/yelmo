# Thermodynamics

## Equations

### Ice column

Prognostic equation:

$$
\frac{\partial T}{\partial t} = \frac{k}{\rho c} \frac{\partial^2 T}{\partial z^2} - u \frac{\partial T}{\partial x} - v \frac{\partial T}{\partial y} - w \frac{\partial T}{\partial z} + \frac{\Phi}{\rho c}
$$

Ice surface boundary condition:

$$
T(z=z_{\rm srf}) = {\rm min}(T_{\rm 2m},T_0)
$$

Ice base (temperate) boundary condition:

$$
T(z=z_{\rm bed}) = T_{\rm pmp}
$$

Ice base (frozen) boundary condition:

$$
k \frac{\partial T}{\partial z} = k_r \frac{\partial T_r}{\partial z}
$$

Note, the following internal Yelmo variables are defined for convenience:

$$
Q_{\rm ice,b} = -k \frac{\partial T}{\partial z}; \quad
Q_{\rm rock} = -k_r \frac{\partial T_r}{\partial z}
$$

### Bedrock column

Prognostic equation:

$$
\frac{\partial T_r}{\partial t} = \frac{k_r}{\rho_r c_r} \frac{\partial^2 T_r}{\partial z^2}
$$

Bedrock surface boundary condition:

$$
T_r(z=z_{\rm bed}) = T(z=z_{\rm bed})
$$

Bedrock base boundary condition:

$$
\frac{\partial T_r}{\partial z} = -\frac{Q_{\rm geo}}{k_r}
$$

## Equilibrium bedrock

In this case, the bedrock temperature profile is prescribed to the equilibrium linear temperature profile. The slope follows:

$$
\frac{\partial T_r}{\partial z} = -\frac{Q_{\rm geo}}{k_r}
$$

and the bedrock surface temperature is given by the ice temperature at its base:

$$
T_r(z=z_{\rm bed}) = T(z=z_{\rm bed})
$$

## Active bedrock

Yelmo calculates the temperature in the lithosphere along with the ice temperature. This can be achieved by assuming equilibrium conditions in the bedrock, i.e., that the temperature profile in the bedrock is always linear with `T_lith_s = T_ice_b` and the slope equal to `dT/dz = -Q_geo / k_lith`. Or, the temperature equation can be solved in the lithosphere together with the temperature in the ice column.

The parameter block `ytherm_lith` controls how the lithosphere is calculated with `ytherm_lith.method=['equil','active']` deciding the two cases above.

### Density of the upper lithosphere

### Heat capacity of the upper lithosphere

In both SICOPOLIS and GRISLI, a value of `cp = 1000.0 [J kg-1 K-1]` is used (referenced in Rogozhina et al., 2012; Greve, 2005; Greve, 1997). This value is adopted in Yelmo as well.

```fortran
cp = 1000.0    ! [J kg-1 K-1] 
```

### Heat conductivity of the upper lithosphere

Note, Yelmo expects input parameter values in units of `[J a-1 m-1 K-1]`, while much literature uses `[W m-1 K-1]`. Given the number of seconds in a year `sec_year = 31536000.0`, `kt [W m-1 K-1] * sec_year = kt [J a-1 m-1 K-1]`.

Rogozhina et al. (2012) use `kt = 2 [W m-1 K-1]` for Greenland:

```fortran
kt = 6.3e7     ! [J a-1 m-1 K-1]
```

This value is supported by Lösing et al. (2020), who perform a Bayesian inversion for GHF in Antarctica. Assuming exponentially decreasing heat production with depth, lower values of `kt` are supported (see Fig. 7b). In a study on the global thermal characteristics of the lithosphere, Cammarano and Guerri (2017) adopt an upper crust thermal conductivity of `kt = 2.5 [W m-1 K-1]`.

To do: This study is also potentially relevant:
[https://link.springer.com/article/10.1186/s40517-020-0159-y](https://link.springer.com/article/10.1186/s40517-020-0159-y). They show ranges of on the order of `kt = 2-3 [W m-1 K-1]` for the Canadian shield.

The above value of `kt = 2 [W m-1 K-1] = 6.3e7 [J a-1 m-1 K-1]` is adopted as the default thermal conductivity of the upper crust in Yelmo. For historical context, see other estimates below.

From Greve (1997) and Greve (2005):

```fortran
kt = 9.46e7    ! [J a-1 m-1 K-1]
```

which is equivalent to `kt = 3 [W m-1 K-1]`. The source of this value is not known.

From GRISLI:

```fortran
kt = 1.04e8    ! [J a-1 m-1 K-1]
```

which is equivalent to `kt = 3.3 [W m-1 K-1]`. The source of this value is not known.
