# Pellet heating model for fast firing

Python GUI tool that simulates the transient heating of a cylindrical ceramic pellet during fast firing. It solves 2D axisymmetric heat conduction with convective and radiative exchange with the furnace and heat transfer through an alumina support, and reports the temperature evolution and the times needed to reach 50 % and 99 % of the furnace temperature at the pellet centre.

> Developed at Universitat Jaume I (UJI). If you use this code, please cite it (see [Citation](#citation)).

## Contents

- [Model](#model)
- [Requirements](#requirements)
- [Usage](#usage)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [Fixed internal parameters](#fixed-internal-parameters)
- [Limitations](#limitations)
- [Citation](#citation)
- [License](#license)

## Model

The pellet is a cylinder of radius $R = D/2$ and thickness $H$. Temperature $T(r, z, t)$ obeys the axisymmetric heat equation with constant diffusivity $\alpha$:

$$\frac{\partial T}{\partial t} = \alpha \left( \frac{\partial^2 T}{\partial r^2} + \frac{1}{r}\frac{\partial T}{\partial r} + \frac{\partial^2 T}{\partial z^2} \right)$$

**Boundary conditions**

- **Axis ($r = 0$):** symmetry, $\partial T / \partial r = 0$.
- **Lateral surface ($r = R$) and top ($z = H$):** convection plus linearised radiation,

$$-\kappa \frac{\partial T}{\partial n} = h\,(T - T_{gas}) + h_{rad}\,(T - T_{rad}), \qquad h_{rad} = \varepsilon \sigma (T + T_{rad})(T^2 + T_{rad}^2)$$

- **Base ($z = 0$):** heat exchange with the support at $T_{base}$ through a contact resistance $R_c$ in series with an alumina layer of thickness $\delta$,

$$-\kappa \frac{\partial T}{\partial n} = h_{base}\,(T - T_{base}), \qquad h_{base} = \left( R_c + \frac{\delta}{\kappa_{Al_2O_3}(T)} \right)^{-1}$$

  where $\kappa_{Al_2O_3}(T) = 85.686\,(T/273.15)^{-1.01}$ W/m/K (Touloukian-type correlation, $T$ in K, lower limit 2 W/m/K).

**Furnace environment.** The furnace is assumed to be already stabilised at $T_{furnace}$. The gas and the radiating walls seen by the pellet approach this temperature with short first-order transients:

$$T_{gas}(t) = T_0 + (T_{furnace} - T_0)\left(1 - e^{-t/\tau_{gas}}\right), \qquad T_{rad}(t) = T_0 + (T_{furnace} - T_0)\left(1 - e^{-t/\tau_{rad}}\right)$$

**Numerics.** Finite differences on a uniform $15 \times 15$ ($r \times z$) grid (method of lines), integrated in time with `scipy.integrate.solve_ivp` (BDF, `rtol = 1e-6`, `atol = 1e-7`).

**Heating times.** $t_{50}$ and $t_{99}$ are the first times at which the centre temperature ($r = 0$, $z \approx H/2$) reaches $T_0 + 0.50\,(T_{furnace} - T_0)$ and $T_0 + 0.99\,(T_{furnace} - T_0)$, respectively.

## Requirements

- Python 3.8 or later (tested with Python X.Y)
- `numpy`, `scipy`, `pandas`, `matplotlib`
- `tkinter` (included with the standard Python installers on Windows and macOS; on Debian/Ubuntu install it with `sudo apt install python3-tk`)

```bash
pip install -r requirements.txt
```

## Usage

```bash
python pellet_fast_firing_gui.py
```

1. Fill in the input window and press **Run simulation**.
2. A window shows the summary of results.
3. You are asked whether to export a `.csv` file with the input conditions and the plotted data.
4. A figure with the temperature evolution is displayed.

Decimal commas are accepted in the input fields.

## Inputs

| Parameter | Symbol | Default | Unit |
|---|---|---|---|
| Thermal diffusivity | α | 0.31e-6 | m²/s |
| Thermal conductivity | κ | 2.0 | W/m/K |
| Initial temperature | T₀ | 25 | °C |
| Furnace temperature | T_furnace | 1200 | °C |
| Base (support) temperature | T_base | 1200 | °C |
| Pellet diameter | D | 10 | mm |
| Pellet thickness | H | 1.25 | mm |
| Convective coefficient | h | 25 | W/m²/K |
| Emissivity | ε | 0.80 | – |
| Simulation time | t_end | 60 | s |

The volumetric heat capacity is implicitly defined by $\rho c_p = \kappa / \alpha$.

## Outputs

- **Summary window:** inputs, final temperatures (centre, mean, base, top, side), final alumina conductivity and base conductance, $t_{50}$ and $t_{99}$.
- **Figure:** temperature versus time at the centre, mean, side, top and base of the pellet, with $t_{50}$ and $t_{99}$ marked.
- **CSV file (optional):** a header block with all input and internal parameters, followed by the columns `time_s`, `T_center_C`, `T_mean_C`, `T_side_C`, `T_top_C`, `T_bottom_C`, `Tfurnace_C`, `k_alumina_W_m_K`, `h_base_W_m2_K`.

## Fixed internal parameters

These values are set in `run_simulation()` and are not exposed in the GUI. Edit them in the code if needed.

| Parameter | Value |
|---|---|
| Gas time constant, τ_gas | 1.0 s |
| Radiation time constant, τ_rad | 0.2 s |
| Base contact resistance, R_c | 3.0 × 10⁻⁴ m²K/W |
| Alumina support thickness, δ | 1.0 mm |
| Grid, N_r × N_z | 15 × 15 |
| Stored time points | 300 |

## Limitations

- Pellet properties (α, κ) are set to be constant and independent of temperature. Experimentally, both can be measured at desired temperatures to improve accuracy.
- The grid is fixed and coarse; check mesh convergence if you change the geometry substantially.
- Radiation is exchanged with a single effective temperature (no view factors) and is linearised.
- The contact resistance and the support thickness are assumed values, not measured ones.
- Densification, shrinkage and reaction heats during sintering are not included.

## Citation

If you use this code in your work, please cite:

```
Ferrer-Nicomedes, S., Mormeneo-Segarra, A., Borrell, A., Vicente-Agut, A., Barba-Juan, A. Pellet heating model for fast firing (v1.0). Zenodo, 2026. https://doi.org/[DOI]
```

## License

Released under the [MIT License](LICENSE).
