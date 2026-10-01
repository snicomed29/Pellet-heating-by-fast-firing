# pellet_heating_gui.py

import numpy as np
import pandas as pd
import tkinter as tk
from tkinter import messagebox, filedialog
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


def k_alumina(T_C):
    """
    Alumina thermal conductivity - Touloukian correlation for dense polycrystalline alumina.
    Monotonically decreasing with temperature, physically correct.
    T_C in °C, k in W/m/K.
    Valid range: 25 - 1200 °C (approx. 78 W/m/K at 25°C, ~16 W/m/K at 1200°C).
    """
    T_K = np.clip(T_C, 0, 2000) + 273.15
    k = 85.686 * (T_K / 273.15) ** (-1.01)
    return max(float(k), 2.0)   # safety floor: 2 W/m/K


class InputDialog:
    def __init__(self, root):
        self.root = root
        self.root.title("Thermal model input data")
        self.values = None

        self.fields = [
            ("Thermal diffusivity α [m²/s]", "0.31e-6"),
            ("Thermal conductivity κ [W/m/K]", "2.0"),
            ("Initial temperature [°C]", "25"),
            ("Furnace temperature [°C]", "1200"),
            ("Base temperature  [°C]", "1200"),
            ("Diameter D [mm]", "10"),
            ("Thickness H [mm]", "1.5"),
            ("Convective coefficient h [W/m²/K]", "25"),
            ("Emissivity ε [-]", "0.80"),
            ("Simulation time [s]", "60"),
        ]

        self.entries = []

        for i, (label, default) in enumerate(self.fields):
            tk.Label(root, text=label, anchor="w").grid(row=i, column=0, sticky="w", padx=8, pady=4)
            entry = tk.Entry(root, width=22)
            entry.insert(0, default)
            entry.grid(row=i, column=1, padx=8, pady=4)
            self.entries.append(entry)

        button_frame = tk.Frame(root)
        button_frame.grid(row=len(self.fields), column=0, columnspan=2, pady=12)

        tk.Button(button_frame, text="Run simulation", command=self.submit).pack(side="left", padx=8)
        tk.Button(button_frame, text="Cancel", command=self.cancel).pack(side="left", padx=8)

    def submit(self):
        try:
            raw = [entry.get().replace(",", ".") for entry in self.entries]

            self.values = {
                "alpha":      float(raw[0]),
                "k":          float(raw[1]),
                "T0_C":       float(raw[2]),
                "Tfurnace_C": float(raw[3]),
                "Tbase_C":    float(raw[4]),
                "D_mm":       float(raw[5]),
                "H_mm":       float(raw[6]),
                "h_conv":     float(raw[7]),
                "epsilon":    float(raw[8]),
                "t_end":      float(raw[9]),

                # Fixed mesh (internal)
                "Nr": 15,
                "Nz": 15,
            }

            self.validate()
            self.root.quit()
            self.root.destroy()

        except Exception as e:
            messagebox.showerror("Input error", f"Please check the input values:\n\n{e}")

    def cancel(self):
        self.values = None
        self.root.quit()
        self.root.destroy()

    def validate(self):
        v = self.values

        if v["alpha"] <= 0:
            raise ValueError("alpha must be > 0")
        if v["k"] <= 0:
            raise ValueError("k must be > 0")
        if v["D_mm"] <= 0:
            raise ValueError("pellet diameter must be > 0")
        if v["H_mm"] <= 0:
            raise ValueError("pellet thickness must be > 0")
        if v["h_conv"] < 0:
            raise ValueError("h_conv must be >= 0")
        if not (0 <= v["epsilon"] <= 1):
            raise ValueError("epsilon must be between 0 and 1")
        if v["t_end"] <= 0:
            raise ValueError("total simulation time must be > 0")


def rhs_pellet_axisym(
    t, y, Nr, Nz, r, dr, dz, alpha, k,
    h_conv, epsilon, sigma,
    Tfurnace, Tbase, T0, tau_gas, tau_rad,
    R_contact, delta_alumina
):
    T = y.reshape((Nr, Nz))
    T = np.maximum(T, 1.0)

    dTdt = np.zeros((Nr, Nz))

    # Short local transient, assuming the furnace is already stabilized
    Tgas = T0 + (Tfurnace - T0) * (1 - np.exp(-t / tau_gas))
    Trad = T0 + (Tfurnace - T0) * (1 - np.exp(-t / tau_rad))

    Tgas = max(Tgas, 1.0)
    Trad = max(Trad, 1.0)
    Tbase = max(Tbase, 1.0)

    for i in range(Nr):
        for j in range(Nz):

            Tij = T[i, j]

            h_rad = epsilon * sigma * (Tij + Trad) * (Tij**2 + Trad**2)

            # =====================================================
            # Radial term
            # =====================================================
            if i == 0:
                radial_term = 4 * (T[1, j] - T[0, j]) / dr**2

            elif i == Nr - 1:
                q_side = h_conv * (Tij - Tgas) + h_rad * (Tij - Trad)
                Tghost_r = T[i - 1, j] - 2 * dr * (q_side / k)
                d2Tdr2 = (Tghost_r - 2 * T[i, j] + T[i - 1, j]) / dr**2
                dTdr_over_r = (1 / r[i]) * (Tghost_r - T[i - 1, j]) / (2 * dr)
                radial_term = d2Tdr2 + dTdr_over_r

            else:
                d2Tdr2 = (T[i + 1, j] - 2 * T[i, j] + T[i - 1, j]) / dr**2
                dTdr_over_r = (1 / r[i]) * (T[i + 1, j] - T[i - 1, j]) / (2 * dr)
                radial_term = d2Tdr2 + dTdr_over_r

            # =====================================================
            # Axial term
            # =====================================================
            if j == 0:
                # Base: Touloukian alumina conductivity + medium contact resistance
                Tref_C = float(Tij) - 273.15
                k_al = k_alumina(Tref_C)
                h_base = 1.0 / (R_contact + delta_alumina / k_al)

                q_bottom = h_base * (Tij - Tbase)
                Tghost_z = T[i, j + 1] - 2 * dz * (q_bottom / k)
                axial_term = (T[i, j + 1] - 2 * T[i, j] + Tghost_z) / dz**2

            elif j == Nz - 1:
                q_top = h_conv * (Tij - Tgas) + h_rad * (Tij - Trad)
                Tghost_z = T[i, j - 1] - 2 * dz * (q_top / k)
                axial_term = (Tghost_z - 2 * T[i, j] + T[i, j - 1]) / dz**2

            else:
                axial_term = (T[i, j + 1] - 2 * T[i, j] + T[i, j - 1]) / dz**2

            dTdt[i, j] = alpha * (radial_term + axial_term)

    return dTdt.ravel()


def run_simulation(params):
    alpha      = params["alpha"]
    k          = params["k"]
    T0_C       = params["T0_C"]
    Tfurnace_C = params["Tfurnace_C"]
    Tbase_C    = params["Tbase_C"]
    D_mm       = params["D_mm"]
    H_mm       = params["H_mm"]
    h_conv     = params["h_conv"]
    epsilon    = params["epsilon"]
    t_end      = params["t_end"]
    Nr         = params["Nr"]
    Nz         = params["Nz"]

    R = (D_mm * 1e-3) / 2
    H = H_mm * 1e-3

    sigma          = 5.670374419e-8
    nTimesToStore  = 300
    tau_gas        = 1.0
    tau_rad        = 0.2
    R_contact      = 3.0e-4    # medium contact resistance [m²K/W]
    delta_alumina  = 1.0e-3    # alumina support thickness [m]

    T0       = T0_C       + 273.15
    Tfurnace = Tfurnace_C + 273.15
    Tbase    = Tbase_C    + 273.15

    r = np.linspace(0, R, Nr)
    z = np.linspace(0, H, Nz)

    dr = r[1] - r[0]
    dz = z[1] - z[0]

    y0 = (T0 * np.ones((Nr, Nz))).ravel()

    t_eval   = np.linspace(0, t_end, nTimesToStore)
    max_step = min(0.10, t_end / 400)

    sol = solve_ivp(
        fun=lambda t, y: rhs_pellet_axisym(
            t, y, Nr, Nz, r, dr, dz, alpha, k,
            h_conv, epsilon, sigma,
            Tfurnace, Tbase, T0, tau_gas, tau_rad,
            R_contact, delta_alumina
        ),
        t_span=(0, t_end),
        y0=y0,
        method="BDF",
        t_eval=t_eval,
        rtol=1e-6,
        atol=1e-7,
        max_step=max_step
    )

    if not sol.success:
        raise RuntimeError(sol.message)

    tSol = sol.t
    ySol = sol.y.T
    nT   = len(tSol)

    Tsol = np.zeros((Nr, Nz, nT))
    for it in range(nT):
        Tsol[:, :, it] = ySol[it, :].reshape((Nr, Nz))

    iz_center = np.argmin(np.abs(z - H / 2))
    T_center  = Tsol[0, iz_center, :] - 273.15

    T_mean   = np.zeros(nT)
    T_side   = np.zeros(nT)
    T_top    = np.zeros(nT)
    T_bottom = np.zeros(nT)

    for it in range(nT):
        T_mean[it]   = np.mean(Tsol[:, :, it]) - 273.15
        T_side[it]   = np.mean(Tsol[-1, :, it]) - 273.15
        T_top[it]    = np.mean(Tsol[:, -1, it]) - 273.15
        T_bottom[it] = np.mean(Tsol[:, 0, it])  - 273.15

    Tgas_hist_C = (T0 + (Tfurnace - T0) * (1 - np.exp(-tSol / tau_gas))) - 273.15
    Trad_hist_C = (T0 + (Tfurnace - T0) * (1 - np.exp(-tSol / tau_rad))) - 273.15

    # Post-process: alumina conductivity and effective base conductance
    k_alumina_hist = np.zeros(nT)
    h_base_hist    = np.zeros(nT)

    for it in range(nT):
        Tref_C = float(np.clip(T_bottom[it], 0, 2000))
        k_al   = k_alumina(Tref_C)
        k_alumina_hist[it] = k_al
        h_base_hist[it]    = 1.0 / (R_contact + delta_alumina / k_al)

    T_target_50 = T0_C + 0.50 * (Tfurnace_C - T0_C)
    T_target_99 = T0_C + 0.99 * (Tfurnace_C - T0_C)

    idx50 = np.where(T_center >= T_target_50)[0]
    idx99 = np.where(T_center >= T_target_99)[0]

    t50 = tSol[idx50[0]] if len(idx50) > 0 else np.nan
    t99 = tSol[idx99[0]] if len(idx99) > 0 else np.nan

    return {
        "params":          params,
        "R":               R,
        "H":               H,
        "r":               r,
        "z":               z,
        "tSol":            tSol,
        "Tsol":            Tsol,
        "T_center":        T_center,
        "T_mean":          T_mean,
        "T_side":          T_side,
        "T_top":           T_top,
        "T_bottom":        T_bottom,
        "Tgas_hist_C":     Tgas_hist_C,
        "Trad_hist_C":     Trad_hist_C,
        "k_alumina_hist":  k_alumina_hist,
        "h_base_hist":     h_base_hist,
        "t50":             t50,
        "t99":             t99,
        "tau_gas":         tau_gas,
        "tau_rad":         tau_rad,
        "R_contact":       R_contact,
        "delta_alumina":   delta_alumina,
        "Tfurnace_C":      Tfurnace_C,
    }


def show_results(results):
    p = results["params"]

    def fmt_time(x):
        return "Not reached" if np.isnan(x) else f"{x:.2f} s"

    text = f"""\
================ RESULTS ================

Pellet diameter          = {p['D_mm']:.4f} mm
Pellet thickness         = {p['H_mm']:.4f} mm
Thermal diffusivity α    = {p['alpha']:.4e} m²/s
Thermal conductivity κ   = {p['k']:.4f} W/m/K
Initial temperature T₀   = {p['T0_C']:.2f} °C
Furnace temperature      = {p['Tfurnace_C']:.2f} °C
Base temperature         = {p['Tbase_C']:.2f} °C
Convective coefficient h = {p['h_conv']:.2f} W/m²/K
Emissivity ε             = {p['epsilon']:.2f}

τ_gas                    = {results['tau_gas']:.2f} s
τ_rad                    = {results['tau_rad']:.2f} s

Base contact             = MEDIUM
Contact resistance R     = {results['R_contact']:.3e} m²K/W
δ_Alumina                = {results['delta_alumina']:.3e} m
κ_Alumina (final)        = {results['k_alumina_hist'][-1]:.3f} W/m/K  [Touloukian]
h_base (final)           = {results['h_base_hist'][-1]:.2f} W/m²/K

Final T_center           = {results['T_center'][-1]:.2f} °C
Final T_mean             = {results['T_mean'][-1]:.2f} °C
Final T_base             = {results['T_bottom'][-1]:.2f} °C
Final T_top              = {results['T_top'][-1]:.2f} °C
Final T_side             = {results['T_side'][-1]:.2f} °C

t50_center               = {fmt_time(results['t50'])}
t99_center               = {fmt_time(results['t99'])}

========================================="""

    messagebox.showinfo("Results", text)


def export_csv(results):
    answer = messagebox.askyesno(
        "Export data",
        "Do you want to export a .csv file with the plotted data and input conditions?"
    )

    if not answer:
        return

    filename = filedialog.asksaveasfilename(
        defaultextension=".csv",
        filetypes=[("CSV files", "*.csv")],
        initialfile="pellet_thermal_evolution_data.csv",
        title="Save thermal evolution data"
    )

    if not filename:
        return

    p = results["params"]

    rows_conditions = [
        ["INPUT CONDITIONS", ""],
        ["alpha_m2_s",               p["alpha"]],
        ["k_pellet_W_m_K",           p["k"]],
        ["T0_C",                     p["T0_C"]],
        ["Tfurnace_C",               p["Tfurnace_C"]],
        ["Tbase_C",                  p["Tbase_C"]],
        ["D_mm",                     p["D_mm"]],
        ["H_mm",                     p["H_mm"]],
        ["h_conv_W_m2_K",            p["h_conv"]],
        ["epsilon",                  p["epsilon"]],
        ["t_end_s",                  p["t_end"]],
        ["Nr",                       p["Nr"]],
        ["Nz",                       p["Nz"]],
        ["tau_gas_s",                results["tau_gas"]],
        ["tau_rad_s",                results["tau_rad"]],
        ["R_contact_m2K_W",          results["R_contact"]],
        ["delta_alumina_m",          results["delta_alumina"]],
        ["k_alumina_correlation",    "Touloukian"],
        ["t50_center_s",             results["t50"]],
        ["t99_center_s",             results["t99"]],
        ["", ""],
    ]

    df_data = pd.DataFrame({
        "time_s":           results["tSol"],
        "T_center_C":       results["T_center"],
        "T_mean_C":         results["T_mean"],
        "T_side_C":         results["T_side"],
        "T_top_C":          results["T_top"],
        "T_bottom_C":       results["T_bottom"],
        "Tfurnace_C":       np.full_like(results["tSol"], results["Tfurnace_C"]),
        "k_alumina_W_m_K":  results["k_alumina_hist"],
        "h_base_W_m2_K":    results["h_base_hist"],
    })

    # newline="" prevents pandas double line endings on Windows
    with open(filename, "w", encoding="utf-8", newline="") as f:
        for row in rows_conditions:
            f.write(f"{row[0]},{row[1]}\n")

        f.write("PLOTTED_DATA\n")
        df_data.to_csv(f, index=False)

    messagebox.showinfo("Export completed", f"CSV file saved:\n\n{filename}")


def plot_figures(results):
    t = results["tSol"]

    plt.figure(figsize=(9, 6))

    plt.plot(t, results["T_center"], linewidth=2,   label="Center")
    plt.plot(t, results["T_mean"],   "--", linewidth=1, label="Mean")
    plt.plot(t, results["T_side"],   "--", linewidth=1, label="Side")
    plt.plot(t, results["T_top"],    "--", linewidth=1, label="Top")
    plt.plot(t, results["T_bottom"], "--", linewidth=1, label="Base")

    if not np.isnan(results["t50"]):
        plt.axvline(results["t50"], color="salmon", linewidth=1.2, linestyle="-", label="t50")

    if not np.isnan(results["t99"]):
        plt.axvline(results["t99"], color="maroon", linewidth=1.2, linestyle="-", label="t99")

    plt.xlabel("Time [s]")
    plt.ylabel("Temperature [°C]")
    plt.title("Pellet thermal evolution")
    plt.grid(True)
    plt.legend()
    plt.minorticks_on()
    plt.tick_params(axis="both", which="major", length=6, width=1)
    plt.tick_params(axis="both", which="minor", length=3, width=0.8)
    plt.tight_layout()
    plt.show()


def main():
    root = tk.Tk()
    dialog = InputDialog(root)
    root.mainloop()

    params = dialog.values

    if params is None:
        return

    try:
        results = run_simulation(params)
    except Exception as e:
        messagebox.showerror("Simulation error", str(e))
        return

    show_results(results)
    export_csv(results)
    plot_figures(results)


if __name__ == "__main__":
    main()
