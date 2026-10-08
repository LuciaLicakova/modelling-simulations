import subprocess
import re
import matplotlib.pyplot as plt

def run_sim(cmd):
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    output = result.stdout + result.stderr
    
    # Parse Potential E: value +/- error
    match_pe = re.search(r"Potential E:\s*([-\d.]+)\s*\+/-\s*([-\d.]+)", output)
    pe = float(match_pe.group(1)) if match_pe else None
    err = float(match_pe.group(2)) if match_pe else None
    return pe, err

# Run Monte Carlo (delta_t / alpha = 0)
print("Running Monte Carlo simulation...")
pe_mc, err_mc = run_sim("cd mc && ./sim N=64 rho=0.3 T=1.0 run")

# Run Brownian dynamics at various dt values with alpha = 10.0
dts = [0.001, 0.002, 0.003, 0.004]
alpha = 10.0
x_vals = [0.0] + [dt / alpha for dt in dts]
y_vals = [pe_mc]
y_errs = [err_mc]

for dt in dts:
    print(f"Running Brownian dynamics with dt={dt}...")
    pe, err = run_sim(f"cd brown && ./sim N=64 rho=0.3 T=1.0 alpha={alpha} deltat={dt} run")
    y_vals.append(pe)
    y_errs.append(err)

# Plotting
plt.figure(figsize=(8, 5))
plt.errorbar(x_vals, y_vals, yerr=y_errs, fmt='o-', capsize=5, color='b', label='Simulations')
plt.axvline(0, color='gray', linestyle='--', alpha=0.5)

plt.xlabel(r'$\Delta t / \alpha$', fontsize=12)
plt.ylabel('Potential Energy', fontsize=12)
plt.title(r'Exercise 5.1: Potential Energy vs $\Delta t / \alpha$', fontsize=14)
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()

plt.savefig('exercise_5_1_plot.png', dpi=300)
print("Simulation complete. Plot saved as exercise_5_1_plot.png.")
plt.show()