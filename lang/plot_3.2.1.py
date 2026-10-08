import matplotlib.pyplot as plt
import numpy as np

alphas = ["0.01", "0.10", "1.00"]

print(f"{'Alpha':<8} | {'10-Block sigma_E':<18} | {'100-Block sigma_E':<18}")
print("-" * 50)

# Dictionaries to store energy time series for plotting
data_10 = {}
data_100 = {}

for alpha in alphas:
  f_10 = f"lang/efile/alpha{alpha}_10bl"
  f_100 = f"lang/efile/alpha{alpha}_100bl"

  d_10 = np.loadtxt(f_10)
  d_100 = np.loadtxt(f_100)

  e_10 = d_10[:, -1] if d_10.ndim > 1 else d_10
  e_100 = d_100[:, -1] if d_100.ndim > 1 else d_100

  data_10[alpha] = e_10
  data_100[alpha] = e_100

  s_10 = np.std(e_10)
  s_100 = np.std(e_100)

  print(f"{alpha:<8} | {s_10:<18.5f} | {s_100:<18.5f}")

# Plotting comparison
fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharey=True)

# 10-Block Plot
for alpha in alphas:
  axes[0].plot(data_10[alpha], label=f"alpha = {alpha}", linewidth=0.8)
axes[0].set_title("10-Block Runs")
axes[0].set_ylabel("Total Energy per Particle (E)")
axes[0].legend()
axes[0].grid(True)

# 100-Block Plot
for alpha in alphas:
  axes[1].plot(data_100[alpha], label=f"alpha = {alpha}", linewidth=0.8)
axes[1].set_title("100-Block Runs")
axes[1].set_xlabel("Sample Index")
axes[1].set_ylabel("Total Energy per Particle (E)")
axes[1].legend()
axes[1].grid(True)

plt.tight_layout()
plt.savefig("lang/alpha_comparison_blocks.png")
plt.show()