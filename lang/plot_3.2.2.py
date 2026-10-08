import matplotlib.pyplot as plt
import numpy as np

alphas = ["0.01", "0.10", "1.00"]

plt.figure(figsize=(10, 6))

for alpha in alphas:
  filename = f"lang/efile/alpha{alpha}_highres"
  data = np.loadtxt(filename)
  energies = data[:, -1] if data.ndim > 1 else data
  time_steps = np.arange(len(energies))

  plt.plot(time_steps, energies, label=f"alpha = {alpha}", alpha=0.7, linewidth=0.5)

plt.xlabel("Measurement Step")
plt.ylabel("Total Energy per Particle (E)")
plt.title("High-Resolution Energy Fluctuations vs Time")
plt.legend()
plt.grid(True)
plt.savefig("lang/energy_vs_time_highres.png")
#plt.show()

colors = ["tab:blue", "tab:orange", "tab:green"]

fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)

for i, alpha in enumerate(alphas):
  filename = f"lang/efile/alpha{alpha}_highres"
  data = np.loadtxt(filename)
  energies = data[:, -1] if data.ndim > 1 else data
  time_steps = np.arange(len(energies))

  axes[i].plot(
      time_steps, energies, color=colors[i], linewidth=0.5, label=f"alpha = {alpha}"
  )
  axes[i].set_ylabel("Total Energy (E)")
  axes[i].legend(loc="upper right")
  axes[i].grid(True)

axes[-1].set_xlabel("Measurement Step")
fig.suptitle(
    "High-Resolution Energy Fluctuations by Friction Coefficient",
    fontsize=12,
    y=0.95,
)

plt.tight_layout()
plt.savefig("lang/energy_vs_time_subplots.png")
plt.show()