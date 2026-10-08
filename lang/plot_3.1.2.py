import matplotlib.pyplot as plt
# the files in efile/ contain time and total energy per particle
filename = "lang/efile/0064_r0.500_T1.000_alpha0.00_dt016"

steps = []
energies = []

with open(filename, "r") as f:
  for i, line in enumerate(f):
    parts = line.split()
    if parts:
      # Time (step)
      steps.append(i)
      # Total energy
      energies.append(float(parts[0]))
      

# Plot E vs t
plt.figure(figsize=(8, 5))
plt.plot(steps, energies, color="r", linewidth=1)
plt.xlim(0, len(steps))
plt.ylim(-0.4, 1.0)
plt.xlabel("Sample / Time Step")
plt.ylabel("Total Energy per Particle (E)")
plt.title("Numerical Instability Test ($\Delta t = 0.016$)")
plt.grid(True)
plt.savefig("instability_plot.png")
plt.show()