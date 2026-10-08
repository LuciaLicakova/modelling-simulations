import matplotlib.pyplot as plt

dt = []
e_mean = []
sigma = []

with open("lang/dt_study.dat", "r") as f:
  for line in f:
    parts = line.split()
    if parts:
      dt.append(float(parts[0]))
      e_mean.append(float(parts[1]))
      sigma.append(float(parts[2]))

# Plot Average Energy vs dt
plt.figure()
plt.plot(dt, e_mean, marker='o', linestyle="-", color="b")
plt.xlabel("Delta t")
plt.ylabel("<E>")
plt.title("Average Energy vs Time Step")
plt.grid(True)
plt.savefig("mean_energy.png")

# Plot Sigma E vs dt
plt.figure()
plt.plot(dt, sigma, marker='o', linestyle="-", color="r")
plt.xlabel("Delta t")
plt.ylabel("Sigma E")
plt.title("Energy Fluctuation vs Time Step")
plt.grid(True)
plt.savefig("sigma_e.png")

plt.show()