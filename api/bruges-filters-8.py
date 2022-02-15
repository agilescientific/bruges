import matplotlib.pyplot as plt
import bruges
w, t = bruges.filters.ricker(0.256, 0.002, 40)
plt.plot(t, w)