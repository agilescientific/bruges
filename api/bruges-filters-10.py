import matplotlib.pyplot as plt
import bruges
w, t = bruges.filters.klauder(0.256, 0.002, [12, 48])
plt.plot(t, w)