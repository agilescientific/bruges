import matplotlib.pyplot as plt
import bruges
w, t = bruges.filters.generalized(0.256, 0.002, 40, u=1.0)
plt.plot(t, w)