import matplotlib.pyplot as plt
import bruges
w, t = bruges.filters.ormsby(0.256, 0.002, [5, 10, 40, 80])
plt.plot(t, w)