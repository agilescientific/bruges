"""
===========
Klauder Waveform
===========

Create a simple plot of the Klauder signal.
"""
import matplotlib.pyplot as plt
import bruges as br

# Data for plotting
wavelet, t = bruges.filters.klauder(0.5, 0.002, [10, 80], return_t=True)

# Note that using plt.subplots below is equivalent to using
# fig = plt.figure and then ax = fig.add_subplot(111)
fig, ax = plt.subplots(t, wavelet)
ax.plot()
ax.set(xlabel='time (s)', ylabel='amplitude',
       title='Klauder wavelet')
ax.grid()
plt.show()

