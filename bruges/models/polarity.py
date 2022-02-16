import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import io
import base64

from bruges.filters import ricker
from bruges.filters import rotate_phase


def _get_colour(cmap, frac):
    """
    Decide whether to make white or black labels.
    """
    cmap = get_cmap(cmap)
    return 'k' if (np.mean(cmap(frac)[:3]) > 0.5) else 'w'


def _make_synthetic(size=256, top=0.4, base=0.6, value=1, freq=25, phase=0):
    """Make a synthetic. Return the wavelet, the model, the RC, and the synthetic.
    """
    v = np.ones(size) - value

    v[int(top*size):int(base*size)] = value
    rc = np.diff(v)

    w = ricker(0.256, 0.001, freq)

    if phase != 0:
        w = rotate_phase(w, phase, degrees=True)

    syn = np.convolve(rc, w, mode='same')
    return w, v, rc, syn


def polarity_cartoon(layer='hard',
                     polarity='normal',
                     freq='med',
                     phase=0,
                     style='vd',
                     cmap=None,
                     fmt='png',
                     ):
    """
    Plot a polarity cartoon.
    """
    freqs = {'vhi': 60, 'hi': 30, 'med': 15, 'lo': 7.5,
             'vhigh': 60, 'high': 30, 'medium': 15, 'low': 7.5,
             'mod': 15, 'mid': 15}
    backgr = 'soft' if layer == 'hard' else 'hard'
    value = 1 if layer == 'hard' else 0
    size, top, base = 256, 0.4, 0.6

    _, v, _, syn = _make_synthetic(size, top, base, value, freq=freqs[freq], phase=phase)

    if polarity.lower() not in ['normal', 'seg', 'usa', 'us', 'canada']:
        syn *= -1

    if style == 'ramp':
        # cbar is a ramp.
        cbar = np.linspace(-1, 1, size).reshape(-1, 1)
    else:
        # cbar is the synthetic.
        cbar = syn.reshape(-1, 1)

    gs = {'width_ratios':[2,2,2,1]}
    fig, axs = plt.subplots(ncols=4,
                            figsize=(6, 4),
                            gridspec_kw=gs,
                            facecolor='w', sharey=True,
                            )

    # Plot the earth model.
    ax = axs[0]
    cmap_ = 'Greys'
    ax.imshow(v.reshape(-1, 1), aspect='auto', cmap=cmap_, vmin=-1.5, vmax=1.5)
    ax.axhline(top*size, c='w', lw=4)
    ax.axhline(base*size, c='w', lw=4)
    ax.axvline(0.55, c='w', lw=6)  # Total hack to adjust RHS.
    ax.text(0, size/4.75, backgr, ha='center', va='center', color=_get_colour(cmap_, (1-value)*256), size=25)
    ax.text(0, size/2+0.75, layer, ha='center', va='center', color=_get_colour(cmap_, (value)*256), size=25)

    # Plot the impedance diagram.
    ax = axs[1]
    cmap_ = 'Greys'
    ax.imshow(v.reshape(-1, 1), aspect='auto', cmap=cmap_, vmin=0, vmax=2)
    ax.axvline(-0.5, c=(0.58, 0.58, 0.58), lw=50)
    ax.text(0.45, 2*size/8, 'imp', ha='right', va='center', color='k', size=25)
    ax.annotate("", xy=(0.33, size/8), xytext=(0, size/8), arrowprops=dict(arrowstyle="->", connectionstyle="arc3"))

    # Plot the waveform.
    ax = axs[2]
    y = np.arange(syn.size)
    ax.plot(syn, y, 'k')
    ax.fill_betweenx(y, syn, 0, where=syn>0, color='k')
    ax.invert_yaxis()
    ax.text(0.65, size/8, '+', ha='center', va='center', size=30)
    ax.text(-0.65, size/8, '–', ha='center', va='center', size=40)

    # Plot the colourbar.
    ax = axs[3]
    cmap = cmap or 'gray'
    frac = 1/8
    top_col = _get_colour(cmap, frac)
    bot_col = _get_colour(cmap, 7*frac)
    ax.imshow(cbar, cmap=cmap, aspect='auto')
    if style == 'ramp':
        ax.text(0, frac*size, '+', ha='center', va='center', color=top_col, size=30)
        ax.text(0, 7*frac*size, '–', ha='center', va='center', color=bot_col, size=40)

    # Make final adjustments to subplots and figure.
    for ax in axs:
        ax.set_axis_off()

    # Hack.
    plt.subplots_adjust(left=0.1)

    return fig, axs
