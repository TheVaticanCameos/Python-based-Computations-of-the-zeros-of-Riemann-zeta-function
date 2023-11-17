import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from check_RH import compute_zeta_AS, compute_Zeta_AS


def draw_zeta_critical_line(begin, end, step=0.1):
    """
    This module draws the dynamical process of zeta(1/2+it) on the critical line.
    :param begin: min of t
    :param end: max of t
    :param step: step of sampling; default: 0.1
    :return:none
    """
    fig = plt.figure(figsize=(10, 10))
    plt.xlabel(r'$\mathfrak{R}\zeta(\frac{1}{2}+it)$')
    plt.ylabel(r'$\mathfrak{I}\zeta(\frac{1}{2}+it)$')
    plt.xlim(-2, 4)
    plt.ylim(-3, 3)
    plt.grid(linestyle=':')
    plt.title('Dynamical '+r'$\zeta(s)$'+' on the critical line')
    x, y = [], []

    def update(n):
        x.append(compute_zeta_AS(n).real)
        y.append(compute_zeta_AS(n).imag)
        plt.plot(x, y, "r-", linewidth=0.5, label=r'$\zeta(\frac{1}{2}+it)$')

    ani = FuncAnimation(fig, update, frames=np.arange(begin, end, step), interval=50, blit=False, repeat=False)
    # ani.save("dynamical_zeta.gif", writer='pillow')
    plt.show()

    return


def compare(begin, end, step=0.1):
    """
    This module compares the absolute value of zeta function and Riemann-Siegel Z function, and visualizes them.
    :param begin: min of t(should be no less than 0.2)
    :param end: max of t
    :param step: step of sampling; default: 0.1
    :return: none
    """
    plt.figure(figsize=(8, 4), dpi=200)
    plt.style.use('seaborn-paper')
    plt.grid(linestyle=':')
    plt.title("Comparison between "+r'$|\zeta(\frac{1}{2}+it)\vert$'+' and '+r'$Z(\frac{1}{2}+it)$')
    x = np.arange(begin, end, step)
    y1 = np.array([compute_Zeta_AS(t) for t in x])
    y2 = np.array([np.abs(compute_zeta_AS(t)) for t in x])
    plt.plot(x, y1, linestyle='-', linewidth=0.5, color='r', alpha=0.7, label=r'$Z(\frac{1}{2}+it)$')
    plt.plot(x, y2, linestyle='--', linewidth=0.5, color='k', alpha=1, label=r'$|\zeta(\frac{1}{2}+it)\vert$')
    plt.xlabel(r'$t=\mathfrak{R}(s)$')
    plt.ylabel('value')
    plt.legend()
    plt.show()
    return


def zeta_critical_line_real_imag(begin, end, step=0.1):
    """
    This module draws the figure of the real and imaginary parts of zeta(1/2+it).
    :param begin: min of t(should be no less than 0.4)
    :param end: max of t
    :param step: step of sampling; default: 0.1
    :return: none
    """
    x = np.arange(begin, end, step)
    real = np.array([compute_zeta_AS(t).real for t in x])
    imag = np.array([compute_zeta_AS(t).imag for t in x])
    color_real = 'b'
    color_imag = 'r'
    title = 'Comparison between ' + r'$\mathfrak{R}\zeta(\frac{1}{2}+it)$' + \
            ' and ' + r'$\mathfrak{I}\zeta(\frac{1}{2}+it)$'
    fig = plt.figure(figsize=(14, 4), dpi=250)
    plt.style.use('seaborn-paper')

    subfigs = fig.subfigures(1, 2, width_ratios=[0.62, 0.37])
    subfigs[0].suptitle(title)
    ax1 = subfigs[0].subplots()
    ax1.set_xlabel('t')
    ax1.set_ylabel(r'$\mathfrak{R}\zeta(\frac{1}{2}+it)$'+' and '+r'$\mathfrak{I}\zeta(\frac{1}{2}+it)$')
    ax1.plot(x, real, linestyle='--', linewidth=0.5, color=color_real, label=r'$\mathfrak{R}\zeta(\frac{1}{2}+it)$')
    ax1.tick_params(axis='x', labelcolor='black')
    ax1.tick_params(axis='y', labelcolor='black')
    ax1.plot(x, imag, linestyle='--', linewidth=0.5, color=color_imag, label=r'$\mathfrak{I}\zeta(\frac{1}{2}+it)$')
    ax1.grid(linestyle=':')
    ax1.legend()

    ax2 = subfigs[1].subplots(2, 1, sharex=True)
    ax2[0].set_title(r'$\mathfrak{R}\zeta(\frac{1}{2}+it)$')
    ax2[0].grid(linestyle=':')
    ax2[0].plot(x, real, linestyle='--', linewidth=0.5, color=color_real, label=r'$\mathfrak{R}\zeta(\frac{1}{2}+it)$')
    ax2[0].yaxis.set_label_position("right")
    ax2[0].set_ylabel('value')
    ax2[0].legend()
    ax2[1].set_title(r'$\mathfrak{I}\zeta(\frac{1}{2}+it)$')
    ax2[1].grid(linestyle=':')
    ax2[1].plot(x, imag, linestyle='--', linewidth=0.5, color=color_imag, label=r'$\mathfrak{I}\zeta(\frac{1}{2}+it)$')
    ax2[1].yaxis.set_label_position("right")
    ax2[1].set_ylabel('value')
    ax2[1].legend()

    plt.show()

    return


# compare(0.2, 50)
# zeta_critical_line_real_imag(0.4, 50)
draw_zeta_critical_line(1, 50)
