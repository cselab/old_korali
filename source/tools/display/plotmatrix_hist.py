# *
# *  plotting.py
# *  PyPi4U
# *
# *  Authors:
# *     Philipp Mueller  - muellphi@ethz.ch
# *     Georgios Arampatzis - arampatzis@collegium.ethz.ch
# *
# *  Copyright 2018 ETH Zurich. All rights reserved.
# *


import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import scipy.interpolate
import argparse


def plot_histogram(ax, theta):
    """Plot histogram of theta to diagonal"""
    num_bins = 50
    for i in range(theta.shape[1]):
        hist, bins, _ = ax[i, i].hist(theta[:, i], num_bins, normed=1,
                                      color=(51/255, 1, 51/255), ec='black')
        if i == 0:

            # Rescale hist to scale of theta -> get correct axis titles
            hist = hist / np.max(hist) * (ax[i, i].get_xlim()[1] -
                                          ax[i, i].get_xlim()[0])
            bottom = ax[i, i].get_xlim()[0]

            widths = np.diff(bins)
            ax[i, i].cla()
            ax[i, i].bar(bins[:-1], hist, widths,
                         color=(51/255, 1, 51/255), ec='black', bottom=bottom)
            ax[i, i].set_ylim(ax[i, i].get_xlim())

            ax[i, i].set_xticklabels([])

        elif i == theta.shape[1] - 1:
            ax[i, i].set_yticklabels([])
        else:
            ax[i, i].set_xticklabels([])
            ax[i, i].set_yticklabels([])
        ax[i, i].tick_params(axis='both', which='both', length=0)


def plot_upper_triangle(ax, theta, lik=None):
    """Plot scatter plot to upper triangle of plots"""
    for i in range(theta.shape[1]):
        for j in range(i + 1, theta.shape[1]):
            if lik is None:
                ax[i, j].plot(theta[:, j], theta[:, i], '.', markersize=1)
            else:
                ax[i, j].scatter(theta[:, j], theta[:, i], marker='o', s=10,
                                 c=lik, facecolors='none', alpha=0.5)
            ax[i, j].set_xticklabels([])
            ax[i, j].set_yticklabels([])


def plot_lower_triangle(ax, theta):
    """Plot 2d histogram to lower triangle of plots"""
    for i in range(theta.shape[1]):
        for j in range(i):
            # returns bin values, bin edges and bin edges
            H, xe, ye = np.histogram2d(theta[:, j], theta[:, i], 8,
                                       normed=True)
            # plot and interpolate data
            ax[i, j].imshow(H.T, aspect="auto", interpolation='spline16',
                            origin='lower', extent=np.hstack((
                                                ax[j, j].get_xlim(),
                                                ax[i, i].get_xlim())),
                                                cmap=plt.get_cmap('jet'))
            if i < theta.shape[1]-1:
                ax[i, j].set_xticklabels([])
            if j > 0:
                ax[i, j].set_yticklabels([])


def plot_theta(file, likelihood=False):
    theta = np.loadtxt(file)
    fig, ax = plt.subplots(theta.shape[1]-2, theta.shape[1]-2)
    plot_histogram(ax, theta[:, :-2])
    if likelihood:
        plot_upper_triangle(ax, theta[:, :-2], theta[:, -2])
    else:
        plot_upper_triangle(ax, theta[:, :-2])
    plot_lower_triangle(ax, theta[:, :-2])
#    plt.jet()
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot generations.')
    parser.add_argument('filename', metavar='filename', help='Select file' +
                        ' for plotting.')
    parser.add_argument("-lik", "--likelihood", action="store_true",
                        help="Plot log-likelihood value")
    args = parser.parse_args()
    plot_theta(args.filename, args.likelihood)
