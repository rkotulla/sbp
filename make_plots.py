#!/usr/bin/env python3


import os
import sys
import pandas
import numpy
import matplotlib
import matplotlib.pyplot as plt

figsize = (5,3)


def make_plot_ellipticity(data, formats=None, pixelscale=1.0, basename='sbp'):

    fig = plt.figure(figsize=figsize, dpi=150)
    ax = fig.add_subplot(111)
    # ax.tight_layout()

    if (pixelscale is None):
        ps = 1.0
        ax.set_xlabel("effective radius [pixels]")
    else:
        ps = pixelscale
        ax.set_xlabel("effective radius [arcsec]")


    # ax.scatter(data['sma']*ps, data['ellipticity'])
    ax.scatter(data['eff_radius'], data['ellipticity'])
    ax.set_xscale('log')

    ax.set_ylabel("ellipticity")

    for f in formats:
        fig.savefig("%s__ellipticity.%s" % (basename, f), bbox_inches='tight')

    pass





def make_plot_positionangle(data, formats=None, pixelscale=1.0, basename='sbp'):

    fig = plt.figure(figsize=figsize, dpi=150)
    ax = fig.add_subplot(111)
    # ax.tight_layout()

    if (pixelscale is None):
        ps = 1.0
        ax.set_xlabel("effective radius [pixels]")
    else:
        ps = pixelscale
        ax.set_xlabel("effective radius [arcsec]")


    ax.errorbar(data['eff_radius'], data['pa'], yerr=data['pa_err'], fmt="o", marker="o")
    ax.scatter(data['eff_radius'], data['pa'], marker="o")

    ax.set_xscale('log')

    ax.set_ylabel("position angle [degrees]")

    for f in formats:
        fig.savefig("%s__position_angle.%s" % (basename, f), bbox_inches='tight')

    pass

def make_plot_intensity(data, formats=None, pixelscale=1.0, basename='sbp'):

    fig = plt.figure(figsize=figsize, dpi=150)
    ax = fig.add_subplot(111)
    # ax.tight_layout()

    if (pixelscale is None):
        ps = 1.0
        ax.set_xlabel("effective radius [pixels]")
    else:
        ps = pixelscale
        ax.set_xlabel("effective radius [arcsec]")

    eff_radius = data['sma'] * numpy.sqrt(1.-data['ellipticity']) * ps
    ax.errorbar(eff_radius, data['intens_bgsub'], yerr=data['intens_err'], fmt="o", marker="o")
    ax.scatter(eff_radius, data['intens_bgsub'], marker="o")

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_ylabel("position angle [degrees]")

    for f in formats:
        fig.savefig("%s__intensity.%s" % (basename, f), bbox_inches='tight')

    pass



def make_plot_surfacebrightness(data, formats=None, pixelscale=1.0, basename='sbp'):

    fig = plt.figure(figsize=figsize, dpi=150)
    ax = fig.add_subplot(111)
    # ax.tight_layout()

    if (pixelscale is None):
        ps = 1.0
        ax.set_xlabel("effective radius [pixels]")
    else:
        ps = pixelscale
        ax.set_xlabel("effective radius [arcsec]")

    ax.errorbar(data['eff_radius'], data['mag_arcsec'], yerr=data['mag_arcsec_err'], fmt="o", marker="o")
    ax.scatter(data['eff_radius'], data['mag_arcsec'], marker="o")

    ax.set_xscale('log')
    ax.set_ylim(ax.get_ylim()[::-1])
    # ax.set_yscale('log')

    ax.set_ylabel(r"surface brightness [mag arcsec$^{-2}$]")

    for f in formats:
        fig.savefig("%s__surface_brightness.%s" % (basename, f), bbox_inches='tight')

    pass



def make_plot(data, formats=None, pixelscale=1.0, zeropoint=0, basename="sbp"):

    if (formats is None):
        formats = ['png', 'pdf']

    background_level = 0.97 * numpy.min(data['intens'])
    background_error = 0.1*background_level

    data['intens_bgsub'] = data['intens'] - background_level
    data['intens_error_full'] = numpy.hypot(data['intens_err'], background_error)

    mag_arcsec = -2.5*numpy.log10(data['intens_bgsub'] / pixelscale**2) + zeropoint
    data['mag_arcsec'] = mag_arcsec
    data['mag_arcsec_err'] = data['intens_error_full'] / data['intens'] * 5
    print(numpy.array(data['mag_arcsec_err']))

    data['eff_radius'] = data['sma'] * numpy.sqrt(1. - data['ellipticity']) * pixelscale

    make_plot_ellipticity(data, formats, pixelscale, basename)
    make_plot_positionangle(data, formats, pixelscale, basename)
    make_plot_intensity(data, formats, pixelscale, basename)
    make_plot_surfacebrightness(data, formats, pixelscale, basename)


if __name__ == "__main__":

    fn = sys.argv[1]

    bn,ext = os.path.splitext(fn)
    df = pandas.read_csv(fn)
    make_plot(df, pixelscale=0.168, zeropoint=26., basename=bn)

