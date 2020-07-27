#!/usr/bin/env python3


import os
import sys
import astropy.io.fits as pyfits
import numpy
import pandas
import photutils
import photutils.isophote
from astropy.io import votable
import astropy.wcs
import argparse

if __name__ == "__main__":

    cmdline = argparse.ArgumentParser()
    #
    # required parameters
    #
    cmdline.add_argument("input_fn", help="input image")
    cmdline.add_argument("catalog_fn", help='catalog filename [must on VOTable format]')
    cmdline.add_argument("source_id", type=int, help='source id from sextractor')
    cmdline.add_argument("segmentation_fn", help="filename of segmentation frame")
    #
    # optional parameters
    #
    cmdline.add_argument("-o", "--output", dest="output_basename", type=str,
                         default="sbptest",
                         help="base filename for all output")
    cmdline.add_argument("-r", "--rerun", dest='rerun', default=False, action='store_true',
                         help="rerun using profile generated in prior run")
    cmdline.add_argument("--srclist", dest='srclist', default=None,
                         help="list of source-ids to include in valid fitting area")
    #
    # optional parameters affecting ellipse fitting and profile generation
    #
    cmdline.add_argument("--minsma", dest='minsma', default=0, type=float,
                         help="minimum semi-major axis length (default: 0.0 pixels)")
    cmdline.add_argument("--maxsma", dest='maxsma', default=None, type=float,
                         help="maximum semi-major axis length (default: no maximum)")

    args = cmdline.parse_args()

    fn = args.input_fn
    catalog_fn = args.catalog_fn
    source_id = args.source_id #[int(id) for id in args.source_id.split(",")]
    segm_fn = args.segmentation_fn
    output_basename = args.output_basename

    # read image
    img_hdu = pyfits.open(fn)
    img = img_hdu[0].data
    img_masked = numpy.ma.array(img)

    # read segmentation frame
    segm_hdu = pyfits.open(segm_fn)
    segm = segm_hdu[0].data
    use_for_fit = ((segm == 0) | (segm == source_id)) & numpy.isfinite(img)
    if (args.srclist is not None):
        srclist = [int(id) for id in args.srclist.split(",")]
        for id in srclist:
            use_for_fit |= (segm == id)

    img_masked.mask = ~use_for_fit

    # read catalog
    vot = votable.parse_single_table(catalog_fn)
    catalog = vot.to_table()
    # print(catalog)
    # print(source_id)
    right_source = (catalog['NUMBER'] == source_id)
    source_data = catalog[right_source][0]
    print("Assuming SOURCE-coord: ", source_data['X_IMAGE'], source_data['Y_IMAGE'])

    # sys.exit(-1)
    n_attempts = 1

    rerun = args.rerun
    surfprofile_csv = output_basename + "_profile.csv"

    if (rerun and os.path.isfile(surfprofile_csv)):
        print("Skipping profile generation, reading instead from %s" % (surfprofile_csv))
        df = pandas.read_csv(surfprofile_csv)
    else:
        print("Preparing ellipse fit")
        dummy_fn = output_basename+"_maskcheck.fits"
        pyfits.PrimaryHDU(data=img_masked.filled(0)).writeto(dummy_fn, overwrite=True)
        print("wrote mask-check image to %s" % (dummy_fn))

        # prepare ellipse fitting
        geo = photutils.isophote.EllipseGeometry(
            x0 = source_data['X_IMAGE'],
            y0 = source_data['Y_IMAGE'],
            sma=10, eps=0.01, pa=0,
        )
        geo.find_center(img_masked)
        # print(geo)
        ellipse = photutils.isophote.Ellipse(img_masked, geometry=geo)
        # print("ellipse:", ellipse)

        sma_start = args.minsma
        for a in range(1): #n_attempts):
            # repeat this until we run out of attmepts or get a complete profile

            # df = pandas.read_csv("surfprofile.csv")
            print("Fitting image (min_sma=%.2f, max_sma=%s)" % (sma_start, "%.2f" % args.maxsma if args.maxsma is not None else "NO_LIMIT"))
            isophot_list = ellipse.fit_image(
                # integrmode='median',
                sclip=3.0, nclip=3, fflag=0.9,
                sma0=sma_start,
                maxsma=args.maxsma,
                minsma=args.minsma,
                step=0.1,
            )
            print("done fitting")
            print(isophot_list.to_table())
            df = isophot_list.to_table().to_pandas()
            df.info()

            df.to_csv(surfprofile_csv+"%d"%a)
            print("Saved surface brightness profile as CSV (%s)" % (surfprofile_csv))

            # check if any of the values are NaN - indicating a problem
            bad_data = df.isnull().any(axis=1)
            print(bad_data)
            if (len(bad_data.index) > 0):
                # now find the index of the first bad entry
                first_bad = numpy.min(df.index[bad_data])
                print(first_bad)
                print(df.iloc[first_bad])

                sma_start = df.iloc[first_bad]['sma'] * 1.1
                print("Restarting profile generation at radius sma=%f" % (sma_start))

            else:
                print("Profile generation completed without issues")
                break

            # break

    #
    # convert all ellipses into ds9 format
    #
    fit_result_colors = {
        0: 'color=green',
        1: 'color=yellow',
        2: 'color=orange',
        3: 'color=orange',
        4: 'color=blue',
        5: 'color=red width=3',
        -1: ''
    }

    ds9_reg_fn = output_basename + "_ellipses.reg"
    with open(ds9_reg_fn, "w") as ds9_reg:
        print("""\
# Region file format: DS9 version 4.1
global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
image""", file=ds9_reg)
        for index, e in df.iterrows():
            # print(e)
            print("ellipse(%f,%f,%f,%f,%f) # %s" % (
                e['x0'], e['y0'],
                e['sma'], (1.-e['ellipticity'])*e['sma'],
                e['pa'],
                fit_result_colors[e['stop_code']]
                ),
                  file=ds9_reg)

    # now write the region file again, this time using WCS instead of pixel information
    hdr = img_hdu[0].header
    try:
        pixelscale = numpy.hypot(hdr['CD1_1'], hdr['CD1_2']) * 3600. # in arcsec/pixel
    except:
        try:
            pixelscale = hdr['CDELT1'] * 3600
        except:
            pixelscale = 1.
            pass

    # convert all x/y into ra/dec
    wcs = astropy.wcs.WCS(hdr)
    (ra0,dec0) = wcs.all_pix2world(df['x0'], df['y0'], 0)
    # print(ra_dec)
    df['ra0'] = ra0 #_dec[:,0]
    df['dec0'] = dec0 #ra_dec[:,1]
    ds9_reg_fn = output_basename + "_ellipses_wcs.reg"
    with open(ds9_reg_fn, "w") as ds9_reg:
        print("""\
# Region file format: DS9 version 4.1
global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
fk5""", file=ds9_reg)
        for index, e in df.iterrows():
            # print(e)
            print('ellipse(%f,%f,%f",%f",%f) # %s' % (
                e['ra0'], e['dec0'],
                e['sma']*pixelscale, (1.-e['ellipticity'])*e['sma']*pixelscale,
                e['pa'],
                fit_result_colors[e['stop_code']]
            ),
                  file=ds9_reg)


