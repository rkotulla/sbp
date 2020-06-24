#!/usr/bin/env python3


import os
import sys
import astropy.io.fits as pyfits
import numpy
import pandas
import photutils
import photutils.isophote
from astropy.io import votable

print(photutils.__path__)

if __name__ == "__main__":

    fn = sys.argv[1]
    catalog_fn = sys.argv[2]
    source_id = int(sys.argv[3])
    segm_fn = sys.argv[4]
    output_basename = sys.argv[5]

    # read image
    img_hdu = pyfits.open(fn)
    img = img_hdu[0].data
    img_masked = numpy.ma.array(img)

    # read segmentation frame
    segm_hdu = pyfits.open(segm_fn)
    segm = segm_hdu[0].data
    use_for_fit = ((segm == 0) | (segm == source_id)) & numpy.isfinite(img)
    img_masked.mask = ~use_for_fit

    # read catalog
    vot = votable.parse_single_table(catalog_fn)
    catalog = vot.to_table()
    print(catalog)
    print(source_id)
    right_source = (catalog['NUMBER'] == source_id)
    source_data = catalog[right_source][0]
    print("SOURCE-coord: ", source_data['X_IMAGE'], source_data['Y_IMAGE'])

    # sys.exit(-1)

    rerun = True
    surfprofile_csv = output_basename + "_profile.csv"

    if (rerun and os.path.isfile(surfprofile_csv)):
        df = pandas.read_csv(surfprofile_csv)
    else:
        print("Preparing ellipse fit")
        dummy_fn = output_basename+"_maskcheck.fits"
        pyfits.PrimaryHDU(data=img_masked.filled(0)).writeto(dummy_fn, overwrite=True)

        # prepare ellipse fitting
        geo = photutils.isophote.EllipseGeometry(
            x0 = source_data['X_IMAGE'],
            y0 = source_data['Y_IMAGE'],
            sma=10, eps=0.01, pa=0,
        )
        # geo.find_center(img_masked)
        print(geo)
        ellipse = photutils.isophote.Ellipse(img_masked, geometry=geo)
        print("ellipse:", ellipse)

        # df = pandas.read_csv("surfprofile.csv")
        print("Fitting image")
        isophot_list = ellipse.fit_image(
            # integrmode='median',
            sclip=3.0, nclip=3, fflag=0.7,
            maxsma=150,
            minsma=5,
        )
        print("done fitting")
        print(isophot_list.to_table())

        df = isophot_list.to_table().to_pandas()
        df.to_csv(surfprofile_csv)
        print("Saved surface brightness profile as CSV")


    #
    # convert all ellipses into ds9 format
    #
    ds9_reg_fn = output_basename + "_ellipses.reg"
    with open(ds9_reg_fn, "w") as ds9_reg:
        print("""\
# Region file format: DS9 version 4.1
global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
image""", file=ds9_reg)
        for index, e in df.iterrows():
            # print(e)
            print("ellipse(%f,%f,%f,%f,%f)" % (e['x0'], e['y0'], e['sma'], (1.-e['ellipticity'])*e['sma'], e['pa']),
                  file=ds9_reg)


