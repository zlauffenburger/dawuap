from __future__ import division
import numpy as np
import shapefile as shp


def add_rr_model_parameters_to_shapefile(shapefile):
    r = shp.Reader(shapefile)
    w = shp.Writer(r.shapeType)
    w.autoBalance = 1

    w.fields = list(r.fields)
    #w.records.extend(r.records())
    w._shapes.extend(r.shapes())

    if 'hbv_ck0' not in w.fields:
        w.field('hbv_ck0', 'C', 40)
    if 'hbv_ck1' not in w.fields:
        w.field('hbv_ck1', 'N', 5, 5)
    if 'hbv_ck3' not in w.fields:
        w.field('hbv_ck2', 'N', 5, 5)

    for ca in r.records():
        print "writing record ", ca
        ca.append(5.4)
        ca.append()
        w.records.append(ca)

    #w.records.append(('hbv_ck0', '1.0'))

    w.save('/Users/marcomaneta/Documents/DataSandBox/MontanaHydrology/test')

add_rr_model_parameters_to_shapefile('/Users/marcomaneta/Documents/DataSandBox/MontanaHydrology/WBDHU8_MT.shp')

