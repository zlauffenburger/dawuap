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
        w.field('hbv_ck0', 'N', 10, 6)
    if 'hbv_ck1' not in w.fields:
        w.field('hbv_ck1', 'N', 10, 6)
    if 'hbv_ck3' not in w.fields:
        w.field('hbv_ck2', 'N', 10, 6)
    if 'hbv_hl1' not in w.fields:
        w.field('hbv_hl1', 'N', 10, 6)
    if 'hbv_perc' not in w.fields:
        w.field('hbv_perc', 'N', 10, 6)

    counter = 1
    for ca in r.records():
        print "Appending parameter to record ", counter
        ca.append(10.0)
        ca.append(50.0)
        ca.append(10000.0)
        ca.append(50.0)
        ca.append(50.0)
        w.records.append(ca)
        counter += 1

    #w.records.append(('hbv_ck0', '1.0'))

    w.save('/Users/marcomaneta/Documents/DataSandBox/MontanaHydrology/test')

add_rr_model_parameters_to_shapefile('/Users/marcomaneta/Documents/workspace/linknodemodel/tests/test_data/HUC8_NetworkLite.shp')

