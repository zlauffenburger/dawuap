from __future__ import division
import numpy as np
import shapefile as shp
import fiona
from shapely.geometry import shape


class ParseNetwork(object):

    def __init__(self, fn_vector):
        self.fn_vector = fn_vector


    def _read_features(self):
        # test it as fiona data source
        features_iter = None
        try:
            with fiona.open(self.fn_vector, 'r') as src:
                assert len(src) > 0

            def fiona_generator(obj):
                with fiona.open(obj, 'r') as src:
                    for feature in src:
                        yield feature

            features_iter = fiona_generator(self.fn_vector)

        except (AssertionError, TypeError, IOError, OSError):
            print "fn_vector does not point to a fiona object, error ", e

        return features_iter

    def parse_network(self):
        feature_iter = self._read_features()

        for i, feats in enumerate(feature_iter):
            geom = shape(feats['geometry'])
            props = shape(feats['properties'])


a = ParseNetwork('test_data/mt_network.geojson')
a.parse_network()


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

