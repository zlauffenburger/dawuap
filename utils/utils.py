from __future__ import division
import numpy as np
import shapefile as shp
import fiona
from shapely.geometry import shape
import rasterio as rio
import pickle

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


def write_array_as_tiff(fn, template, array):
    """
    Utility ot write numpy arrays as tiff file using metadata from a template map
    :param fn: filename of output tiff file, string
    :param template: filename of valid image to obtain metadata for output image
    :param array: numpy array to be converted to tiff
    :return: None
    """
    with rio.open(template, 'r') as src:
        profile = src.profile
        if not isinstance(array, np.ndarray) or array.shape != src.shape:
            raise ValueError("Shape mismatch!. Template shape is %s. Arrays shape is %s" % (src.shape, array.shape))
        profile.update(driver='GTiff', count=1, dtype='float64')
        with rio.open(fn, 'w', **profile) as dst:
            dst.write(array.astype('float64'), 1)


def write_structured_parameter_array(filenames, shape):
    tp = zip(filenames.keys(), ['float32' for i in range(len(filenames))], [shape for i in range(len(filenames))])
    pars = np.empty(len(filenames), dtype=tp)

    for name, items in filenames.items():
        with rio.open(items, 'r') as src:
            pars[name] = src.read()



    #np.zeros(filenames.items(), dtype=[('', '', )])

def add_rr_model_parameters_to_shapefile(shapefile, outshp=''):
    """
    Add default field and parameter values to shapefile. A bug in the pyShp library 
    handles incorrectly Date fields. If date fields are contained in the shapefile 
    the library needs to be patched including the follwing lines in shapefile.py:

    from datetime import date
    comment out line 525 and replace by:
     value = date(y,m,d).strftime('%Y%m%d')
    :param shapefile: Polygon (catchment) shapefile to which parameter fields will be added
    :param outshp: optional, output shapefile filename. Overwrite original if empy 
    :return: None 
    """
    r = shp.Reader(shapefile)
    w = shp.Writer(r.shapeType)
    w.autoBalance = 1

    w.fields = list(r.fields)
    #w.records.extend(r.records())
    w._shapes.extend(r.shapes())

    fldnames = [fld[0] for fld in w.fields]

    if 'hbv_ck0' not in fldnames:
        w.field('hbv_ck0', 'N', 10, 6)
    if 'hbv_ck1' not in fldnames:
        w.field('hbv_ck1', 'N', 10, 6)
    if 'hbv_ck2' not in fldnames:
        w.field('hbv_ck2', 'N', 10, 6)
    if 'hbv_hl1' not in fldnames:
        w.field('hbv_hl1', 'N', 10, 6)
    if 'hbv_perc' not in fldnames:
        w.field('hbv_perc', 'N', 10, 6)
    if 'hbv_pbase' not in fldnames:
        w.field('hbv_pbase', 'N', 10, 0)

    counter = 1
    for ca in r.records():
        print "Appending parameter to record ", counter
        ca.append(10.0)
        ca.append(50.0)
        ca.append(10000.0)
        ca.append(50.0)
        ca.append(50.0)
        ca.append(5)
        w.records.append(ca)
        counter += 1

    #w.records.append(('hbv_ck0', '1.0'))
    if not outshp:
        outshp=shapefile

    w.save(outshp)
    r = w = None

#add_rr_model_parameters_to_shapefile('test_data/HUC8_NetworkLite.shp')

