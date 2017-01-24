from __future__ import division
import numpy as np
import shapefile as shp
import rasterio as rio
import pickle


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

#add_rr_model_parameters_to_shapefile('test_data/HUC8_NetworkLite.shp')

