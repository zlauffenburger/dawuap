from __future__ import division
import numpy as np
import pandas as pd
import shapefile as shp
import fiona
from shapely.geometry import mapping, shape
import rasterio as rio


class ReadVector(object):

    def __init__(self, fn_vector):
        self.fn_vector = fn_vector
        # Next three variables are updated after the call to calc_connectivity_matrix()
        self.crs = None
        self.schema = None
        self.driver = None
        self._read_features()

    def _read_features(self):
        # test it as fiona data source
        features_iter = None
        try:
            with fiona.open(self.fn_vector, 'r') as src:
                assert len(src) > 0
                self.crs = src.crs
                self.schema = src.schema.copy()
                self.driver = src.driver

            def fiona_generator(obj):
                with fiona.open(obj, 'r') as src2:
                    for feature in src2:
                        yield feature

            features_iter = fiona_generator(self.fn_vector)

        except (AssertionError, TypeError, IOError, OSError):
            print "fn_vector does not point to a fiona object, error "

        return features_iter


class ParseNetwork(ReadVector):

    def __init__(self, fn_vector):
        super(ParseNetwork, self).__init__(fn_vector)
        self.conn_matrix = self.calc_connectivity_matrix()

    def calc_connectivity_matrix(self):
        feature_iter = self._read_features()

        lst_node_connections = []

        for i, feats in enumerate(feature_iter):
            geom = shape(feats['geometry'])
            props = feats['properties']
            lst_node_connections.append((props["FROM_NODE"], props["TO_NODE"]))

        lstNodes = set(zip(*lst_node_connections)[0])
        df = pd.DataFrame(0, index=lstNodes, columns=lstNodes)
        # drop the connections that go out of the basin to node 0
        lst_node_connections = [i for i in lst_node_connections if i[1]>0]
        for link in lst_node_connections:
            df.loc[link] = 1

        return df


class ParameterIO(ReadVector):
    def __init__(self, fn_vector):
        super(ParameterIO, self).__init__(fn_vector)

    def _write_fiona_object(self, fn, crs=None, driver=None, schema=None, params=None):
        """
        Writes a vector file using fn_vector as template
        :param fn: outfile name
        :param params: dictionary
        :return: None
        """
        if crs is None:
            crs = self.crs
        if driver is None:
            driver = self.driver
        if schema is None:
            schema = self.schema

        # TODO: solve problems writing projection information in shapefiles and geojson
        feature_iter = self._read_features()
        with fiona.open(fn, 'w', crs=crs, driver=driver, schema=schema) as sink:

            for i, feats in enumerate(feature_iter):
                geom = shape(feats['geometry'])
                if params is not None:
                    props = params[i]['properties']
                else:
                    props = feats['properties']
                sink.write({'properties': props, 'geometry': mapping(geom)})

    def write_muskingum_parameters(self, outfn, params):
        # type: (str, list) -> None
        """
        Adds or updates the vector network file with the Muskingum-Cunge
        parameters.
        :param outfn: filename for updated vector network
        :param params: list of parameter dictionaries with format [{'ARCID': ID, 'e': value, 'ks': value},{}]
        :return: None
        """
        schema = self.schema.copy()
        schema['properties']['e'] = 'float'
        schema['properties']['ks'] = 'float'

        lstDicts = []
        feature_iter = self._read_features()
        for i, feats in enumerate(feature_iter):
            arc_id = feats['properties']['ARCID']
            try:
                val = (item for item in params if item['ARCID'] == arc_id).next()
            except:
                val = {}
            feats['properties']['e'] = val.get('e', 0.35)
            feats['properties']['ks'] = val.get('ks', 82400)
            lstDicts.append(feats)

        self._write_fiona_object(outfn, schema=schema, params=lstDicts)

    def write_hvb_parameters(self, outfn, params):
        # type: (str, list) -> None
        """
        Adds or updates the vector file of model subwatersheds with the hbv RR model parameters.

        :param outfn: filename for updated vector network
        :param params: list of parameter dictionaries with format [{'ARCID': ID, 'e': value, 'ks': value},{}]
        :return: None
        """
        schema = self.schema.copy()
        schema['properties']['hbv_ck0'] = 'float'
        schema['properties']['hbv_ck1'] = 'float'
        schema['properties']['hbv_ck2'] = 'float'
        schema['properties']['hbv_hl1'] = 'float'
        schema['properties']['hbv_perc'] = 'float'
        schema['properties']['hbv_pbase'] = 'int'

        lstDicts = []
        feature_iter = self._read_features()
        for i, feats in enumerate(feature_iter):
            arc_id = feats['properties']['GRIDCODE']
            try:
                val = (item for item in params if item['GRIDCODE'] == arc_id).next()
            except:
                val = {}

            feats['properties']['hbv_ck0'] = val.get('hbv_ck0', 10.)
            feats['properties']['hbv_ck1'] = val.get('hbv_ck1', 50.)
            feats['properties']['hbv_ck2'] = val.get('hbv_ck2', 10000)
            feats['properties']['hbv_hl1'] = val.get('hbv_hl1', 50)
            feats['properties']['hbv_perc'] = val.get('hbv_perc', 50)
            feats['properties']['hbv_pbase'] = val.get('hbv_pbase', 5)
            lstDicts.append(feats)

        self._write_fiona_object(outfn, schema=schema, params=lstDicts)






def write_array_as_tiff(fn, template, array):
    """
    Utility to write numpy arrays as tiff file using metadata from a template map
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


def add_muskingum_model_parameters_to_network(vectorNetwork, outshp=''):
    """
    Add default muskingum cunge parameters to streamflow network
    :param vectorNetwork:
    :param outshp:
    :return:
    """

    def _read_features(inFile):
        # test it as fiona data source
        features_iter = None
        try:
            with fiona.open(inFile, 'r') as src:
                assert len(src) > 0

            def fiona_generator(obj):
                with fiona.open(obj, 'r') as src:
                    for feature in src:
                        yield feature

            features_iter = fiona_generator(inFile)

        except (AssertionError, TypeError, IOError, OSError):
            print "fn_vector does not point to a fiona object, error "

        return features_iter

    feature_iter = _read_features(vectorNetwork)
    lst_node_connections = []

    for i, feats in enumerate(feature_iter):
        geom = shape(feats['geometry'])
        props = feats['properties']
        lst_node_connections.append((props["FROM_NODE"], props["TO_NODE"]))


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

