import argparse
import utils
import os, sys

def main(argc):
    """This script writes default rainfall-runoff model parameters to a polygon vector file.
    """
    basinFile = utils.ModelVectorDatasets(fn_subsheds=argc.in_basins, fn_network=argc.in_rivers)
    basinFile.write_hvb_parameters(argc.out_basins, argc.basin_params)
    basinFile.write_muskingum_parameters(argc.out_rivers, argc.routing_params)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Writes rainfall-runoff model parameters to polygon and polyline vector'
                                                 ' file as new attributes.',
                                     epilog="""Default operation is to read a Fiona vector files provided in the 
                                     input argument and add attributes to each feature of the file. Each feature  
                                     will contain a default value. 
                                     
                                     Example: {} -subBasinVectorFile ./data/subs.shp""".format(os.path.basename(sys.argv[0])))
    parser.add_argument('-subBasinsVectorFile', dest='in_basins', type=str, default=None, help="Path to sub-basin polygon file to be modified. Any format supported by Fiona.")
    parser.add_argument('-streamVectorFile', dest='in_rivers', type=str, default=None,
                        help="Path to stream network vector file to be modified. Any format supported by Fiona.")
    parser.add_argument('-subBasinOutFile', dest='out_basins',type=str, default='', help="(Optional) Path to modified subBasinsVector file. If none is given the"
                                                               " input vector file is overwritten.")
    parser.add_argument('-streamOutFile', dest='out_rivers', type=str, default='',
                        help="(Optional) Path to modified subBasinsVector file. If none is given the"
                             " input vector file is overwritten.")
    parser.add_argument('-basinParams', dest='basin_params', type=str, default=None, help="(Optional) path to list of python "
                                                                               "dictionaries containing the parameter "
                                                                               "values for each polygon in subBasinVectorFile."
                                                                               " If none is given, default parameters"
                                                                               "are written. See Documentation for"
                                                                               " further information.")
    parser.add_argument('-routingParams', dest='routing_params', type=str, default=None,
                        help="(Optional) path to list of python "
                             "dictionaries containing routing parameter "
                             "values for each polygon in streamVectorFile."
                             " If none is given, default parameters"
                             "are written. See Documentation for"
                             " further information.")
    args = parser.parse_args()
    main(args)
