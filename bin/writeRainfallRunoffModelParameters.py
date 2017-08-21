import argparse
import utils
import os, sys

def main(argc):
    """This script writes default rainfall-runoff model parameters to a polygon vector file.
    """

    vecFile = utils.ModelVectorDatasets(fn_subsheds=argc.vectorFile)
    vecFile.write_hvb_parameters(argc.outFile, argc.params)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Writes rainfall-runoff model parameters to polygon vector file as'
                                                 ' new attributes.',
                                     epilog="""Default operation is to read a Fiona vector file provided in the 
                                     mandatory argument and add attributes to each polygon of the file with 
                                     default values. 
                                     
                                     Example: {} ./data/subs.shp""".format(os.path.basename(sys.argv[0])))
    parser.add_argument('vectorFile', help="Path to input vector file to be modified. Any format supported by Fiona.")
    parser.add_argument('outFile', nargs='?', default='', help="(Optional) Path to output vector file. If none is given the"
                                                               " input vector file is overwritten.")
    parser.add_argument('-params', dest='params', type=str, default=None, help="(Optional) path to list of python "
                                                                               "dictionaries containing the parameter "
                                                                               "values for each polygon in vectorFile."
                                                                               " If none is given, default parameters"
                                                                               "are written. See Documentation for"
                                                                               " further information.")
    args = parser.parse_args()
    main(args)
