import argparse
import utils


def main(argc):
    """This script writes default rainfall-runoff model parameters to a polygon vector file.

    :param vectorFile: 'vector polygon file with model catchments'
    :param outFile:  'optional output filename with modified vectorFile'
    :param --params: 'optional dictionary with parameter values
    """
    vecFile = utils.ModelVectorDatasets(fn_subsheds=argc.vectorFile)
    vecFile.write_hvb_parameters(argc.outFile, argc.params)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Writes rr model parameters to polygon shapefile')
    parser.add_argument('vectorFile')
    parser.add_argument('outFile', nargs='?', default='')
    parser.add_argument('-params', dest='params', type=str, default=None)
    args = parser.parse_args()
    main(args)
