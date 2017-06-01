import argparse
import utilsVector


def main(argc):
    """
    This script writes default rainfall-runoff model parameters to a polygon vector file.
     
    :param vectorFile: 'vector polygon file with model catchments'
    :param outFile:  'optional output filename with modified vectorFile'
    """

    utilsVector.add_rr_model_parameters_to_shapefile(argc.vectorFile, argc.outFile)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Writes rr model parameters to polygon shapefile')
    parser.add_argument('vectorFile')
    parser.add_argument('outFile', nargs='?', default='')
    args = parser.parse_args()
    main(args)
