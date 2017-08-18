try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
        'description': 'Hydrologic Engine for daWUAP model',
        'author': 'Marco Maneta',
        'url': 'https://bitbucket.org/umthydromodeling/dawuaphydroengine',
        'download_url': 'https://bitbucket.org/umthydromodeling/dawuaphydroengine/get/a88ac8f68b01.zip',
        'author_email': 'marco.maneta@umontana.edu',
        'version': '0.1',
        'install_requires': ['numpy', 'nose', 'rasterstats', 'rasterio', 'fiona', 'shapely', 'cdo'],
        'packages': ['hydroengine', 'utils', 'hydrovehicle'],
        'scripts': ['bin/hydrovehicle.py', 'bin/dataCollectionThredds.py', 'bin/writeRainfallRunoffModelParameters.py'],
        'name': 'daWUAP'
}

setup(**config)
