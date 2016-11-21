try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
        'description': 'Hydrologic Engine for daWUAP model',
        'author': 'Marco Maneta',
        'url': 'http://www.',
        'download_url': 'http://www.',
        'author_email': 'marco.maneta@umontana.edu',
        'version': '0.1',
        'install_requires': ['numpy', 'nose', 'rasterstats', 'rasterio'],
        'packages': ['hydroengine'],
        'scripts': [],
        'name': 'daWUAP'
}

setup(**config)
