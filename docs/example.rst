=============
Usage example
=============

1. open a terminal and ``cd`` to the example directory
2. Retrieve climate data from Idaho's Metdata Thredds server from Sept 1, 2012 to August 31, 2013::

    dataCollectionThredds.py -ds 2012-09-01 -de 2013-08-31 -a precip --flip  vectorFile Domain_latlong.shp
    dataCollectionThredds.py -ds 2012-09-01 -de 2013-08-31 -a tempmax --flip  vectorFile Domain_latlong.shp
    dataCollectionThredds.py -ds 2012-09-01 -de 2013-08-31 -a tempmin --flip  vectorFile Domain_latlong.shp

3. 