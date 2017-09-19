=============
Usage example
=============

1. open a terminal and ``cd`` to the example directory
2. Retrieve climate data from Idaho's Metdata Thredds server from Sept 1, 2012 to August 31, 2013::

    dataCollectionThredds.py -ds 2012-09-01 -de 2013-08-31 -a precip --flip  vectorFile Domain_latlong.shp
    dataCollectionThredds.py -ds 2012-09-01 -de 2013-08-31 -a tempmax --flip  vectorFile Domain_latlong.shp
    dataCollectionThredds.py -ds 2012-09-01 -de 2013-08-31 -a tempmin --flip  vectorFile Domain_latlong.shp

3. Write default parameters for the HBV model and the Muskingum routing model on the basin shapefiles::

    writeVectorModelParameters.py -subBasinsVectorFile subs1.shp -subBasinOutFile subsout.shp
    -streamVectorFile riv1.shp -streamOutFile rivout.shp
    
4. Write default parameters for the raster components of the rainfall runoff model::

    writeRasterModelParameters.py precip_F2012-09-01_T2013-08-31.nc

5. Run hydrologic model::

    hydrovehicle.py 08/31/2016 precip_F2012-09-01_T2013-08-31.nc tempmin_F2012-09-01_T2013-08-31.nc
     tempmax_F2012-09-01_T2013-08-31.nc param_files_test.json rivout.shp subsout.shp

