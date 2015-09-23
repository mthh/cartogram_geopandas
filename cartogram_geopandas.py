# -*- coding: utf-8 -*-
"""
 cartogram_geopandas v0.0.0c:

    Easy construction of continuous cartogram on a Polygon/MultiPolygon
    GeoDataFrame (modify the geometry in place or create a new GeoDataFrame).

    Code adapted to fit the geopandas.GeoDataFrame datastructure from
    Carson Farmer's code (https://github.com/carsonfarmer/cartogram : former
    code in use in 'Cartogram' QGis python plugin). Carson Farmer's code is
    partially related to 'pyCartogram.py' from Eric Wolfs.

    Algorithm itself based on
        { Dougenik, J. A, N. R. Chrisman, and D. R. Niemeyer. 1985.
          "An algorithm to construct continuous cartograms."
          Professional Geographer 37:75-81 }

    No warranty concerning the result.
    Copyright (C) 2013 Carson Farmer, 2015  mthh
"""
from cycartogram import Cartogram


def make_cartogram(geodf, field_name, iterations=5, inplace=False):
    """
    Make a continuous cartogram on a geopandas.GeoDataFrame collection
    of Polygon/MultiPolygon (wrapper to call the core functions
    written in cython).

    :param geopandas.GeoDataFrame geodf: The GeoDataFrame containing the
        geometry and a field to use for the transformation.

    :param string field_name: The name of the field (Series) containing the
        value to use.

    :param integer iterations: The number of iterations to make.
        [default=5]

    :param bool inplace: Append in place if True is set. Otherwhise return a
        new GeoDataFrame with transformed geometry.
        [default=False]
    """
    if inplace:
        crtgm = Cartogram(geodf, field_name, iterations)
        crtgm.make()
    else:
        crtgm = Cartogram(geodf.copy(), field_name, iterations)
        return crtgm.make()
