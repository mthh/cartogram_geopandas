# -*- coding: utf-8 -*-
"""
cartogram_geopandas
"""

from distutils.core import setup
from Cython.Build import cythonize

setup(
    name='cartogram_geopandas',
    version='0.0.0c',
    description='Fast and convenient cartogram creation on a GeoDataFrame',
    author='mthh',
    py_modules=['cartogram_geopandas'],
    ext_modules=cythonize("cycartogram.pyx"),
    license='GPL v2',
    )
