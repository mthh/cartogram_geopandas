# -*- coding: utf-8 -*-
"""
 cartogram_geopandas v0.0.0b:

    Easy construction of continuous cartogram on a Polygon/MultiPolygon
    GeoDataFrame (modify the geometry in place or create a new GeoDataFrame).

    Code adapted from Carson Farmer code
    (https://github.com/carsonfarmer/cartogram : code used in 'Cartogram' QGis
    python plugin before some changes were made) to fit the
    geopandas.GeoDataFrame datastructure. Carson Farmer's code is partially
    related to 'pyCartogram.py' from Eric Wolfs.

    Algorithm itself based on
        { Dougenik, J. A, N. R. Chrisman, and D. R. Niemeyer. 1985.
          "An algorithm to construct continuous cartograms."
          Professional Geographer 37:75-81 }

    No warranty concerning the result.
    Copyright (C) 2013 Carson Farmer, 2015  mthh
"""

import math
from shapely.geometry import (
    Point, LineString, MultiLineString, Polygon, MultiPolygon
    )


def make_cartogram(geodf, field_name, iterations=5, inplace=False):
    crtgm = Cartogram(geodf, field_name, iterations, inplace=inplace)
    crtgm.make()
    if not inplace:
        return crtgm.result


def transform_geom(aLocal, dForceReductionFactor,
                   geom, featCount, sqrt=math.sqrt):
    new_geom = []
    if isinstance(geom, Polygon):
        geom = [geom]
    for single_geom in geom:
        boundarys = single_geom.boundary
        tmp_bound = []
        if not isinstance(boundarys, MultiLineString):
            boundarys = [boundarys]
        for single_boundary in boundarys:
            line_coord = []
            line_add_pt = line_coord.append
            for x, y in zip(single_boundary.coords.xy[0],
                            single_boundary.coords.xy[1]):
                x0, y0 = x, y
                # Compute the influence of all shapes on this point
                for i in range(featCount):
                    lf = aLocal[i]
                    cx = lf.ptCenter_x
                    cy = lf.ptCenter_y
                    # Pythagorean distance
                    distance = sqrt((x0 - cx) ** 2 + (y0 - cy) ** 2)

                    if (distance > lf.dRadius):
                        # Calculate the force on verteces far away
                        # from the centroid of this feature
                        Fij = lf.dMass * lf.dRadius / distance
                    else:
                        # Calculate the force on verteces far away
                        # from the centroid of this feature
                        xF = distance / lf.dRadius
                        Fij = lf.dMass * (xF ** 2) * (4 - (3 * xF))
                    Fij = Fij * dForceReductionFactor / distance
                    x = (x0 - cx) * Fij + x
                    y = (y0 - cy) * Fij + y
                line_add_pt(Point(x, y))
            line = LineString(line_coord)
            tmp_bound.append(line)

        if len(tmp_bound) == 1:
            poly = Polygon(tmp_bound[0])
        else:
            poly = MultiPolygon(
                [Polygon(sl) for sl in MultiLineString(tmp_bound)]
                )
        new_geom.append(poly)

    if len(new_geom) > 1:
        return MultiPolygon(new_geom)
    elif len(new_geom) == 1:
        return new_geom[0]


class Holder(object):
    def __init__(self):
        self.lFID = -1
        self.ptCenter_x = -1
        self.ptCenter_y = -1
        self.dValue = -1
        self.dArea = -1
        self.dMass = -1
        self.dRadius = -1


class Cartogram(object):
    def __init__(self, geodf, field_name, iterations, inplace=False):
        allowed = {'MultiPolygon', 'Polygon'}
        geom_type = {i for i in geodf.geom_type}
        if not geom_type.issubset(allowed):
            raise ValueError(
                "Geometry type doesn't match 'Polygon'/'MultiPolygon"
                )
        if inplace:
            self.geodf = geodf
        else:
            self.geodf = geodf.copy()
        self.iterations = iterations
        self.index_field = [i for i, j in enumerate(list(self.geodf.columns))
                            if field_name in j]
        self.result = None
        self.total_features = len(self.geodf)
        self.pi = math.pi
        self.sqrt = math.sqrt

    def make(self):
        res_geom, it = self.cartogram()
        assert it == self.iterations
        self.geodf.set_geometry(res_geom, inplace=True)
        self.result = self.geodf

    def cartogram(self):
        iterations = self.iterations
        temp_geo_serie = self.geodf.geometry.copy()
        for i in range(iterations):
            # Useless on the first iteration :
            (aLocal, dForceReductionFactor) = self.getInfo(self.index_field[0])
            new_geo_serie = temp_geo_serie.copy()

            for fid, geom in temp_geo_serie.items():
                newgeom = transform_geom(aLocal, dForceReductionFactor,
                                         geom, self.total_features)
                new_geo_serie[fid] = newgeom
            temp_geo_serie = new_geo_serie
        return temp_geo_serie, i+1

    # Gets the information required for calcualting size reduction factor
    def getInfo(self, index):
        sqrt = self.sqrt
        pi = self.pi
        featCount = self.total_features
        aLocal = []
        cx = 0
        cy = 0
        area_total = self.geodf.area.sum()
        value_total = self.geodf.iloc[:, index].sum()
        for fid, geom in self.geodf.geometry.items():
            lfeat = Holder()
            lfeat.dArea = geom.area  # save area of this feature
            lfeat.lFID = fid  # save id for this feature
            lfeat.dValue = self.geodf.iloc[fid, index]  # save 'area' value for this feature
            wkt_centroid = geom.centroid.wkt
            cx, cy = wkt_centroid.replace(
                'POINT (', ''
                ).replace(')', '').split(' ')
            lfeat.ptCenter_x = float(cx)  # save centroid x for this feature
            lfeat.ptCenter_y = float(cy)  # save centroid y for this feature
            aLocal.append(lfeat)

        dFraction = area_total / value_total
        dSizeErrorTotal = 0
        dSizeError = 0
        for i in range(featCount):
            lf = aLocal[i]  # info for current feature
            dPolygonValue = lf.dValue
            dPolygonArea = lf.dArea
            if (dPolygonArea < 0):  # area should never be less than zero
                dPolygonArea = 0
            # this is our 'desired' area...
            dDesired = dPolygonValue * dFraction
            # calculate radius, a zero area is zero radius
            dRadius = sqrt(dPolygonArea / pi)
            lf.dRadius = dRadius
            tmp = dDesired / pi
            if tmp > 0:
                # calculate area mass, don't think this should be negative
                lf.dMass = sqrt(dDesired / pi) - dRadius
            else:
                lf.dMass = 0
            # both radius and mass are being added to the feature list for
            # later on...
            # calculate size error...
            dSizeError = \
                max(dPolygonArea, dDesired) / min(dPolygonArea, dDesired)
            # this is the total size error for all polygons
            dSizeErrorTotal = dSizeErrorTotal + dSizeError
        # average error
        dMean = dSizeErrorTotal / featCount
        # need to read up more on why this is done
        dForceReductionFactor = 1 / (dMean + 1)

        return (aLocal, dForceReductionFactor)
