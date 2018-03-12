
import shapely
shapely.__file__
import fiona
import sys

sys.path.append('../../alg')
import distance_buffer as db

shp = fiona.open("../data/YukonAlaska/YukonRegionCRS.shp")
print shp.schema
l = iter(shp)

first = next(l)
shp_geom = shapely.geometry.shape(first['geometry'])

region = shp_geom
print len(region.exterior.coords)
print len(region.interiors)
print region.area
print type(region)
regionring = region.exterior

buffer1 = regionring.buffer(15.0)
print type(buffer1)
print type(region)
print len(buffer1.interiors)
print len(region.exterior.coords)