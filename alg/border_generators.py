import shapely
from shapely.geometry import Point, Polygon, MultiPolygon
import fiona
import osmnx as ox, geopandas as gpd

#the following code block comes from : https://github.com/gboeing/urban-data-science/blob/master/19-Spatial-Analysis-and-Cartography/rtree-spatial-indexing.ipynb
def getplace(nameofplace):
	# get the boundary of some city
	gdf = ox.gdf_from_place(nameofplace)


	# make the geometry a multipolygon if it's not already
	geometry = gdf['geometry'].iloc[0]
	if isinstance(geometry, Polygon):
		geometry = MultiPolygon([geometry])

	return gdf, geometry

def parse(filepath):

	configfile = open(filepath, "r")
	flag = True
	dat = []
	while flag:
		l = configfile.readline()
		if (len(l)==0):
			flag = False
		elif(not(l[0]=="#" or l[0]=="\n")):
			dat.append([int(i) for i in l.split()])


	#num sides
	border=[]
	for i in range(0,dat[0[0]]):
		border.append(dat[i])

	return border

