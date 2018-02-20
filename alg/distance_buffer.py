import sys

from shapely.geometry import Polygon
from shapely.geometry import LinearRing
from shapely.geometry import LineString
from shapely.geometry import Point

from math import sqrt
from math import floor
from matplotlib import pyplot
from matplotlib.patches import Circle
from descartes.patch import PolygonPatch

#from shapely.figures import BLUE, SIZE, set_limits, plot_coords, color_isvalid
  



def plot_coords(ax, ob, color='#999999', zorder=1, alpha=1):
    x, y = ob.xy
    ax.plot(x, y, 'o', color=color, zorder=zorder, alpha=alpha)


def gen_dist_buffers(rad, area):
	mult = -1.5
	temp = area
	buffers = []
	#i = 1
	# est max num of buffers
	# idea: targetArea area / area of (mult*rad) , sqrt result , divide by 2
	# idk if this always works. this is for squares(?) ... vvv
	#nb = int(floor(sqrt((targetArea.area/(mult*rad*mult*rad)))/2))
	#print nb 
	for i in range(1,13):
		buffers.append(temp)
		temp = area.buffer(i*(mult*rad))
		#i=i+1

	return buffers

def plot_dist_buffers(ax, buffers):
	cnt = len(buffers)
	print "num of buffers "+str(cnt-1)
	print type(buffers)
	print type(buffers[0])

	for i in range(1, cnt):
		buf = buffers[i]
		if buf.is_valid:
			temppatch = PolygonPatch(buf, facecolor='#000000', alpha=(1.0/cnt))
			ax.add_patch(temppatch)

def eucdist(x1,y1,x2,y2):
	xdist = x2-x1
	ydist = y2-y1
	
	dist = sqrt(pow(xdist,2)+pow(ydist,2))
	return [dist, xdist, ydist]
	

def find_cand_points(ring, rad):
	candpts = []
	coords = list(ring.coords)
	numcoords=len(coords)
	cnt = 0
	x1 = coords[cnt][0]
	y1 = coords[cnt][1]
	
	while cnt < numcoords:

		x2 = coords[cnt+1][0]
		y2 = coords[cnt+1][1]

		# 1 - find length of current segment
		dist = eucdist(x1,y1,x2,y2)
		# if its on the 1st segment we use this dist, on current segment
		if dist[0] > sqrt(3)*rad:
			newx = x1 + eucdist[1]
			newy = y1 + eucdist[2]
			#cnt +=1
			candpts.append((newx,newy))
			x1 = newx
			y1 = newy

		# 2 - find segment that has that distance away on it
		while dist[0] <= sqrt(3)*rad:
			cnt += 1
			x1 = x2
			y1 = y2

			x2 = coords[cnt+1][0]
			y2 = coords[cnt+1][1]

			dist[0] += eucdist(x1,y1,x2,y2)[0]

			if dist[0] > sqrt(3)*rad:
				# this is the segment the dist is on -- 
				# next take the eucdist from point to cand point 
				print "poop"

	return pts




## distance fig

fig = pyplot.figure(1, dpi=90)

ax = fig.add_subplot(111)

targetArea = Polygon([(0,0),(0,400),(400,400),(400,0)],[[(101,74),(87,93),(99,119),(119,95)],[(200,70),(200,50),(300,20),(330,45),(370,45)]])
print targetArea.area

plot_coords(ax, targetArea.interiors[0])
plot_coords(ax, targetArea.exterior)
patch = PolygonPatch(targetArea, facecolor='#6699cc', edgecolor='#6699cc', alpha=0.1, zorder=2)
ax.add_patch(patch)

#pyplot.show()

## buffer fig
buffers = gen_dist_buffers(10.0, targetArea)
plot_dist_buffers(ax, buffers)

pts_buffer1 = find_cand_points(buffers[1],10.0)
print pts_buffer1

# experiment with buffer polys
#for i in range(0, len(buffers)):
#	print "buffer "+str(i)
#	if buffers[i].instanceof(shapely.geometry.polygon.MultiPolygon):
#		for j in range(0, len(buffers[i])):
#			print "linestring length "+str(j)+" "+buffers[i][j].length
#	else:
#		print "exterior length "+str(buffers[i].exterior.length)
#		print "interior length "+" "+str(buffers[i].interiors[0].length)


#show plot
pyplot.show()







#### FINISH vvv

def generateTargetArea(targetborder):
	numsides = targetborder[0[0]]
	outerborder = toTuples(targetborder[0])
	#holes = [for i in range(0, targetBorder[1[0]]): toTuples(targetborder[])


def toTuples(tuplelist):
	numtuples = tuplelist[0]
	tuples = []
	for i in range(1,numtuples):
		if not(i%2==0):
			tuples.append((tuplelist[i],tuplelist[i+1]))

	return tuples