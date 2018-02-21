import sys

from shapely.geometry import Polygon
from shapely.geometry import LinearRing
from shapely.geometry import LineString
from shapely.geometry import Point

import numpy as np

from scipy.spatial.distance import euclidean

from math import sqrt
from math import floor
from matplotlib import pyplot
from matplotlib.patches import Circle
from descartes.patch import PolygonPatch

#from shapely.figures import BLUE, SIZE, set_limits, plot_coords, color_isvalid
  

def plot_coords_list(ax, pts, dotsize=0.25, color="#999999", zorder=1, alpha=1):
	temp = zip(*pts)
	ax.scatter(temp[0], temp[1], dotsize, color=color, zorder=zorder, alpha=alpha)

def plot_coords(ax, ob, color='#999999', zorder=1, alpha=1):
    x, y = ob.xy
    ax.plot(x, y, 'o', color=color, zorder=zorder, alpha=alpha)


def gen_dist_buffers(rad, area):
	mult = -1.5
	temp = area
	buffers = []
	#i = 1
	# est max num of buffers
	# 
	#print nb 
	for i in range(1,13):
		buffers.append(temp)
		temp = area.buffer(i*(mult*rad))
		#i=i+1

	return buffers

def plot_dist_buffers(ax, buffers):
	cnt = len(buffers)
	print "num of buffers "+str(cnt-1)
	#print type(buffers)
	#print type(buffers[0].exterior)
	#print type(buffers[0].interiors)
	#print type(buffers[2])

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
	
def line_eq(x1, y1, x2, y2):
	m  = 0
	b = 0
	t = 0
	try:
		num = (float(y1)-y2)
		if num == 0: # horizontal line
			print "horizontal line y = "+str(y1)
			t = 1 # horizontal == 1
		else:
			m = num/(float(x1)-x2)
			b = y1 - m*x1
			print ("y = "+str(m)+"*x + "+str(b))
	except ZeroDivisionError: #vertical line
		print "vertical line x = "+str(x1)
		t = 2 # vertical == 2

	return [m,b,t]


# Given the equation of a line segment with form y=mx+b between (x1,y1)
#	and (x2,y2), calculate point on line that is dist away from (x1,y1)
def pt_on_line(lineeq, dist, x1, y1, x2, y2):
	flag = 0
	newcoords = [0.,0.]
	if lineeq[2] == 1: #horizontal
		print "horizontal"+ str(x1)+ " "+ str(dist)
		newcoords[0] = x1 + dist
		if newcoords[0] > x2:
			print " SET FLAG around corner x"
			flag = 1
			# TODO: this is the case that it bends around the corner,
			# so find distances from center of potl sites 
		newcoords[1] = y1
	elif lineeq[2] == 2: # vertical
		print "vertical "+ str(y1)+ " "+ str(dist)
		newcoords[1] = y1 + dist
		if newcoords[1] > y2:
		 	print " SET FLAG around corner y"
		 	flag = 1
		 	# TODO: this is the case that it bends around the corner,
			# so find distances from center of potl sites 
		newcoords[0] = x1
	else: # sloped line
		v = np.array([x2,y2]) - np.array([x1,y1])
		u = v / np.linalg.norm(v)
		newcoords = np.array([x1,y1]) + dist*u

	return newcoords, flag

def find_cand_points(poly, rad):
	candpts = []
	optdist = sqrt(3)*rad
	# TODO : have to deal with interior and exterior rings
	coords = list(poly.exterior.coords)
	numcoords=len(coords)
	cnt = 0
	x1 = coords[cnt][0]
	y1 = coords[cnt][1]
	
	while cnt < numcoords:
		if cnt+1>=numcoords-1:
			print "closing ring/polygon at cnt = "+str(cnt)
			x2 = coords[0][0]	#loop back around to the first coord 
			y2 = coords[0][1]
			cnt = numcoords
		else:
			x2 = coords[cnt+1][0]
			y2 = coords[cnt+1][1]
			cnt+=1

		# 1 - find length of current segment
		dist = eucdist(x1,y1,x2,y2)
		print str(dist)

		# if its on the 1st segment we use this dist, on current segment
		if dist[0] > optdist:
			#l = polyfit([x1,y1],[x2,y2],1)
			#print str(l)
			num_cand_sites = int(floor(dist[0]/optdist)) # num of cand sites on this segment
			print num_cand_sites
			line = line_eq(x1,y1,x2,y2)
			#for i in range(0, num_cand_sites):
			f = 0 
			while f == 0:
				candpt, f  = pt_on_line(line, optdist,x1,y1,x2,y2)
				if f == 0:
					x1 = candpt[0]
					y1 = candpt[1]
					#cnt +=1
					candpts.append((x1,y1))
					print str(candpt)
					#x1 = newx
					#y1 = newy
				else:
					#TODO: handle corner calculations
					print "FLAG bends around corner"
					x1 = x2
					y1 = y2
					break

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

	return candpts




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

pts_buffer1 = find_cand_points(buffers[0],10.0)
print pts_buffer1
plot_coords_list(ax, pts_buffer1, color="#000000", zorder=3)

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