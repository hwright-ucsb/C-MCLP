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

def plot_radii(ax, pts, rad, color="#999999", zorder=1, alpha=1):
	circles = []
	for pt in pts:
		c = Circle(pt, rad, color=color, zorder=zorder, alpha=alpha)
		circles.append(c)
		ax.add_patch(c)


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



# Given the eq. of line 1 (l1) and the endpt of line 2 (l2) that are joined by a pt,
#	find the pt on l2 s.t. the straight line distance from that point to the original
#	pt on l1 is the given distance
def pt_on_bend(l1, optdist, x1, y1, x2, y2, x3, y3):
	r = optdist
	l1_len = eucdist(x1,y1,x2,y2)
	l2_len_all = eucdist(x2,y2,x3,y3)
	l2 = line_eq(x2,y2,x3,y3)
	# 1 - find angle b/w segments (alpha) w/ slopes m1 and m2
	alpha = np.arctan((l1[0]-l2[0])/(1+l1[0]*l2[0]))
	# 2 - solve triangle for missing angle mu
	sinalpha = np.sin(alpha)
	mu = np.arcsin((2.0*l1_len[0]*sinalpha)/(3.0*r))
	# 3 - Subtract mu and alpha to get nu, the missing angle of the triangle
	nu = np.pi - (mu+alpha)
	# 4 - use law of sines to solve for missing triangle length l2
	l2_len = (3.0*nu*r)/(2.0*sinalpha)

# Given the line segment and the circle, determine their intersection, if any
#	src: stackexchange .. line-segment-to-circle-collision-algorithm
def pt_seg_circle(optdist, x1, y1, x2, y2, x3, y3):
	l1_len = eucdist(x1,y1,x2,y2)
	l2_len = eucdist(x2,y2,x3,y3)

	# "center" of the virtual circle
	Q = np.array([x1,y1])
	r = 10.0
	# start of line segment
	P1 = np.array([x2,y2])
	# vector along line segment
	V = np.array([x3,y3]) - P1
	a = V.dot(V)
	b = 2.0 * V.dot(P1 - Q)
	c = P1.dot(P1) + Q.dot(Q) - 2.0 * P1.dot(Q) - r**2
	disc = b**2 - 4*a*c
	pt = np.array([np.NaN, np.NaN])
	f = 0
	if disc < 0:
		f = 1
		print "no intercept"
	else:
	    sqrt_disc = sqrt(disc)
	    t1 = (-b + sqrt_disc)/(2*a)
	    t2 = (-b - sqrt_disc)/(2*a)
	    print t1, t2
	    if (t1>0 and t1<1):
	        pt = P1 + t1*V
	    elif (t2>0 and t2<1):
	        pt = P1 + t2*V
	    else:
	    	f = 1
	        print "no intercept (there would be one if the line segment cont'd)"  
	            
	return pt,f



# Given the equation of a line segment with form y=mx+b between (x1,y1)
#	and (x2,y2), calculate point on line that is dist away from (x1,y1)
def pt_on_line(lineeq, xydist, optdist, x1, y1, x2, y2):
	flag = 0
	newcoords = [0.,0.]
	print str(x1)+" "+str(y1)+" "+str(x2)+" "+str(y2)
	print (xydist)
	if xydist[1]<0 or xydist[2]<0:
		optdist = -optdist


	if lineeq[2] == 1: #horizontal
		print "horizontal"+ str(x1)+ " "+ str(optdist)
		newcoords[0] = x1 + optdist
		if (optdist>0 and newcoords[0] > x2) or (optdist<0 and newcoords[0]<x2):
			print " SET FLAG around corner x"
			flag = 1
			# TODO: this is the case that it bends around the corner,
			# so find distances from center of potl sites 
		newcoords[1] = y1
	elif lineeq[2] == 2: # vertical
		print "vertical "+ str(y1)+ " "+ str(optdist)
		newcoords[1] = y1 + optdist
		if (optdist>0 and newcoords[1] > y2) or (optdist<0 and newcoords[1] < y2):
		 	print " SET FLAG around corner y"
		 	flag = 1
		 	# TODO: this is the case that it bends around the corner,
			# so find distances from center of potl sites 
		newcoords[0] = x1
	else: # sloped line
		v = np.array([x2,y2]) - np.array([x1,y1])
		u = v / np.linalg.norm(v)
		newcoords = np.array([x1,y1]) + optdist*u
		#if newcoords[0]

	return newcoords, flag

# Take in list of polygons (which have 1 exterior and any # of interiors)
def find_cand_points(list_of_poly, rad):
	candpts = []
	for i in range(0,len(list_of_poly)):
		tpoly= list_of_poly[i]
		candpts.append(find_cand_points_list(list(tpoly.exterior.coords),rad))
		for j in range(0, len(tpoly.interiors)):
			candpts.append(find_cand_points_list(list(tpoly.interiors[j].coords), rad))

	return candpts


# Takes in a coordinate sequence and finds the cand pts
#	Note: a polygon +/or buffers is made of rings (exterior & interior)
def find_cand_points_list(coords, rad):
	candpts = []
	optdist = sqrt(3)*rad
	# TODO: fix this naive method of guessing how many cand pts there are on the seg.
	numcoords=len(coords)
	cnt = 0
	x1 = coords[cnt][0]
	y1 = coords[cnt][1]
	
	while cnt < numcoords:
		if cnt+1>=numcoords-1:
			# TODO : make sure this part works 
			print "closing ring/polygon at cnt = "+str(cnt)
			x2 = coords[0][0]	#loop back around to the first coord 
			y2 = coords[0][1]
			cnt = numcoords
		else:
			cnt+=1
			x2 = coords[cnt][0]
			y2 = coords[cnt][1]
			

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
			print line
			f = 0 
			while f == 0:
				candpt, f  = pt_on_line(line, dist, optdist,x1,y1,x2,y2)
				if f == 0:
					x1 = candpt[0]
					y1 = candpt[1]
					#cnt +=1
					candpts.append((x1,y1))
					print str(candpt)
					#x1 = newx
					#y1 = newy
				else:
					# 1 - calculate point around bend with procedure
					# 2 - set new cand pt as start of next line and cont. w/ normal procedure
					#newpt = pt_on_bend()
					#TODO: handle corner calculations
					print "FLAG bends around corner"
					#x1 = x2
					#y1 = y2
					break

		# 2 - find segment that has that distance away on it
		# TODO : need to take the last candidate pt and get the dist from that
		if f ==1: # TODO : investigate if it should be 1.5*r everywhere!!!
			x1=x2
			y1=y2
			cnt=+1
			x2 = coords[cnt][0]
			y2 = coords[cnt][1]

			lastcandpt = candpts[len(candpts)-1]
			potlpt, f = pt_seg_circle(optdist, lastcandpt[0],lastcandpt[1], x1,y1,x2,y2)
			if f == 0:
				candpts.append((potlpt[0],potlpt[1]))

	

	return candpts






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