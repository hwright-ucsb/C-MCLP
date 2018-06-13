import sys
import math
import shapely

from shapely.geometry import Polygon
from shapely.geometry import MultiPolygon
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
  
def unit_vector(vector):
	return (vector / np.linalg.norm(vector))

#plots points with scatter function
def plot_coords_list(ax, pts, dotsize=0.25, color="#999999", zorder=1, alpha=1):
	temp = zip(*pts)
	ax.scatter(temp[0], temp[1], dotsize, color=color, zorder=zorder, alpha=alpha)

#plots  pts....
def plot_coords(ax, ob, dotsize=0.01, color='#999999', zorder=1, alpha=1):
    x, y = ob.xy
    ax.plot(x, y, 'o', ms=dotsize, color=color, zorder=zorder, alpha=alpha)

#plots lines
def plot_line(ax, ob, ls='-', color='#999999', zorder=1, alpha=1):
	x, y = ob.xy
	ax.plot(x, y, linestyle=ls, color=color, zorder=zorder, alpha=alpha)


def plot_radii(ax, pts, rad, color="#999999", zorder=1, alpha=1):
	circles = []
	for pt in pts:
		c = Circle(pt, rad, color=color, zorder=zorder, alpha=alpha)
		circles.append(c)
		ax.add_patch(c)

def plot_fig(area):
    fig = pyplot.figure(1, figsize=(20,20), dpi=90)
    ax = fig.add_subplot(222)
    ax.axis('equal')
    interiors = area.interiors
    if len(interiors) > 0:
    	plot_coords(ax, area.interiors[0])
    plot_coords(ax, area.exterior)
    patch = PolygonPatch(area, facecolor='yellow', edgecolor='#6699cc', alpha=0.1, zorder=2)
    ax.add_patch(patch)
    return ax, fig


def gen_dist_buffers(rad, area, capstyle):
	# TODO : determine if this is the correct multiplier
	mult = -1.5*rad
	temp = area
	buffers = []
	i = 1

	while not temp.area == 0:
		if not temp.is_empty:
			buffers.append(temp)
		temp = area.buffer(i*(mult), capstyle)
		print temp.area
		i=i+1

	return buffers

def plot_dist_buffers(ax, buffers, alpha=1.00):
	cnt = len(buffers)
	print "num of buffers "+str(cnt-1)
	#print type(buffers)
	#print type(buffers[0].exterior)
	#print type(buffers[0].interiors)
	#print type(buffers[2])

	for i in range(1, cnt):
		buf = buffers[i]
		if buf.is_valid:
			temppatch = PolygonPatch(buf, facecolor='#000000', alpha=(alpha/cnt))
			ax.add_patch(temppatch)

def gen_initial_candpt(r, pt, deg):
	# multiplier - depends on the multiplier that is used to generate the dist. buffers
	mult = 1.5*r
	a = math.radians(deg) # deg = 60 or 120 
	newx = mult*math.cos(a)+pt[0]
	newy = mult*math.sin(a)+pt[1]
	return [newx,newy]

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
	print "ENTERING PT SEG CIRCLE ENTERING PT SEG CIRCLE ENTERING PT SEG CIRCLE ENTERING PT SEG CIRCLE "
	l1_len = eucdist(x1,y1,x2,y2)
	l2_len_all = eucdist(x2,y2,x3,y3)

	# "center" of the virtual circle
	Q = np.array([x1,y1])
	r = optdist
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
	    if (t1>=0 and t1<1):
	        potpt = P1 + t1*V
	        l2_len = eucdist(x2,y2,potpt[0], potpt[1])
	        if l2_len[0] <= l2_len_all[0]:
	        	pt=potpt
	        else:
	        	f=1
	        	print "pt extends past the end of the segment"
	    elif (t2>=0 and t2<1):
	        potpt = P1 + t2*V
	        l2_len = eucdist(x2,y2,potpt[0],potpt[1])
	        if l2_len[0] <= l2_len_all[0]:
	        	pt=potpt
	        else:
	        	f=1
	        	print "pt extends past the end of the segment"
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


	if lineeq[2] == 1: #horizontal
		print "horizontal"+ str(x1)+ " "+ str(optdist)
		if xydist[1]<0 or xydist[2]<0:
			optdist = -optdist
		newcoords[0] = x1 + optdist
		if (optdist>0 and newcoords[0] > x2) or (optdist<0 and newcoords[0]<x2):
			print " SET FLAG around corner x"
			flag = 1
			# TODO: this is the case that it bends around the corner,
			# so find distances from center of potl sites 
		newcoords[1] = y1
	elif lineeq[2] == 2: # vertical
		if xydist[1]<0 or xydist[2]<0:
			optdist = -optdist
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
		checkdist = eucdist(x1,y1,newcoords[0],newcoords[1])

		if checkdist[0] > xydist[0]:
			print "DISTS "+ str(checkdist[0]) + "      "+str(xydist[0])
			flag = 1
			print "the new coord was past the end of the line"


	return newcoords, flag

# always to - from
def ang_bw_vectors(to1, from1, to2, from2):
	v1_u = unit_vector(np.array(to1)-np.array(from1))
	v2_u = unit_vector(np.array(to2)-np.array(from2))
	return (np.arccos(np.clip(np.dot(v1_u,v2_u),-1.0, 1.0))) # we do the clipping extra step to handle when the vecs are in the same or opposite directions


#find the inflection point of the coord list 
def find_infl_pt(newcandpt, lastcandpt, coords, ind, r):
	ptnotfound = True
	pt = (np.NAN, np.NAN)
	n = len(coords)-1
	#cnt = ind+2
	from1 = lastcandpt
	to1 = newcandpt
	#from2 = lastcandpt
	#to2 = eol

	eol = coords[ind]
	cntt = 0
	cnt = ind
	flag = False
	maxdist = 0.0
	d1 = [sqrt(3)*r + 1]
	while d1[0] > (sqrt(3)*r)/2:
		# try additive dist first
		d1 = eucdist(lastcandpt[0], lastcandpt[1], eol[0], eol[1])
		d2 = eucdist(eol[0], eol[0], newcandpt[0], newcandpt[1])
		print d1
		print d2
		d = d1[0]+d2[0]
		if d > maxdist:
			pt = eol
			maxdist = d
			flag = True
			cntt=cntt+1
		
		cnt = cnt+1
		eol=coords[cnt]
		print "MAXDIST: "+str(maxdist)+"   "+str(cntt)







 #FIX THIS 
	# while ptnotfound and cnt <=n:
	# 	deg = np.rad2deg(ang_bw_vectors(from1, to1, from1, to2))
	# 	if deg > 30.0: #and deg < 150.0:
	# 		pt = eol
	# 		ptnotfound = False
	# 		flag = True
	# 	#pt1 = pt2
	# 	to1 = pt3
	# 	to2 = coords[cnt]
	# 	if cnt+1 > n:
	# 		ptnotfound = False
	# 	cnt = cnt+1
	
	return pt, flag




# Take in list of polygons (which have 1 exterior and any # of interiors)
def find_cand_points(list_of_poly, rad):
	candpts = []
	for i in range(0,len(list_of_poly)):
		tpoly= list_of_poly[i]
		candpts.append(find_cand_points_list(list(tpoly.exterior.coords),rad))
		for j in range(0, len(tpoly.interiors)):
			candpts.append(find_cand_points_list(list(tpoly.interiors[j].coords), rad))

	return candpts


# REFACTORING of find_cand_points_list - try using circles ?
def find_cand_points_list_ref(coords, rad):
	inflpts = []
	candpts = []
	places = []
	optdist = sqrt(3)*rad
	numcoords = len(coords)

	x1 = coords[0][0]
	y1 = coords[0][1]
	lastcandpt = [x1,y1]
	lastsegend = [x1,y1]
	cnt = 1

	while cnt < numcoords:
		# if we're at the end of the coord list, close the loop 
		if cnt>numcoords-1:
			# TODO: handle case where we need to put an extra cand pt where we're closing the loop
			print "closing ring/polygon at cnt = "+str(cnt)
			x2 = coords[0][0]	#loop back around to the first coord 
			y2 = coords[0][1]
			cnt = numcoords
			break # TODO:  is this necessary??
		else: # not at the end of the loop, keep looking for cand pts.
			x2 = coords[cnt][0]
			y2 = coords[cnt][1]

			# choose the first coord as the first cand pt for now, later
			# this will be more sophisticated (perhaps randomized)
			print "LINE EQ: "
			line = line_eq(x1,y1,x2,y2)
			f = 0 # what is this flag for again?

			# cant remember if it is undefd. for horizontal or vert lines - check
			while f==0:
				if not line[2]==6: #horizontal 2=vertical
					candpt, f = pt_seg_circle(optdist, lastcandpt[0],lastcandpt[1],x1,y1,x2,y2)
					if f==0: #cand pt found
						# NEW 3/12 - check the angle b/w line formed w/ candpts; if above a threshold, there is a sharp curve ahead
						lastsegcnt = cnt
						b = len(candpts)-1
						try:
							lastcandpt = candpts[b]
							ang = np.rad2deg(ang_bw_vectors(np.array(lastcandpt), np.array(candpt), np.array(lastcandpt), np.array(coords[lastsegcnt])))
						
							if ang > 30.0 :#and ang < 150.0: # angle between two candpt sites and the end of the current line
								print "THIS IS A PLACE"
								places.append(lastcandpt)
								inflpt, iflag = find_infl_pt(candpt, lastcandpt, coords, lastsegcnt, rad)
								if iflag == True:
									inflpts.append(inflpt)
								print str(ang)
								print str(candpt)+" "+str(lastcandpt)+" "+str(x1)+", "+str(y1)
						except IndexError:
							print "INDEXERROR: ANGLE NOT CALCULATED"
						candpts.append(candpt)
						lastcandpt = [candpt[0],candpt[1]]

						x1 = candpt[0]
						y1 = candpt[1]
						# dont change p1 and p2 because they are not changed yet; were still on this line
					else: #cand pt not found
						# no interception with this line, get the next line
						cnt+=1
						if cnt>numcoords-1:
							break
						x1 = x2
						y1 = y2
						x2 = coords[cnt][0]
						y2 = coords[cnt][1]
					
					#n = len(candpts)-1
					#x1 = candpts[n][0]
					#y1 = candpts[n][1]


				else:
					# do the other procedure
					break # for now

	return candpts, inflpts, places


def find_cand_points_polys(poly, r):
    candpts = []
    inflpts = []
    places = []
    candpt, inflpt, place = find_cand_points_list_ref(list(poly.exterior.coords),r)
    candpts.append(candpt)
    inflpts.append(inflpt)
    places.append(place)
    #candpts.append(find_cand_points_list_ref(list(poly.exterior.coords),r))
    for i in range(0, len(poly.interiors)):
        candpt, inflpt, place = find_cand_points_list_ref(list(poly.interiors[i].coords),r)
        inflpts.append(inflpt)
        candpts.append(candpt)
        places.append(place)
        #candpts.append(find_cand_points_list_ref(list(poly.interiors[i].coords),r))
    
    return candpts, inflpts, places

def find_cand_points_buffers(buffers, r):
    candpts = []
    inflpts = []
    places = []
    for i in range(0, len(buffers)):
        if isinstance(buffers[i], Polygon):
            print "poly"
            candpt, inflpt, place = find_cand_points_polys(buffers[i], r)
            candpts.append(candpt)
            inflpt.append(inflpt)
            places.append(place)
            #candpts.append(find_cand_points_polys(buffers[i], r))
        elif isinstance(buffers[i], MultiPolygon):
            print "multi"
            for j in range(0, len(buffers[i])):
            	candpt, inflpt, place = find_cand_points_polys(buffers[i].geoms[j], r)
            	candpts.append(candpt)
            	inflpts.append(inflpt)
            	places.append(place)
                #candpts.append(find_cand_points_polys(buffers[i].geoms[j], r))
    return candpts, inflpts, places





# Takes in a coordinate sequence and finds the cand pts
#	Note: a polygon +/or buffers is made of rings (exterior & interior)
def find_cand_points_list(coords, rad):
	candpts = []
	optdist = sqrt(3)*rad
	numcoords=len(coords)

	x1 = coords[0][0]
	y1 = coords[0][1]
	cnt = 1

	while cnt < numcoords:
		print " FIRST CNT "+str(cnt)
		if cnt>numcoords-1:
			# TODO : make sure this part works 
			print "closing ring/polygon at cnt = "+str(cnt)
			x2 = coords[0][0]	#loop back around to the first coord 
			y2 = coords[0][1]
			cnt = numcoords
			break # TODO:  is this necessary??
		else:
			x2 = coords[cnt][0]
			y2 = coords[cnt][1]
			cnt+=1

		# 1 - find length of current segment
		dist = eucdist(x1,y1,x2,y2)
		print str(dist)

		# if its on the 1st segment we use this dist, on current segment
		if dist[0] > optdist:
			#l = polyfit([x1,y1],[x2,y2],1)
			#print str(l)
			# TODO: fix this naive method of guessing how many cand pts there are on the seg.
			num_cand_sites = int(floor(dist[0]/optdist)) # num of cand sites on this segment
			print num_cand_sites
			line = line_eq(x1,y1,x2,y2)
			print line
			f = 0 
			while f == 0:
				if cnt>numcoords-1: # CHECK
					break
				candpt, f  = pt_on_line(line, dist, optdist,x1,y1,x2,y2)
				if f == 0:
					x1 = candpt[0]
					y1 = candpt[1]
					#cnt +=1
					candpts.append((x1,y1))
					print str(candpt)
					print " SEC CNT "+str(cnt)
					cnt+=1
					#x1 = newx
					#y1 = newy
				else:
					# 1 - calculate point around bend with procedure
					# 2 - set new cand pt as start of next line and cont. w/ normal procedure
					#newpt = pt_on_bend()
					#TODO: handle corner calculations
					print "FLAG bends around corner"
					x1=x2
					y1=y2

					if cnt >=numcoords:
						break
					x2 = coords[cnt][0]
					y2 = coords[cnt][1]

					print " THIRD CNT "+str(cnt)
					print "x1 "+str(x1)+" y1 "+str(y1)+" x2 "+str(x2)+" y2 "+str(y2)
					
					cnt+=1
					lastcandpt = candpts[len(candpts)-1]
					potlpt, f = pt_seg_circle(optdist, lastcandpt[0],lastcandpt[1], x1,y1,x2,y2)
					if f == 0:
						candpts.append((potlpt[0],potlpt[1]))
					break

		# 2 - find segment that has that distance away on it
		# TODO : need to take the last candidate pt and get the dist from that
		# while f ==1: # TODO : investigate if it should be 1.5*r everywhere!!!
		# 	x1=x2
		# 	y1=y2
		# 	cnt=+1
		# 	x2 = coords[cnt][0]
		# 	y2 = coords[cnt][1]

		# 	lastcandpt = candpts[len(candpts)-1]
		# 	potlpt, f = pt_seg_circle(optdist, lastcandpt[0],lastcandpt[1], x1,y1,x2,y2)
		# 	if f == 0:
		# 		candpts.append((potlpt[0],potlpt[1]))

	

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