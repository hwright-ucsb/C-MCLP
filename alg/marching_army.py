from shapely.geometry import LineString, Point, Polygon
from shapely import geometry
from numpy.linalg import norm
from numpy import array
from numpy.random import random
from math import pi, sqrt

import distance_buffer as db


def gen_rand_pt(lo, hi):
    return ((lo - hi) * random() + hi)

def gen_starting_line(bounds):
	# randomly pick 2 pts in area then get unit vector of that line
    pt1 = (gen_rand_pt(bounds[0],bounds[2]),gen_rand_pt(bounds[1],bounds[3]))
    pt2 = (gen_rand_pt(bounds[0],bounds[2]),gen_rand_pt(bounds[1],bounds[3]))
    # this switch doesnt really make a difference
    temp = pt2
    pt2 = array(pt1)
    pt1 = array(temp)
    v = pt2 - pt1
    line = v/norm(v)
    # TODO : vvv change this to a while loop
    #line = []
    #if not ((pt1[0]==pt2[0]) and (pt1[1]==pt2[1])):
    #    v = pt2 - pt1
    #    line = v / np.linalg.norm(v)
    #else:
    #    line = gen_starting_line(bounds)
    return line, pt1,pt2


# line[0] is starting pt, line[1] endpt; poly is a ring
def intersect_line_region(pt1,linevec, poly):
    
    intpt = []
    i = 1
    while len(intpt) < 2:
    #    # the line is a unit vector, get it long enough to make sure it intersects before we return the pt
        temp1 = pt1+ (100000.0*i)*linevec
        temp2 = pt1 - (100000.0*i)*linevec
        newline = LineString([temp1,temp2])
        intpt = newline.intersection(poly)
        i+=1
    
    return LineString(intpt)


# Note : this method will just pick the first rand pt found in the region
def gen_rand_pt_on_line(pt1, pt2, vec, region):
    pt = []
    mag = db.eucdist(pt1[0], pt1[1], pt2[0], pt2[1])
    randdist = gen_rand_pt(0, mag[0])
    randptp = pt1 + randdist*vec
    randptn = pt1 - randdist*vec
    randptp2 = pt2 + randdist*vec
    randptn2 = pt2 - randdist*vec

    if region.contains(Point(randptp)):
        pt = randptp
    elif region.contains(Point(randptn)):
        pt = randptn
    elif region.contains(Point(randptp2)):
        pt = randptp2
    elif region.contains(Point(randptn2)):
        pt = randptn2
    else:
        print "contains none...."
    #randptn = pt1 - randdist*vec
    #if 
    return pt

def gen_next_line_startpt(startpt, linevec, r):
    perpvec = [-linevec[1], linevec[0]]
    perpdist = (1.5)*r
    pardist = (sqrt(3)/2.0)*r
    temppt = startpt + perpdist*array(perpvec)
    newstartpt = temppt + pardist*linevec
    return newstartpt

def gen_cand_pts_on_line(poly, linevec, startpt, r):
    #main algorithm, start from the randomly generated intersecting line & then 
    # keeps track of running missed pts; stop once we get to a thresh
    
    eps = 0.01 # epsilon for tolerance of intersection of region and coverage area(r)
    covarea = pi*r*r
    nummissedpts = 0
    pt = startpt
    candpts = []
    tempregion = poly
    regiondiff = poly.difference(Point(startpt).buffer(r))
    bordercands = []

    startpts = []    
    lastpt = startpt 

    optdist = sqrt(3)*r
    startpt = startpt + optdist*linevec
    #if regiondiff.area < tempregion.area:
    #    candpts.append(startpt)
    #    startpts.append(startpt)
   
            
    tempregion = regiondiff
    while nummissedpts < 50:
            temppt = pt + optdist*linevec
            pt = temppt
            regiondiff = tempregion.difference(Point(pt).buffer(r))
            rdarea = regiondiff.area
            trarea = tempregion.area
            if rdarea < trarea:
                candpts.append(temppt)
                tempregion = regiondiff
                nummissedpts = 0
                lastpt = temppt
                #if (trarea - rdarea) < covarea:
                    # this means the new candpt's covarea is on a border/hole
                #    bordercands.append(temppt)
            else:
                nummissedpts+=1
                
    
    
    #print "MOVE TO NEG DIR PT LOOP "+ str(len(candpts))
    
    
    nummissedpts = 0
    pt = startpt
    tempregion = regiondiff
    while nummissedpts < 50:
            temppt = pt - optdist*linevec
            pt = temppt
            regiondiff = tempregion.difference(Point(pt).buffer(r))
            rdarea = regiondiff.area
            trarea = tempregion.area
            if rdarea < trarea:
                candpts.append(temppt)
                tempregion = regiondiff
                nummissedpts = 0
                lastpt = temppt
                #if (trarea - rdarea) < covarea:
                    # this means the new candpt's covarea is on a border/hole
                #    bordercands.append(temppt)
            else:
                nummissedpts+=1

    return candpts, bordercands, startpts#, tempregion

def gen_cand_pts(poly, startvec, startpt, r):
    nummissedlines = 0
    candpts = []
    borderpts = []
    starts = []
    nextstpt = startpt
    #tempregions = []
    tempregion = poly
    #positive direction loop
    while nummissedlines < 100:
        temppts, tempbpts, stpts = gen_cand_pts_on_line(tempregion, startvec, nextstpt, r)
        if len(temppts) > 1:
            candpts.append(temppts)
            borderpts.append(tempbpts)
            starts.append(stpts)
            #tempregions.append(tempregion)
            nextstpt = gen_next_line_startpt(nextstpt, startvec, r)
            nummissedlines=0
        else:
            nummissedlines+=1
            
    nummissedlines=0
    
    #print "MOVE TO NEG DIR LINE LOOP "+ str(len(tempregions))
    nextstpt = gen_next_line_startpt(startpt, startvec, -1*r)
    while nummissedlines < 100:
        temppts, tempbpts, stpts = gen_cand_pts_on_line(tempregion, startvec, nextstpt, r)
        if len(temppts) > 1:
            candpts.append(temppts)
            borderpts.append(tempbpts)
            starts.append(stpts)
            #tempregions.append(tempregion)
            nextstpt = gen_next_line_startpt(nextstpt, startvec, -1*r)
            nummissedlines=0
        else:
            nummissedlines+=1     
    return candpts, borderpts, stpts#, tempregions