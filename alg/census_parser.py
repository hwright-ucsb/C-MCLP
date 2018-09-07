import shapely
from shapely import geometry
from shapely.ops import transform
from descartes import PolygonPatch
from matplotlib.patches import Arc
from shapely.geometry import Point, Polygon, MultiPolygon, LineString
import sys
sys.path.append('../../../alg')
import distance_buffer as db, border_generators as bg, marching_army as ma


import math
import pyproj
import fiona

import matplotlib.pyplot as plt, pandas as pd, geopandas as gpd, numpy as np
from matplotlib import pyplot
from functools import partial
from scipy.stats import truncnorm

#transform to meters so we can do math on the RoI
proj2 = partial(pyproj.transform, pyproj.Proj(init='epsg:4269'),pyproj.Proj(init='epsg:3857'))


# filter collection based on a places commonly known name
# filter_value : String - the value by which you wish to filter property prop by 
# prop: String - the property in the "properties" that you want to filter fcollection by
# fcollection : a shapely collection of data 
# RETURNS: list - of matching entries found in fcollection
def filter_by_property(filter_value, prop, fcollection):
    matches = []
    for i in range(0, len(fcollection)):
        if fcollection[i]['properties'][prop] == filter_value:
            matches.append(fcollection[i])
  
    return matches


def remove_duplicate_blocks(collec):
    temp_data = []
    dups = []
    listofids = []
    for i in range(0, len(collec)):
        cur = collec[i]
        tid = str(cur['properties']['BLOCKID10'])
        if tid not in listofids:
            listofids.append(tid)
            temp_data.append(cur)
        else:
            dups.append(cur)
            
    # MAKE SURE NONE OF THE BLOCKS YOU ELIMINATED HAVE A DIFFERENT POP
    for i in range(0, len(dups)):
        thedup = str(dups[i]['properties']['BLOCKID10'])
        for j in range(0, len(temp_data)):
            if thedup == str(temp_data[j]['properties']['BLOCKID10']):
                p1 = dups[i]['properties']['POP10']
                p2 = temp_data[j]['properties']['POP10']
                if not p1==p2:
                    print(p1)
                    print(p2)
                    
    return temp_data


# generate N random coords within the SQUARE BOUNDS bds [minx,miny,maxx,maxy]
# RETURNS: xs - x coordinate list of rand coords
#          ys - y coord list of rand coords 
def gen_n_random_coords(bds, n):
    xs = []
    ys = []
    for i in range(0,n):
        xs.append(ma.gen_rand_pt(bds[0], bds[2]))
        ys.append(ma.gen_rand_pt(bds[1],bds[3]))
    return xs, ys

#generate points for the whole RoI
def filter_random_pts_by_RoI(x_rands, y_rands, RoI_bg_geoms, pop_data_schema):

    # turn x list and y list into coords
    random_pts_all = gpd.GeoDataFrame(data={'x':x_rands, 'y':y_rands})
    random_pts_all['geometry'] = random_pts_all.apply(lambda row: Point((row['x'], row['y'])), axis=1)
    # create index for all random pts
    all_pts_sindex = random_pts_all.sindex

    #now find ones that intercept the region by BG

    pts_per_bg = []
    pops_per_bg= []


    for i in range(0,len(pop_data_schema)):
        tbg_pop = pop_data_schema[i]['properties']['POP10']
        pops_per_bg.append(tbg_pop)
        all_precise_matches = []
        if tbg_pop > 0: # only check for point intersections if the population >0
            poly = RoI_bg_geoms[i]
            polybds = poly.bounds 
        
            all_possible_matches_index = list(all_pts_sindex.intersection(polybds))
            all_possible_matches = random_pts_all.iloc[all_possible_matches_index]
            all_precise_matches = all_possible_matches[all_possible_matches.intersects(poly)]
        
            while len(all_precise_matches)==0:
            #gen new random points focused on the BG at hand
                tx, ty = gen_n_random_coords(polybds, 10)
                trandom_pts = gpd.GeoDataFrame(data={'x':tx, 'y':ty})
                trandom_pts['geometry'] = trandom_pts.apply(lambda row: Point((row['x'], row['y'])), axis=1)
                trandom_pts_sindex = trandom_pts.sindex
            
                tpossible_matches_index = list(trandom_pts_sindex.intersection(polybds))
                tpossible_matches = trandom_pts.iloc[tpossible_matches_index]
                all_precise_matches = tpossible_matches[tpossible_matches.intersects(poly)]
            
        
            # make sure there's not more points than people in the bg
            if len(all_precise_matches) > tbg_pop:
                all_precise_matches = all_precise_matches[0:tbg_pop]
        
        pts_per_bg.append(all_precise_matches)
        
    return pts_per_bg#, pops_per_bg

def assign_pops_to_pts(pts_per_bg, pops_per_bg):
    pop_partitions = []
    
    for i in range(0, len(pts_per_bg)):
        bgpop = pops_per_bg[i]['properties']['POP10']
        bgpts = pts_per_bg[i]
        numpts = len(bgpts)
        pop_remains = bgpop
        pt_remains = numpts
        #print "NEW IIIII"
        
        partition = []
        partition_not_found = True
        # only partition pts if nonzero pop
        if bgpop == 0:
            #print "0 pop"
            pop_partitions.append(partition)
            partition_not_found = False
        elif bgpop == numpts: # simple partition each pt gets 1
            #print "each part gets 1"
            partition = [1]*numpts
            pop_partitions.append(partition)
            partition_not_found = False
        elif numpts > bgpop: # too many pts for the total pop, shave off some pts and assign 1's
            #print "too many pops for the pt"
            pts_per_bg[i] = pts_per_bg[i][0:bgpop]
            partition = [1]*bgpop
            pop_partitions.append(partition)
            partition_not_found = False
            
        # non simple partition case; now find the partition
        while partition_not_found: # loop until we find the partition
            #print "IN WHILE"
            #print i
            #print "bgpop "+str(bgpop)
            #print "numpts "+str(numpts)
            #print "pop_remains "+str(pop_remains)
            #print "pt_remains "+str(pt_remains)
                
            if pop_remains == pt_remains:
                #print "pop remain is pt remain"
                # if the # pts left ever =='s the # pop left, 
                # just assign 1's to the rest of the pts
                partition.extend(([1]*pop_remains))
                pop_partitions.append(partition)
                partition_not_found = False
                continue
            elif pop_remains < pt_remains: # not enough pop left to cover the remaining pts; start over
                pop_remains = bgpop
                pt_remains = numpts
                partition = []
                
            if pt_remains == 1:
                #print "last one,"
                partition.append(pop_remains)
                pop_partitions.append(partition)
                partition_not_found = False
                continue
            
                    
            randnum = np.random.randint(1,pop_remains)
            partition.append(randnum)
            # decrease remaining pts and pops accordingly
            pop_remains = pop_remains-randnum
            pt_remains -= 1
    
    #print pop_partitions
    return pop_partitions


def test_pop_pt_agreement(pops_bg, popassigs, pts_bg):
    for i in range(0, len(popassigs)):
        t = popassigs[i]
        pts = pts_bg[i]
        #print t
        tsum = 0
        ptcnt = 0
        for j in range(0, len(t)):
            tsum +=t[j]
        for j in range(0, len(pts)):
            ptcnt += 1
        if not tsum == pops_bg[i]['properties']['POP10']:
            print("POP DISAGREEMENT")
            print("I "+str(i))
            print(tsum)
            print(pops_bg[i])
        if not ptcnt == len(t):
            print("PTCNT DISAGREEMENT")
            print(ptcnt)
            print(len(t))

# DATA PROCESSING CODE
#get the geom and pop data for blockgroups for a PLACE
def get_place_blockgroup_data(bg_data_str, places_data_str, faces_data_str, countycode, placename):
    state_bg_data = fiona.open(bg_data_str)
    state_places_data = fiona.open(places_data_str)
    county_faces_data = fiona.open(faces_data_str)

    place_bg_data = filter_by_property(placename, 'NAME10', state_places_data)

    # get the shapely geometry of the BGs for place
    place_geoms = geometry.shape(place_bg_data[0]['geometry'])
    t_place_geoms = transform(proj2,place_geoms)
    # turn it into a multipolygon for easier processing
    place_geoms_poly = MultiPolygon([place_geoms])

    placecode = str(place_bg_data[0]['properties']['PLACEFP10'])

    place_faces = filter_by_property(placecode, 'PLACEFP10', county_faces_data)
    county_bg_data = filter_by_property(countycode, 'COUNTYFP10', state_bg_data)

    # now we have all the BG info we need to get the populations in PLACE
    place_pop_data = []

    for i in range(0, len(place_faces)):
        cur = place_faces[i]
        matchstr = str(cur['properties']['STATEFP10'])+ str(cur['properties']['COUNTYFP10'])+ str(cur['properties']['TRACTCE10'])+ str(cur['properties']['BLOCKCE10'])
        place_pop_data.append(filter_by_property(matchstr, 'BLOCKID10', county_bg_data)[0])

    temp_place_pop_data = remove_duplicate_blocks(place_pop_data)
    place_pop_data = temp_place_pop_data

    t_place_pop_geoms = []
    for bg in place_pop_data:
        t_place_pop_geoms.append(transform(proj2, geometry.shape(bg['geometry'])))

    place_geoms = MultiPolygon(t_place_pop_geoms)
    
    return {"place_geom": place_geoms, "place_pop_data": place_pop_data }




def get_place_blockgroup_data_old(bg_data_str, places_data_str, faces_data_str, countycode, placename):

    state_bg_data = fiona.open(bg_data_str) # this has the pop. data
    state_places_data = fiona.open(places_data_str)
    county_faces_data = fiona.open(faces_data_str)

    # filter data (inc. pop.) to county of interest
    county_bg_data = filter_by_property(countycode, 'COUNTYFP10', state_bg_data)
    
    place_bg_data = filter_by_property(placename, 'NAME10', state_places_data)

    # get the shapely geometry of the BGs for place
    place_geoms = geometry.shape(place_bg_data[0]['geometry'])
    # turn it into a multipolygon for easier processing
    place_geoms_poly = MultiPolygon([place_geoms])

    # match the PLACECODE to the data in the FACES file to get place BG data
    placecode = county_faces_data[0]['properties']['PLACEFP10']

    place_faces_data = filter_by_property(placecode,'PLACEFP10', county_faces_data)

    # now we have all the BG info we need to get the populations in PLACE
    place_pop_data = []

    for i in range(0, len(place_faces_data)):
        cur = place_faces_data[i]
        matchstr = str(cur['properties']['STATEFP10'])+ str(cur['properties']['COUNTYFP10'])+ str(cur['properties']['TRACTCE10'])+ str(cur['properties']['BLOCKCE10'])
        place_pop_data.append(filter_by_property(matchstr, 'BLOCKID10', county_bg_data)[0])

    # remove BG's with duplicate BLOCKIDs
    temp_place_pop_data = remove_duplicate_blocks(place_pop_data)
    place_pop_data = temp_place_pop_data

    # turn coords into actual shapley shapes
    place_pop_geom_data  = []
    for i in range(len(temp_place_pop_data)):
        place_pop_geom_data.append(geometry.shape(temp_place_pop_data[i]['geometry']))

    place_pop_geom_data = MultiPolygon(place_pop_geom_data)

    #we have to transform everything to meters so we can operate on the coords
    t_place_pop_geom_data = []
    for i in range(len(place_pop_geom_data)):
    	t_place_pop_geom_data.append(transform(proj2, place_pop_geom_data[i]))

    t_place_pop_geom_data = MultiPolygon(t_place_pop_geom_data)

    return {"trans_place_pop_geom_data": t_place_pop_geom_data,
    		"place_pop_geom_data": place_pop_geom_data,
    		"place_pop_data": place_pop_data}


   



