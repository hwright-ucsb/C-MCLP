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


def gdf_to_list(gdf):
    tmp = []
    for key,item in gdf['geometry'].iteritems():
        tmp.append(item.coords[0])
        
    return tmp

def list_of_gdf_to_list_of_coords(lgdf):
	tmp = []
	for l in lgdf:
	    if len(l)>0:
	        for key,item in l['geometry'].iteritems():
	            tmp.append(item.coords[0])
	return tmp
            

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
        
    return pts_per_bg


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

#make sure all the sums and pops are equal
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


################################################################
################################################################
################################################################
				# GRID PT GENERATION #


def get_largest_diameter(bds):
	# get the largest diagonal dist and this is what we will make 
	#our regtangle of points on the grid based on (so when we shift/rotate the whole 
	#shape will be covered still)

	[west, south, east, north] = bds
	bounds = [west,south,east,north]
	diags = [db.eucdist(west,north,east,south)[0], db.eucdist(west,south,east,north)[0]]
	maxdiag = max(diags)
	return maxdiag

def gen_triangular_gridpts(bds, rad):
	#generate the "two" grids

	bounds = bds
	[west,south,east,north] = bounds
	maxdiag = get_largest_diameter(bds)
	r = rad
	x_shift = (r/2.0)*math.sqrt(3)
	y_shift = (3.0/2.0)*r

	x_ints = int(math.ceil(maxdiag/x_shift))
	y_ints = int(math.ceil(maxdiag/y_shift))

	stopx = bounds[0]+(x_ints*x_shift)
	stopy = bounds[1]+(y_ints*y_shift)

	# create the equally spaced coord points for the "grid" cand pts
	#xs_all = np.linspace(west,stopx,x_ints)
	#ys_all = np.linspace(south,stopy,y_ints)
	xs_all = []
	ys_all = []
	cur_x = south
	cur_y = west
	for i in range(0,x_ints):
		xs_all.append(cur_x)
		cur_x = cur_x+x_shift
	for i in range(0,y_ints):
		ys_all.append(cur_y)
		cur_y = cur_y+y_shift


	# now take odds with odds and evens with evens ("2" grids)

	xs_odd = []
	xs_even = []
	ys_odd = []
	ys_even = []

	for i in range(0,len(xs_all)):
	    if i%2 == 0:
	        xs_even.append(xs_all[i])
	    else:
	        xs_odd.append(xs_all[i])
	for i in range(0,len(ys_all)):
	    if i%2 == 0:
	        ys_even.append(ys_all[i])
	    else:
	        ys_odd.append(ys_all[i])

	xs_odd = np.array(xs_odd)
	xs_even = np.array(xs_even)
	ys_odd = np.array(ys_odd)
	ys_even = np.array(ys_even)

	x_grid1, y_grid1 = np.meshgrid(xs_odd, ys_odd)
	x_grid2, y_grid2 = np.meshgrid(xs_even, ys_even)
	# create the data dicts that will create the geodataframe
	# create the data dicts that will create the geodataframe
	data1={'x':x_grid1.flatten(), 'y':y_grid1.flatten()}
	data2={'x':x_grid2.flatten(), 'y':y_grid2.flatten()}
	data_all = data1.copy()   # start with x's keys and values
	data_all.update(data2)
	#grid1gdf = gpd.GeoDataFrame(data={'y':x_grid1.flatten(), 'x':y_grid1.flatten()})
	#grid2gdf = gpd.GeoDataFrame(data={'y':x_grid2.flatten(), 'x':y_grid2.flatten()})
	grid1gdf = gpd.GeoDataFrame(data={ 'x':y_grid1.flatten(), 'y':x_grid1.flatten()})
	grid2gdf = gpd.GeoDataFrame(data={'x':y_grid2.flatten(), 'y':x_grid2.flatten()})

	gridgdf = gpd.GeoDataFrame( pd.concat( [grid1gdf,grid2gdf], ignore_index=True) )
	# transform the grid pts for indexing
	gridgdf['geometry'] = gridgdf.apply(lambda row: Point((row['x'], row['y'])), axis=1)

	#create the grid index
	grid_index = gridgdf.sindex

	#possible_matches_index = list(grid_index.intersection(bounds))
	#possible_matches = gridgdf.iloc[possible_matches_index]
	#precise_matches = possible_matches[possible_matches.intersects(wc_census10)]
	#points_within_geometry = precise_matches
	#points_outside_geometry = gridgdf[~gridgdf.isin(points_within_geometry)]

	return gridgdf




