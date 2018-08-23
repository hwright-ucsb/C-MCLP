# THIS SOLVER USES THE GUROBI OPTIMIZATION SOFTWARE

from gurobipy import *

import shapely
import math
import pyproj
import sys
import fiona
sys.path.append('../../../alg')

import matplotlib.pyplot as plt, pandas as pd, geopandas as gpd, numpy as np
import distance_buffer as db, border_generators as bg, marching_army as ma

from matplotlib import pyplot
from functools import partial
from shapely import geometry
from shapely.ops import transform
from descartes import PolygonPatch
from matplotlib.patches import Arc
from shapely.geometry import Point, Polygon, MultiPolygon, LineString
from scipy.stats import truncnorm

def todict(keys, vals):
    return dict(zip(keys,vals))

def clean_stage_solns(stagesoln, facilitypt_coords):
    site_coords = []
    for stage in stagesoln:
        if stage['coveredpop'] >0:
            tmp_coords = []
            for ind in stage['coveringsites']:
                tmp_coords.append(facilitypt_coords[ind])
            site_coords.append(tmp_coords)
            
    return site_coords

#might not need this
def distance(a,b):
    dx = a[0] - b[0]
    dy = a[1] - b[1]
    return math.sqrt(dx*dx + dy*dy)

def generate_set_N(demandSites, potlFacilitySites, S, tol):
		# TODO: need to use rtree or something to find the N set(might have done?)
	# N[i] is the covering set of i (all j that can cover i)
	numDemandSites = len(demandSites)
	numPotlFacilitySites = len(potlFacilitySites)
	maxdist = S+(S*tol)

	N = {}
	for i in range(numDemandSites):
	    N[i] = []
	    for j in range(numPotlFacilitySites):
	        if distance(demandSites[i], potlFacilitySites[j]) <= maxdist:
	            N[i].append(j)

	return N

def cmclp_solve(facsites, dmdpts, dmdwts, N, P, S, objparams, relgaptol=0.0004):
	output = {}

	#Problem Data
	demandSites = dmdpts
	demandWt = dmdwts 
	potlFacilitySites = facsites 

	numDemandSites = len(demandSites)
	numPotlFacilitySites = len(potlFacilitySites)
	filestr = "DMD-"+str(numDemandSites)+"-FAC-"+str(numPotlFacilitySites)+".log"

	m = Model('C-MCLP-v0.0')
	#m.setParam(GRB.Param.OutputFlag, 1)
	m.setParam(GRB.Param.LogToConsole, 0)
	m.setParam(GRB.Param.LogFile, filestr)
	
	#SET GAP TOLERANCE
	m.setParam(GRB.Param.MIPGap, relgaptol)

	# X : POTL FACILITY SITES
	x = m.addVars(range(numPotlFacilitySites), vtype=GRB.BINARY, name="x")
	e = m.addVars(range(numPotlFacilitySites), vtype=GRB.BINARY, name="e")
	y = m.addVars(range(numDemandSites), vtype=GRB.CONTINUOUS, name="y")

	# Add Constraints

	#TODO - not sure what this does/where it should go
	m.update()

	for i in range(numDemandSites):
	    m.addConstr(quicksum(x[j] + y[i] for j in N[i]) >= 1) #c1
	    m.addConstr(quicksum(x[j] + e[j] for j in N[i]) >= 1) #c3
	    
	m.addConstr(x.sum() == P) #c2

	for i in range(numDemandSites): #c6
	    m.addConstr(y[i] >= 0)

	# Set global sense for ALL objectives
	m.ModelSense = GRB.MINIMIZE

	# Set Objectives
	    #expr: New alternative objective.
	    #index: Index for new objective. If you use an index of 0, this routine will change the primary optimization objective.
	    #priority: Priority for the alternative objective.
	    #weight: Weight for the alternative objective.
	    #abstol: Absolute tolerance for the alternative objective. 
	    #reltol: Relative tolerance for the alternative objective. 
	    #name: Name of the alternative objective. 

	    
	#o1- min amt of demand that is NOT covered
	m.setObjectiveN(quicksum(y[i]*demandWt[i] for i in range(numDemandSites)),
	                    objparams["index"][0], objparams["priority"][0], objparams["weight"][0], 
	                    reltol=objparams["reltol"][0])#, abstol[0], reltol[0])

	#o2 - (while) min'g w/e you need to add to cover everything
	m.setObjectiveN(quicksum(e[j] for j in range(numPotlFacilitySites)),
	                    objparams["index"][1], objparams["priority"][1], objparams["weight"][1], 
	                    reltol=objparams["reltol"][1])#, abstol[1], reltol[1])

		#%tb
	siteschosen = []
	try:
	    m.optimize()

	    # Status checking
	    status = m.Status
	    if status == GRB.Status.INF_OR_UNBD:
	        output["status"]= "Inf or unbd"
	        return [], {}
	    elif status == GRB.Status.INFEASIBLE:
	        output["status"]= "infeasible"
	        return [], {}
	    elif status == GRB.Status.UNBOUNDED:
	        output["status"]= "unbounded"
	        return [], {}
	    elif status == GRB.Status.OPTIMAL:
	        output["status"]= "optimal"
	    else:
	        output["status"]= str(status)
	        return [], {}
	        
	    # record chosen sites (soln)
	    for j in range(numPotlFacilitySites):
	        if x[j].X > 0.9:
	            siteschosen.append(j)
	        

	    #  number of solutions 
	    nSolutions = m.SolCount
	    output["numsolns"] = nSolutions

	    # Print objective values of solutions
	    if nSolutions > 10:
	        nSolutions = 10
	    output["objvals"] = []
	    for i in range(2): #num of objectives
	        m.setParam(GRB.Param.ObjNumber, i)
	        output["objvals"].append([])
	        for e in range(nSolutions):
	            m.setParam(GRB.Param.SolutionNumber, e)
	            output["objvals"][i].append(m.ObjNVal)

		output["errorstr"] = "NONE"
	except GurobiError as e:
	    output["errorstr"] = str(e.errno) + ": " + str(e)
	except AttributeError as e:
	    output["errorstr"] = 'ATTR ERR: ' + str(e)

	return siteschosen, output


def cmclp_stages(demandpts, demandwts, facilitypts, stageP, stageS, stopcond=2, disttol=0.0, objparams={"index": [0,1], "priority":[0,0], "weight":[1.,1.],
                      "abstol":[0.,0.], "reltol":[0.,0.]}, relgaptol=0.0004):
	# TODO: find out  - how do we get the chosen sites from the nonopt solns?
	#       (i.e. X # solns found...) - also, do we care? (why would we?)

	curdemandpts = demandpts
	curdemandwts = demandwts
	curpotlsites = facilitypts
	curnumdemandpts = len(curdemandpts)
	pt_wt_dict = todict(curdemandpts, curdemandwts)

	P = stageP
	S = stageS
	objectiveparams = objparams

	# {"sites":, "coveredpts":, "coveredpop":, "output":, "coverror":}
	stagesolns = []

	# TESTING PURPOSES
	#k = 0
	#cnt # times we've had the same # demand pts; stop after STOP
	repeatdemandptcnt = 0
	stop = stopcond
	tol = disttol

	while (repeatdemandptcnt < stop):
	    print len(curdemandpts)
	    print len(curpotlsites)
	    cur_N = generate_set_N(curdemandpts, curpotlsites, S, tol)

	    cursoln, curstatus = cmclp_solve(curpotlsites, curdemandpts, curdemandwts, cur_N, P, S, objectiveparams, relgaptol)
	    
	    #return the ACTUAL covered demand pts and the ACTUAL covered pop
	    curcoveredpts, coveringsites, curcoveredpop = actual_coverage(cursoln, cur_N, curdemandwts)
	    #these are the demand pts NOT covered by the soln
	    coverage_error = list(set(cursoln)-set(curcoveredpts))
	    
	    # we will use the the ACTUALLY covered pts, get those coords
	    tmp_demandpts = []
	    for ind in curcoveredpts:
	        tmp_demandpts.append(curdemandpts[ind])
	        
	    nextdemandpts = list(set(curdemandpts)-set(tmp_demandpts))
	    
	    nextdemandwts = []
	    for pt in nextdemandpts:
	        nextdemandwts.append(pt_wt_dict[pt])
	        
	    stagesolns.append({"sites":cursoln, "coveredpts":curcoveredpts, "coveringsites": coveringsites,
	                       "coveredpop":curcoveredpop, "output":curstatus, "coverror":coverage_error})

	    if len(curdemandpts) == len(nextdemandpts):
	        repeatdemandptcnt += 1
	    else:
	        repeatdemandptcnt = 0
	        
	    curdemandwts = nextdemandwts
	    curdemandpts = nextdemandpts

	return stagesolns
    

def actual_coverage(chosen_sites, N, wts):
    #map coverage to a dict for easy accessing 
    #chosen_sites = siteschosen
    #wts = tmp_wts
    #{facility_site_index}: [list of demand pts covered (index)]
    if not len(chosen_sites) >0:
    	return [],0

    tmpdict = {}
    for i in range(0,len(N)):
        tmp = N[i]
        for j in range(0,len(tmp)):
            if tmp[j] not in tmpdict.keys():
                tmpdict[tmp[j]] = [i]
            else:
                tmpdict[tmp[j]].append(i)

    covered_pts = []
    covering_sites = []
    for i in range(0,len(chosen_sites)):
    	if chosen_sites[i] in tmpdict:
	        tmp_cov = tmpdict[chosen_sites[i]]
	        covering_sites.append(chosen_sites[i])

	        for j in range(0,len(tmp_cov)):
	            if tmp_cov[j] not in covered_pts:
	                covered_pts.append(tmp_cov[j])

    cov_pop = 0
    cov_wts = []
    for i in range(0,len(covered_pts)):
        curwt = wts[covered_pts[i]]
        cov_pop += curwt
        cov_wts.append(curwt)

    return covered_pts, covering_sites, cov_pop


def mclp_solve(tmp_facilitypts, tmp_demandpts, tmp_wts, N, P, S):
    
    #Problem Data

    output = {}
    demandSites = tmp_demandpts
    demandWt = tmp_wts 
    potlFacilitySites = tmp_facilitypts 

    numDemandSites = len(demandSites)
    numPotlFacilitySites = len(potlFacilitySites)
    filestr = "MCLP-DMD-"+str(numDemandSites)+"-FAC-"+str(numPotlFacilitySites)+".log"


    m = Model('MCLP-v0.0')
    m.setParam(GRB.Param.LogToConsole, 0)
    m.setParam(GRB.Param.LogFile, filestr)


        # X : POTL FACILITY SITES
    x = m.addVars(range(numPotlFacilitySites), vtype=GRB.BINARY, name="x")
    y = m.addVars(range(numDemandSites), vtype=GRB.BINARY, name="y")

    #TODO - not sure what this does/where it should go
    m.update()

    for i in range(numDemandSites): #c1
        m.addConstr(quicksum(x[j] + y[i] for j in N[i]) >= 1)

    m.addConstr(x.sum() == P) #c2


        #Set global sense for ALL objectives
    m.ModelSense = GRB.MINIMIZE

    m.setObjective(quicksum(y[i]*demandWt[i] for i in range(numDemandSites)))#, abstol[0], reltol[0])


    #%tb
    siteschosen_mclp = []
    try:
        m.optimize()
        m.setParam(GRB.Param.OutputFlag, 0)

        # Status checking
        status = m.Status
        if status == GRB.Status.INF_OR_UNBD:
            output["status"]= "Inf or unbd"
            return [], {}
        elif status == GRB.Status.INFEASIBLE:
            output["status"]= "infeasible"
            return [], {}
        elif status == GRB.Status.UNBOUNDED:
            output["status"]= "unbounded"
            return [], {}
        elif status == GRB.Status.OPTIMAL:
            output["status"]= "optimal"
        else:
            output["status"]= str(status)
            return [], {}

        # Print best selected set
        for j in range(numPotlFacilitySites):
            if x[j].X > 0.9:
                siteschosen_mclp.append(j)


        # Print number of solutions stored
        nSolutions = m.SolCount
        output["numsolns"] = nSolutions

        # Print objective values of solutions
        if nSolutions > 10:
            nSolutions = 10
        output["objvals"] = []
        for i in range(1): #num of objectives
            m.setParam(GRB.Param.ObjNumber, i)
            output["objvals"].append([])
            for e in range(nSolutions):
                m.setParam(GRB.Param.SolutionNumber, e)
                output["objvals"][i].append(m.ObjNVal)

	output["errorstr"] = "NONE"
    except GurobiError as e:
        output["errorstr"] = str(e.errno) + ": " + str(e)

    except AttributeError as e:
        output["errorstr"] = 'ATTR ERR: ' + str(e)

    return siteschosen_mclp, output


def mclp_stages(demandpts, demandwts, facilitypts, stageP, stageS, stopcond=2, disttol=0.0):
    curdemandpts = demandpts
    curdemandwts = demandwts
    curpotlsites = facilitypts
    curnumdemandpts = len(curdemandpts)
    pt_wt_dict = todict(curdemandpts, curdemandwts)

    P = stageP
    S = stageS

    # {"sites":, "coveredpts":, "coveredpop":, "output":, "coverror":}
    stagesolns = []

    #cnt # times we've had the same # demand pts; stop after STOP
    repeatdemandptcnt = 0
    stop = stopcond
    tol = disttol

    while (repeatdemandptcnt < stop):
    
        cur_N = generate_set_N(curdemandpts, curpotlsites, S, tol)
        cursoln, curstatus = mclp_solve(curpotlsites, curdemandpts, curdemandwts, cur_N, P, S)

        #return the ACTUAL covered demand pts and the ACTUAL covered pop
        curcoveredpts, coveringsites, curcoveredpop = actual_coverage(cursoln, cur_N, curdemandwts)
        #these are the demand pts NOT covered by the soln
        coverage_error = list(set(cursoln)-set(curcoveredpts))

        # we will use the the ACTUALLY covered pts, get those coords
        tmp_demandpts = []
        for ind in curcoveredpts:
            tmp_demandpts.append(curdemandpts[ind])

        nextdemandpts = list(set(curdemandpts)-set(tmp_demandpts))

        nextdemandwts = []
        for pt in nextdemandpts:
            nextdemandwts.append(pt_wt_dict[pt])

        stagesolns.append({"sites":cursoln, "coveredpts":curcoveredpts, "coveringsites": coveringsites,
                           "coveredpop":curcoveredpop, "output":curstatus, "coverror":coverage_error})

        if len(curdemandpts) == len(nextdemandpts):
            repeatdemandptcnt += 1
        else:
            repeatdemandptcnt = 0

        curdemandwts = nextdemandwts
        curdemandpts = nextdemandpts

    return stagesolns