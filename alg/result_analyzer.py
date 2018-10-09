############################################
############ result_analyzer.py ############
############################################

import shapely
from shapely import geometry
from shapely.ops import transform
from shapely.ops import cascaded_union
from shapely.geometry import Point, Polygon, MultiPolygon, LineString
from descartes import PolygonPatch
import sys
sys.path.append('../../../alg')
import distance_buffer as db, census_parser as cp, cmclp_solver as solver, data_generator as datagen


from rtree import index

import math
import pyproj
import fiona

import matplotlib.pyplot as plt, pandas as pd, geopandas as gpd, numpy as np

from matplotlib import pyplot
from functools import partial

from matplotlib.patches import Arc
from scipy.stats import truncnorm


def generate_covering_circles(r,coords):
    return [Point(coord).buffer(r) for coord in coords]



def site_indexes_to_coords(results, site_coords):
    """
    Takes the stage result dicts and  gets their coordinate lists
    
    Parameters
    ----------
      results: list[dict['coveringsites']] int
      site_coords: list[(float,float)] ; list of all potl facility sites

    Returns
    -------
      coords: list[list[(float,float)]]
      
    Note: 
      - eliminates empty stages also
    """
    coords = []
    for i in range(0,len(results)):
        cursites = results[i]['coveringsites']
        if len(cursites) > 0: #clean out empty end stages
            tmpcoords = []
            for ind in cursites:
                tmpcoords.append(site_coords[ind])
            coords.append(tmpcoords)
    return coords


def coalesce_geoms(list_of_geoms):
    one_shape = list_of_geoms[0]
    for i in range(1,len(list_of_geoms)):
        one_shape = one_shape.union(list_of_geoms[i])
    return one_shape


def calculate_overlap(circles):
    """
    
    Parameters
    ----------
      circles: 

    Returns
    -------
      intersection: 
      
    Note: 
      - src: https://gis.stackexchange.com/questions/237053/getting-intersection-of-multiple-polygons-efficiently-in-python?rq=1
      - USES RTREE PKG
    """
    intersections = []
    idx = index.Index()
    for pos, circle in enumerate(circles):
        idx.insert(pos, circle.bounds)

    for circle in circles:
        merged_circles = cascaded_union([circles[pos] for pos in idx.intersection(circle.bounds) if circles[pos] != circle])
        intersections.append(circle.intersection(merged_circles))

    intersection = cascaded_union(intersections)
    return intersection


def generate_graph_data(mclpdata, cmclpdata, filterstr):
    """
    
    Parameters
    ----------
      mclpdata: list[dict{...}]
      cmclpdata: list[dist{...}]
      filterstr: str

    Returns
    -------
      lines: list[list[]]
      
    Note: 
    """
    lines = []
    tline = []
    for i in range(0,len(mclpdata)):
        tline.append(mclpdata[i][filterstr])
    lines.append(tline)
    tline = []
    for i in range(0,len(cmclpdata)):
        tline.append(cmclpdata[i][filterstr])
    lines.append(tline)
    return lines
        

 ### PLOTTING METHODS

def plot_data(ax, numstages, data, styles):
    for i in range(0,len(data)):
        ax.plot(list(np.arange(numstages)),data[i][0:numstages], styles[i])

def plot_overlap(ax, numstages, data, styles=['r-','b*','go']):
    """
    Plots overlapping area in ?? meters^2 ?? for each stage
    
    Parameters
    ----------
      ax : matplotlib axis
      numstages : int
      data : list[list[float]]
          List of data to plot for each stage. len(data) == number of lines
      style: list[str]
          
    Returns
    -------
    """
    plot_data(ax,numstages,data,styles)
    


def plot_coverage(ax, numstages, data, styles=['r-','b*','go']):
    """
    Plots ACTUAL and CALC'D pop. coverage for each stage
    
    Parameters
    ----------
      ax : matplotlib axis
      numstages : int
      data : list[list[int]]
          List of data to plot for each stage. len(data) == number of lines
      style: list[str]
          
    Returns
    -------
    """
    plot_data(ax,numstages,data,styles)


def plot_geoms(geoms):
    fig1, ax1 = plt.subplots(figsize=(10,10))
    for i in range(0,len(geoms)):
        xs, ys = geoms[i].exterior.xy
        ax1.fill(xs, ys, alpha=0.5, fc='r', ec='none')
    pyplot.axis('scaled')
    plt.show() 


