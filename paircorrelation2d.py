# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 16:01:32 2020

@author: Antoine Bérut
"""

import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, Point, MultiPoint
from shapely.affinity import translate
import time

def pcf2d(array_positions,bins_distances,coord_border=None,coord_holes=None,fast_method=False,show_timing=False,plot=False,full_output=False):
    r"""
    Computes the 2D radial distribution function (pair correlation function) 
    g(r) for a set of points, corrected to take account of the boundary effects.
    
    Parameters:
        
        - array_positions (numpy array, shape Nx2): 
            the (x,y) coordinates of N points.
                
        - bins_distances (numpy array, shape Mx1): 
            a monotonically increasing array of bin edges defining the values 
            of r for which g(r) is going to be computed.
        
    Optional parameters :
        
        - coord_border (numpy array, shape Lx2): 
            the (x,y) coordinates of L points, defining the boundary enclosing 
            the area of interest to compute the g(r).
            Points in array_positions that are outside the area of interest are 
            automatically excluded.
            !! The list of coordinates must be "valid" in the sense used by the
            shapely library: linking all the points in order should result in a
            simple polygon with no line intersecting each other. !!
            
            Ex: Assuming one wants a square border of area 1
                [[0,0],[0,1],[1,1],[1,0]] is valid (square shape)
                [[0,0],[0,1],[1,0],[1,1]] is not valid (bow tie shape)
                
            If no value is provided, the smallest convex Polygon containing all 
            the points in array_positions is computed using convex_hull. 
            (default value: None)
            
        - coord_holes (list of numpy arrays):
            a list of the (x,y) coordinates of points forming "holes" in the
            area of interest (useful if one is using a geometry with obstacles 
            or exclusion zones where no particles can be found).
            It is possible to define several holes (as long as they do not 
            intersect each other).
            
            Ex: Assuming the area of interest is a square and one wants to 
                to remove a smaller square at the center
                coord_border=[[0,0],[0,1],[1,1],[1,0]]
                coord_holes=[[[0.2,0.2],[0.2,0.8],[0.8,0.8],[0.8,0.2]]]
            
            If no value is provided, the area of intereste will be a simple
            polygon with no hole.
            (default value: None)
        
        - fast_method (bool):
            if True, all the points whose distance to the boundary is less than
            the longest distance in bins_distance are excluded from the g(r)
            computation, and only the points sufficiently far away from the 
            boundary are considered. This method is faster, but might exclude a
            lot of points.
            if False, the code computes for each point its distance to the 
            boundary, then computes a normalization factor for the points that 
            are too close to the boundary. This is the default method that 
            correctly takes account of all points, but might be time consuming.
            (default value: False)
            
        - show_timing (bool):
            if True, the code will print outputs showing timings at different
            stages of the computation (to let one know what operations are the
            most time consuming).
            (default value: False)

        - plot (bool): 
            if True, shows the points that were kept, and the computed g(r). 
            (default value: False).
            
        - full_output (bool):
            if True, the function also returns also "raw" distribution of 
            distances between the points PDF(r), the new array of coordinates 
            with only the points that were considered for computing g(r), the 
            distance of each point to the closest boundary of the area of 
            interest, the normalization factors, and the estimated density of 
            points in the area of interest.
            (default value: False). 
    
    Outputs:
    
        - g(r): a 1x(M-1) numpy array (where M is the length of bin_distances)
        - r: a 1x(M-1) numpy array (where M is the length of bin_distances)
    
    Optional output:
        
        - PDF(r): a 1x(M-1) numpy array
        - new_array_positions: a 2xN numpy array
        - distance_to_boundary: a 1xN numpy array
        - normalization_factor: a Nx(M-1) numpy array
        - estimated_density: a float
        
        (where N in the number of points in the area of interest and M is the 
         length of bin_distances)
    
    """    
    
    if show_timing==True:
        t0=time.time()
    
    if coord_border is None:
        positions_of_interest=MultiPoint(array_positions) #all the points are considered
        area_of_interest=positions_of_interest.convex_hull #the boundary is the convex_hull of all the points
    else:
        if coord_holes is None:
            area_of_interest=Polygon(shell=coord_border) #definition of the boundary
        else:
            area_of_interest=Polygon(shell=coord_border,holes=coord_holes) #definition of the boundary
        
        if not area_of_interest.is_valid:
            print('The list of coordinates your provided for the border is not valid (see help for the definition of valid coord_border).')
            return
        
        #redefinition of "array_positions" with only the points inside the area of interest
        positions_of_interest=area_of_interest.intersection(MultiPoint(array_positions)) #only the points inside the area of interest are considered
        
        array_positions=np.zeros((len(positions_of_interest.geoms),2))
        for i_points in range(len(positions_of_interest.geoms)):
            array_positions[i_points,:]=np.asarray(positions_of_interest.geoms[i_points].coords)
    
    if show_timing==True:
        t1=time.time()-t0
        print('Creating boundary polygon and array of points inside took %f s'%t1)
    
    nb_part=len(array_positions) #number of particles
    densite=nb_part/(area_of_interest.area) #average density of particles
    border_area=area_of_interest.boundary #the boundary (line) of the area of interest (polygon)
    
    rings=[[] for i in range(len(bins_distances)-1)]
    ring_area=np.zeros(len(bins_distances)-1) #true ring areas
    ring_area_approx=np.zeros(len(bins_distances)-1) #approximate ring areas (useful for normalization calculation)
    
    #For each distance bin, defines the ring (difference between the disk of radius r[jj+1] and the disk of radius r[jj]) 
    #and computes the area of those rings.
    for jj in range(len(bins_distances)-1):
        
        inner_circle=Point([0,0]).buffer(bins_distances[jj])
        outer_circle=Point([0,0]).buffer(bins_distances[jj+1])
        rings[jj]=outer_circle.difference(inner_circle)
        
        ring_area_approx[jj]=rings[jj].area
        ring_area[jj]=np.pi*(bins_distances[jj+1]**2-bins_distances[jj]**2)
        
    if show_timing==True:    
        t2=time.time()-t0
        print('Creating all ring polygons took %f s'%(t2-t1))
    
    g_of_r=np.zeros(len(bins_distances)-1)
    g_of_r_normalized=np.zeros(len(bins_distances)-1)
    
    #For each point, computes its distance to the boundary, and for each bin computes the normalization factor
    #(the area of "the intersection of the ring and the area of interest", divided by the ring area)
    normalisation=np.ones((nb_part,len(bins_distances)-1)) #normalization factors to take account of the boundaries
    dist_to_border=np.zeros(nb_part)
    
    if fast_method==True:
        
        for ii in range(nb_part):
            dist_to_border[ii]=positions_of_interest.geoms[ii].distance(border_area) #distance of current point to boundary
        
        far_enough=np.where(dist_to_border>bins_distances[-1])[0] #indexes of points far enough from the boundary
        array_positions=array_positions[far_enough,:] #points too close to the boundary are excluded
        nb_part=len(array_positions) #the new number of points
        
        if full_output==True:
            dist_to_border=dist_to_border[far_enough] #points too close to the boundary are excluded
            normalisation=np.ones((nb_part,len(bins_distances)-1)) #just so that the matrix has the right size
        
    else:
        
        for ii in range(nb_part):
            dist_to_border[ii]=positions_of_interest.geoms[ii].distance(border_area) #distance of point ii to boundary
            
            if dist_to_border[ii]<=bins_distances[0]: #special case the point is too close to the boundary for every r
                for jj in range(len(bins_distances)-1):
                    normalisation[ii,jj]=(area_of_interest.intersection(translate(rings[jj],xoff=positions_of_interest.geoms[ii].xy[0][0],yoff=positions_of_interest.geoms[ii].xy[1][0])).area)/ring_area_approx[jj]
            else:
                for jj in (np.where(bins_distances>dist_to_border[ii])[0]-1): #the normalization factor needs only to be computed for a subset of r
                    normalisation[ii,jj]=(area_of_interest.intersection(translate(rings[jj],xoff=positions_of_interest.geoms[ii].xy[0][0],yoff=positions_of_interest.geoms[ii].xy[1][0])).area)/ring_area_approx[jj]
    
    if show_timing==True:
        t3=time.time()-t0
        print('Computing normalization factors took %f s'%(t3-t2))
    
    #For each point, computes the distance to others, then compute the g(r) by binning
    for ii in range(nb_part):
        #coordinates of the current point
        x_loc=array_positions[ii,0]
        y_loc=array_positions[ii,1]
        
        dist_to_loc=np.sqrt((array_positions[:,0]-x_loc)**2+(array_positions[:,1]-y_loc)**2) #distance of the current point to each other point
        dist_to_loc[ii]=np.inf #the distance from a point to itself is always zero (it is therefore excluded from the computation)
        
        g_of_r=g_of_r+np.histogram(dist_to_loc,bins=bins_distances)[0] #computes the histogram of distances to the current point
        g_of_r_normalized=g_of_r_normalized+np.histogram(dist_to_loc,bins=bins_distances)[0]/(ring_area*normalisation[ii,:]) #computes the histogram of distances to the current point normalized by the area of the intersection of the ring and the area of interest

    if show_timing==True:
        t4=time.time()-t0
        print('Computing g(r) took %f s'%(t4-t3))
    
    g_of_r=g_of_r/nb_part #computes PDF(r)
    g_of_r_normalized=g_of_r_normalized/(nb_part*densite) #computes g(r)
    
    radii=(bins_distances[1::]+bins_distances[0:-1])/2 #computes the values of "r"
    
    if plot==True:
        plt.figure()
        plt.scatter(array_positions[:,0],array_positions[:,1])
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Points kept to compute g(r)')
        
        plt.figure()
        plt.plot(radii,g_of_r_normalized)
        plt.xlabel('r')
        plt.ylabel('g(r)')
        plt.title('Radial Distribution Function')
    
    if full_output==True:
        results=(g_of_r_normalized, radii, g_of_r, array_positions, dist_to_border, normalisation, densite)
    else:
        results=(g_of_r_normalized, radii)
    
    if show_timing==True:
        t5=time.time()-t0
        print('Total time: %f s for %i points '%(t5,nb_part))
    
    return results
