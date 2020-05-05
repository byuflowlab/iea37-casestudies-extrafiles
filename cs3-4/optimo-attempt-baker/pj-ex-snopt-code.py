from pyoptsparse import Optimization, SNOPT, pyOpt_solution, NSGA2
import numpy as np
import matplotlib.pyplot as plt
from position_constraints import *
from windRoses import *
from grid_param_test import *
from aep_calc import *
import scipy as sp
import os
import sys
sys.dont_write_bytecode = True
def obj_func_grid(xdict):
    global rotorDiameter
    global turbineZ
    global windDirections
    global windSpeeds
    global windFrequencies
    global shearExp
    global minSpacing
    global nTurbs
    global boundaryVertices
    global boundaryNormals
    global circle_radius
    global rf
    global wakemodel
    global nRows
    global turbs_per_row
    global x_start
    global anchor_x
    global anchor_y
    global nCalls
    nCalls += 1
    dx = xdict['dx']
    dy = xdict['dy']
    offset = xdict['offset']
    rotate = xdict['rotate']
    # scale = xdict['scale']
    scale = 1.0
    turbineX,turbineY = makeGrid_centered(dx,dy,offset,rotate,turbs_per_row,x_start)
    if nCalls == 1:
        print turbineX
        print turbineY
    show = False
    if show == True:
        plt.figure(1)
        plt.clf()
        for i in range(nTurbs):
            circ = plt.Circle((turbineX[i],turbineY[i]), rotorDiameter[i]/2.,facecolor="blue",edgecolor="blue",alpha=0.2)
            plt.gca().add_patch(circ)
        # circ = plt.Circle((0.,0.), circle_radius,facecolor="None",edgecolor="black",alpha=0.8)
        # plt.plot(np.array([-1800.,-1800.,1800.,1800.,-1800]),np.array([-1800.,1800.,1800.,-1800.,-1800.]),'-k')
        # plt.gca().add_patch(circ)
        # plt.xlim(-2200.,2200.)
        xb = boundaryVertices[:,0]
        yb = boundaryVertices[:,1]
        xb = np.append(xb,xb[0])
        yb = np.append(yb,yb[0])
        plt.plot(xb,yb,'--k')
        plt.axis('equal')
        plt.axis('off')
        plt.draw()
        plt.pause(0.001)
    funcs = {}
    AEP = fast_calc_AEP(turbineX, turbineY, turbineZ, rotorDiameter, windDirections,
                windSpeeds, windFrequencies, wakemodel=wakemodel,relaxationFactor=rf)
    print -AEP/1.E5
    funcs['obj'] = -AEP/1.E5
    funcs['sep'] = SpacingConstraint(turbineX, turbineY, rotorDiameter, minSpacing=minSpacing)/1.E5
    bounds = arbitraryBoundary(turbineX, turbineY, boundaryVertices, boundaryNormals)/1.E3
    b = np.zeros(np.shape(bounds)[0])
    for i in range(len(b)):
        b[i] = min(bounds[i])
    funcs['bound'] = b
    fail = False
    return funcs, fail
## Main
if __name__ == "__main__":
    global rotorDiameter
    global turbineZ
    global windDirections
    global windSpeeds
    global windFrequencies
    global shearExp
    global minSpacing
    global nTurbs
    global boundaryVertices
    global boundaryNormals
    global circle_radius
    global rf
    global wakemodel
    global nRows
    global turbs_per_row
    global x_start
    global anchor_x
    global anchor_y
    global nCalls
    wakemodel = "gaussian"
    nTurbs = 100
    rose = 'ukiah'
    windDirections, windFrequencies, windSpeeds = ukiahRose(30)
    wind_angle = windDirections[np.argmax(windFrequencies)]
    windDirections, windFrequencies, windSpeeds = ukiahRose(30,nSpeeds=8)
    windDirections -= wind_angle
    turbineZ = np.ones(nTurbs)*100.
    rotorDiameter = np.ones(nTurbs)*130.
    shearExp = 0.15
    minSpacing = 2.0
    maxAEP = 0.
    spacing = 4.
    side_length = (np.sqrt(nTurbs)-1.)*rotorDiameter[0]*spacing
    a = side_length**2
    circle_radius = np.sqrt(a/np.pi)
    # folder = 'grid_square_%s_%s_%s_%s'%(nTurbs,spacing,rose,wakemodel)
    # folder = 'grid4_square_%s_%s_%s_%s'%(nTurbs,spacing,rose,wakemodel)
    #
    # if not os.path.exists(folder):
    #     os.makedirs(folder)
    """circle boundary"""
    nBounds = 20
    # circle_radius = 5280.
    xBounds = np.zeros(nBounds)
    yBounds = np.zeros(nBounds)
    theta = np.linspace(0.,2.*np.pi-2.*np.pi/float(nBounds),nBounds)
    for i in range(nBounds):
        xBounds[i] = circle_radius*np.cos(theta[i])
        yBounds[i] = circle_radius*np.sin(theta[i])
    x = np.zeros_like(xBounds)
    x[:] = xBounds[:]
    y = np.zeros_like(yBounds)
    y[:] = yBounds[:]
    xBounds = x*np.cos(np.deg2rad(wind_angle)) - y*np.sin(np.deg2rad(wind_angle))
    yBounds = x*np.sin(np.deg2rad(wind_angle)) + y*np.cos(np.deg2rad(wind_angle))
    """square boundary rotated 30 deg from dominant wind direction"""
    # nBounds = 4
    # x = np.array([-side_length/2.,side_length/2.,side_length/2.,-side_length/2.])
    # y = np.array([-side_length/2.,-side_length/2.,side_length/2.,side_length/2.])
    # xBounds = x*np.cos(np.deg2rad(30.)) - y*np.sin(np.deg2rad(30.))
    # yBounds = x*np.sin(np.deg2rad(30.)) + y*np.cos(np.deg2rad(30.))
    # locations = np.zeros((nBounds,2))
    # locations[:, 0] = xBounds
    # locations[:, 1] = yBounds
    # boundaryVertices, boundaryNormals = calculate_boundary(locations)
    """amalia boundary"""
    # locations = np.loadtxt('layout_amalia.txt')
    # xBounds = locations[:, 0]
    # yBounds = locations[:, 1]
    # xBounds = xBounds - min(xBounds) - (max(xBounds)-min(xBounds))/2.
    # yBounds = yBounds - min(yBounds) - (max(yBounds)-min(yBounds))/2.
    # locations[:, 0] = xBounds
    # locations[:, 1] = yBounds
    # boundaryVertices, boundaryNormals = calculate_boundary(locations)
    # xBounds = boundaryVertices[:, 0]
    # yBounds = boundaryVertices[:, 1]
    nBounds = len(xBounds)
    points = np.zeros((nBounds,2))
    points[:, 0] = xBounds
    points[:, 1] = yBounds
    # hull = sp.spatial.ConvexHull(points)
    # area = hull.volume
    # area_ratio = area/(np.pi*circle_radius**2)
    # xBounds = xBounds/np.sqrt(area_ratio)
    # yBounds = yBounds/np.sqrt(area_ratio)
    # points[:, 0] = xBounds
    # points[:, 1] = yBounds
    # boundaryVertices, boundaryNormals = calculate_boundary(points)
    xBounds = np.append(xBounds,xBounds[0])
    yBounds = np.append(yBounds,yBounds[0])
    dx_start,dy_start,offset_start,rotate_start,turbs_per_row,x_start = make_start_grid_test(nTurbs,boundaryVertices,boundaryNormals)
    factors = np.array([3.0,2.75,2.5,2.25,2.0,1.75,1.5,1.25,1.0])
    rf = 1.0
    num = 1
    nCalls = 0
    for i in range(num):
        print 'iterations: ', i
        # dx = float(np.random.rand(1))*dx_start+200.
        # dy = float(np.random.rand(1))*dy_start+200.
        # offset = float(np.random.rand(1))*25.+10.
        # rotate = float(np.random.rand(1))*30.-15.
        dx = dx_start
        dy = dy_start
        offset = offset_start
        rotate = rotate_start
        input = {'dx':dx,'dy':dy,'offset':offset,'rotate':rotate}
        funcs,_ = obj_func_grid(input)
        AEPstart = funcs['obj']
        print AEPstart
        nCalls = 0
        """Optimization"""
        optProb = Optimization('Wind_Farm_AEP', obj_func_grid)
        optProb.addObj('obj')
        optProb.addVar('dx', type='c', lower=0., upper=None, value=dx)
        optProb.addVar('dy', type='c', lower=0., upper=None, value=dy)
        optProb.addVar('offset', type='c', lower=None, upper=None, value=offset)
        optProb.addVar('rotate', type='c', lower=None, upper=None, value=rotate)
        num_cons_sep = (nTurbs-1)*nTurbs/2
        optProb.addConGroup('sep', num_cons_sep, lower=0., upper=None)
        optProb.addConGroup('bound', nTurbs, lower=0., upper=None)
        opt = SNOPT()
        opt.setOption('Scale option',0)
        opt.setOption('Iterations limit',1000000)
        opt.setOption('Summary file','summary_grid.out')
        opt.setOption('Major optimality tolerance',1.e-5)
        opt.setOption('Major feasibility tolerance',1.e-6)
        res = opt(optProb)
        dx_f = res.xStar['dx']
        dy_f = res.xStar['dy']
        o_f = res.xStar['offset']
        r_f = res.xStar['rotate']
        # print 'dx_f: ', dx_f
        # print 'dy_f: ', dy_f
        # print 'offset_f: ', o_f
        # print 'rotate_f: ', r_f
        print nCalls
        input = {'dx':dx_f,'dy':dy_f,'offset':o_f,'rotate':r_f}
        funcs,_ = obj_func_grid(input)
        separation = min(funcs['sep'])
        boundary = min(funcs['bound'])
        if separation > -1.E-4 and boundary > -1.E-4 and -funcs['obj'] > maxAEP:
            maxAEP = -funcs['obj']
        if separation > -1.E-4 and boundary > -1.E-4:
            print 'AEP opt: ', -funcs['obj']
            # file = open('%s/AEP.txt'%folder, 'a')
            # file.write('%s'%(-funcs['obj']) + '\n')
            # file.close()
            #
            # file = open('%s/funcCalls.txt'%folder, 'a')
            # file.write('%s'%nCalls + '\n')
            # file.close()
        print 'maxAEP: ', maxAEP
    print 'final maxAEP: ', maxAEP
    plt.show()