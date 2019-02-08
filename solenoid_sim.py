# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 23:09:27 2019

@author: Dylan
"""

import vpython as vp
import numpy as np


class Particle():
    '''A particle moving in magnetic field.
    
    Attributes:
        x: 3-vector position of particle
        v: 3-vector velocity of particle
        charge: charge of particle
        mass: mass of particle
    '''
    
    def __init__(self, x0, v0, charge, mass):
        self.x = np.array(x0)
        self.v = np.array(v0)
        self.charge = charge
        self.mass = mass
        
        


def main():
    pars, consts, objs = init()
    simulate(pars, consts, objs)
    print('donzo')
    
    
    
def init():
    pars = initPars()
    consts = initConsts()
    objs = initObjs(pars, consts)
    
    return(pars, consts, objs)


def initPars():
    pars = {'visual':True,
            't0':0.0, #s Start time
            'tf':1.0e-8, #s End time
            'dt':1.0e-14, #s Time step
            'ballR':0.1, #cm Radius of balls representing electrons.
            'ballColors':[vp.color.yellow, vp.color.green], #vp.color Color of electrons.
            'ballTrail':True, #bool Make trail behind electrons if true.
            'd':0.24, #cm Radius of detector opening
            'L':106.4, #cm Distance from source to detector opening.
            'B':np.array([0.0, 0.0, 5.027*11]), #gauss Magnetic field, assumed constant and homogenous.
            'polarAngles':[12.0, 12.0], #degrees Each angle will initialize electron with given angle as initial velocity relative to B-field axis.
            'azimuthAngles':[0.0, 180.0], #degrees 
            'energies':np.array([1.0e-7, 2.0e-8]) + 8.1871053618112e-7, #ergs Each energy will initialize electron with given energy. Rest mass added.
            'sourcePos':np.array([0.0, 0.0, 0.0]) #cm Location of source where particles originate
            }
    
    return(pars)


def initConsts():    
    consts = {'e':4.80320427e-10, #Fr Electron charge in weird Gaussian units.
              'c':2.99792458e10, #cm/s Speed of light
              'm_e':9.10938e-28 #g Mass of electron
              }
    
    return(consts)



def initObjs(pars, consts):
    particles = []
    balls = []
    for i in range(len(pars['polarAngles'])):
        x = pars['sourcePos']
        charge = -consts['e']
        mass = consts['m_e']
        v_mag = vFromP(pFromE(pars['energies'][i], mass, consts['c']), mass, consts['c'])
        v = calcV(pars['polarAngles'][i], pars['azimuthAngles'][i], v_mag)
        
        particles.append(Particle(x, v, charge, mass))
        balls.append(vp.sphere(radius=pars['ballR'], make_trail=pars['ballTrail'], color=pars['ballColors'][i]))
        
    
    detector = vp.cylinder(pos=vp.vector(0,0,pars['L']), axis=vp.vector(0,0,5), radius=pars['d'])
    source = vp.box(pos=vp.vector(0,0,0), size=vp.vector(1,1,1))
    axis = vp.cylinder(pos=vp.vector(0,0,0), axis=vp.vector(0,0,pars['L']), radius=0.1, color=vp.color.red)
        
    vp.scene.camera.follow(balls[0])
        
    objs = {'particles':particles,
            'balls':balls,
            'detector':detector,
            'source':source,
            'axis':axis}
        
    return(objs)



def simulate(pars, consts, objs):
    t = pars['t0']
    while(t <= pars['tf']):
        for i in range(len(objs['particles'])):
            
            particle = objs['particles'][i]
            v = particle.v
#            print(v)
            m = particle.mass
            E = eFromP(pFromV(np.linalg.norm(v), m, consts['c']), m, consts['c'])
#            print(E)

#            print(particle.x)
            
            dx = v * pars['dt']
            particle.x += dx
            
            dv = calcDV(v, E, pars['B'], consts['c'], particle.charge, pars['dt'])
#            print(dv)
            particle.v += dv
            
            objs['balls'][i].pos = vp.vector(particle.x[0], particle.x[1], particle.x[2])
            
        t += pars['dt']


def calcV(polar, azimuth, v_mag):
    polar = radFromDeg(polar)
    azimuth = radFromDeg(azimuth)
    x = v_mag * np.sin(polar) * np.cos(azimuth)
    y = v_mag * np.sin(polar) * np.sin(azimuth)
    z = v_mag * np.cos(polar)
    
    return(np.array([x, y, z]))


def calcDV(v, E, B, c, charge, dt):
#    print(v, E, B, c, charge, dt)
    return(charge * c / E * np.cross(v, B) * dt)

def pFromE(E, m, c):
    return((E**2 / c**2 - m**2 * c**2)**0.5)

def vFromP(p, m, c):
    return((1.0 / (1.0 / c**2 + m**2 / p**2))**0.5)

def pFromV(v, m, c):
    return(m * v / (1.0 - v**2 / c**2)**0.5)

def eFromP(p, m, c):
    return(c * (p**2 + m**2 * c**2)**0.5)

def radFromDeg(deg):
    return(deg / 180.0 * np.pi)




if(__name__ == '__main__'):
    main()