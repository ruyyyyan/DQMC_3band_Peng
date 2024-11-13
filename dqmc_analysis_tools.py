'''
Subroutines for analysing the data from DQMC simulation
'''
import subprocess
import os
import sys
import time
import shutil

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy

from pylab import *
from scipy import *
from scipy import interpolate

import re
import mmap
from linecache import getline

######### writting functions

def write_data_2cols(fname, xs, ys):
    f = open(fname,'w',1) 
    for i in range(len(xs)):
        f.write('{:.6e}\t{:.6e}\n'.format(float(xs[i]),float(ys[i])))

def write_data_3cols(fname, xs, ys, zs):
    f = open(fname,'w',1) 
    for i in range(len(xs)):
        f.write('{:.6e}\t{:.6e}\t{:.6e}\n'.format(float(xs[i]),float(ys[i]),float(zs[i])))

def write_data_4cols(fname, xs, ys, zs, ws):
    f = open(fname,'w',1) 
    for i in range(len(xs)):
        f.write('{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\n'.format(float(xs[i]),float(ys[i]),float(zs[i]),float(ws[i])))

def write_data_add_2cols(fname, x, y):
    f = open(fname,'a',1) 
    f.write('{:.6e}\t{:.6e}\n'.format(float(x),float(y)))
        
def write_data_add_3cols(fname, x, y, z):
    f = open(fname,'a',1) 
    f.write('{:.6e}\t{:.6e}\t{:.6e}\n'.format(float(x),float(y),float(z)))
    
def write_data_TODO(fname, xvar, yvar):
    '''
    NOT finished routine !!!
    write a line of data consisting of xvar, yvar regardless of their length into fname
    '''
    f = open(fname,'a') 
    f.write('{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\n'.format(float(V),float(tperp),float(saf2inf),float(saf1inf)))
    
def write_line2data(fname,line):
    '''
    write a line of any file into data
    '''
    f = open(fname,'a') 
    f.write(line)
    
def find_lines(phrase, ff, dline, Nline, foutname):
    '''
    find the corresonding line x including phrase in ff file, 
    and then starting from line x+dline, write next Nline lines into a file called foutname for further analysis
    '''
    #nums = []
    for i, line in enumerate(ff,1):  # ,1 means acount number from 1 instead of 0 default in python
        if phrase in line:
            #print i, line
            #nums.append(i)
            for j in range(dline, dline+Nline):
                l = getline(ff.name, i+j)
                write_line2data(foutname,l)
    #return nums
      
def Get_dataline(fname, phrase, dataname):
    '''
    Get the line including phrase in file fname and store that line in dataname
    then reopen that temporary file dataname to replace the phrase by empty space
    to get the values
    Note: phrase should be long enough to include comma because it will be replaced
    Normally, dataname is only temporary so can be like './data/tmp'
    '''
    if os.path.isfile(dataname):
        os.remove(dataname)
            
    f = open(fname, 'r')
    find_lines(phrase, f, 0, 1, dataname) 
    
    file = open(dataname, "r")
    text = file.read()
    file.close()
    text = text.replace(phrase, str(' '))
    file = open(dataname, "w")
    file.write(text)
    file.close()
    
    a = loadtxt(dataname, skiprows=0)
    return a

def interpolation(xs, ys, xwant):
    #f = interpolate.Akima1DInterpolator(n_c[:,i,iV1,iV2], n_c1/n_c)
    f = interpolate.interp1d(xs, ys, kind='cubic')
    x = np.linspace(min(xs), max(xs), num=101, endpoint=True)
    y = f(xwant)
    return y
               
def datafit(x,y,deg,xval):
    coef = polyfit(x, y, deg)
    p = poly1d(coef)
    yval = p(xval)
    return yval
        
    
def datafit_curve(x,y,deg):
    coef = polyfit(x, y, deg)
    p = poly1d(coef)
    
    xmin = min(x)
    xmax = max(x)
    
    xp = linspace(xmin, xmax, 100)
    
    return xp, p(xp)

def del_existing_file(file):
    if os.path.isfile(file):
        os.remove(file)
        
        
##################################################################
# The following get different types of data from DQMC output files
# Get_sign
# Get_Etot
# Get_den_avg
# Get_den_orb
# Get_spincorre_r0_3orb
# Get_spincorre_nn
# Get_Saf_orb
# Get_Scdw_orb
# Get_Pair_orb
#===========================
# Then for tdm measurement
#===========================
# Get_beta_Gbeta2
# Get_Gtau_orb
# Get_chi_r0_orb
# Get_chi_q_orb
# Get_chi_q_orb_stackedPAM
# Get_Pd_orb
#==============================
# Finally some useful examples
#==============================
##########################################################
# Below for static (no Gtau) measurements
##########################################################
def Get_sign(data):
    '''
    Obtain avg sign from .out files
    '''
    sign = zeros(2)
    sign[0], sign[1] = Get_dataline(data, ' Avg sign : ', './data/tmp')
    return sign

def Get_Etot(data):
    '''
    Obtain total energy from .out files
    '''
    Etot = zeros(2)
    Etot[0], Etot[1] = Get_dataline(data, ' Total energy : ', './data/tmp')
    return Etot

def Get_den_avg(data):
    '''
    Obtain average density from .out files
    '''
    den = zeros(2)
    den[0], den[1] = Get_dataline(data, ' Density : ', './data/tmp')
    return den

def Get_den_orb(data, norb, Nline):
    '''
    Obtain average density for each orbital from .out files
    '''
    dens = zeros((norb,2))
    
    phrase = ' Mean Equal time Green'
    
    if os.path.isfile('./data/tmp'):
        os.remove('./data/tmp')
            
    f = open(data, 'r')
    find_lines(phrase, f, 1, Nline, './data/tmp') 
 
    a = loadtxt('./data/tmp', skiprows=0)   

    for i in range(len(a)):   
        if a[i,0]==a[i,1] and abs(a[i,3])<1.e-4 and abs(a[i,4])<1.e-4 and abs(a[i,5])<1.e-4:
            io = int(a[i,0])
            dens[io,0] = 2.*(1.0-a[i,6])
            dens[io,1] = a[i,7]

    return dens

def Get_spincorre_r0_3orb(data, mode, Ncorr, Nline):
    '''
    Obtain static spin XX and ZZ correlation between orbitals from .out files
    Only applies for 3 orbital case right now because Ncorr denotes orb=(0,1),(0,2),(1,2) in order !!!
    TODO: more general case
    '''
    Scorr = zeros((Ncorr,2))
    
    phraseX = ' XX Spin correlation function'
    phraseZ = ' ZZ Spin correlation function'
    
    if os.path.isfile('./data/tmpX'):
        os.remove('./data/tmpX')
    if os.path.isfile('./data/tmpZ'):
        os.remove('./data/tmpZ')
            
    f = open(data, 'r')
    find_lines(phraseX, f, 1, Nline, './data/tmpX') 
    f = open(data, 'r')
    find_lines(phraseZ, f, 1, Nline, './data/tmpZ') 
    
    a = loadtxt('./data/tmpX', skiprows=0)
    b = loadtxt('./data/tmpZ', skiprows=0)
    assert(len(a)==len(b))
    
    for i in range(0,len(a)):   
        # orb=(0,1) 
        if a[i,0]==0 and a[i,1]==1 and abs(a[i,3])<1.e-4 and abs(a[i,4])<1.e-4 and abs(a[i,5]-1.0)<1.e-4:
            xavg = a[i,6]; xerr = a[i,7]
            zavg = b[i,6]; zerr = b[i,7]
            
            if mode=='minerrbar':
                if xerr<zerr:
                    Scorr[0,0] = xavg; Scorr[0,1] = xerr
                else:
                    Scorr[0,0] = zavg; Scorr[0,1] = zerr
            elif mode=='averageXZ':
                Scorr[0,0] = (2.*xavg+zavg)/3.
                Scorr[0,1] = (2.*xerr+zerr)/3.

        # orb=(0,2) 
        if a[i,0]==0 and a[i,1]==2 and abs(a[i,3])<1.e-4 and abs(a[i,4])<1.e-4 and abs(a[i,5]-2.0)<1.e-4:
            xavg = a[i,6]; xerr = a[i,7]
            zavg = b[i,6]; zerr = b[i,7]
            
            if mode=='minerrbar':
                if xerr<zerr:
                    Scorr[1,0] = xavg; Scorr[1,1] = xerr
                else:
                    Scorr[1,0] = zavg; Scorr[1,1] = zerr
            elif mode=='averageXZ':
                Scorr[1,0] = (2.*xavg+zavg)/3.
                Scorr[1,1] = (2.*xerr+zerr)/3.

        # orb=(1,2) 
        if a[i,0]==1 and a[i,1]==2 and abs(a[i,3])<1.e-4 and abs(a[i,4])<1.e-4 and abs(a[i,5]-1.0)<1.e-4:
            xavg = a[i,6]; xerr = a[i,7]
            zavg = b[i,6]; zerr = b[i,7]
            
            if mode=='minerrbar':
                if xerr<zerr:
                    Scorr[2,0] = xavg; Scorr[2,1] = xerr
                else:
                    Scorr[2,0] = zavg; Scorr[2,1] = zerr
            elif mode=='averageXZ':
                Scorr[2,0] = (2.*xavg+zavg)/3.
                Scorr[2,1] = (2.*xerr+zerr)/3.
                
    return Scorr

def Get_spincorre_nn(data, mode, norb, Nline):
    '''
    Obtain static nearest neighbor spin XX and ZZ correlation for each orbital from .out files
    '''
    Scorr = zeros((norb,2))
    
    phraseX = ' XX Spin correlation function'
    phraseZ = ' ZZ Spin correlation function'
    
    if os.path.isfile('./data/tmpX'):
        os.remove('./data/tmpX')
    if os.path.isfile('./data/tmpZ'):
        os.remove('./data/tmpZ')
            
    f = open(data, 'r')
    find_lines(phraseX, f, 1, Nline, './data/tmpX') 
    f = open(data, 'r')
    find_lines(phraseZ, f, 1, Nline, './data/tmpZ') 
    
    a = loadtxt('./data/tmpX', skiprows=0)
    b = loadtxt('./data/tmpZ', skiprows=0)
    assert(len(a)==len(b))
    
    for i in range(0,len(a)):   
        if a[i,0]==a[i,1]:
            io = int(a[i,0])
            
            if abs(a[i,3]-1.)<1.e-4 and abs(a[i,4])<1.e-4 and abs(a[i,5])<1.e-4:
                xavg = a[i,6]; xerr = a[i,7]
                zavg = b[i,6]; zerr = b[i,7]
            
                if mode=='minerrbar':
                    if xerr<zerr:
                        Scorr[io,0] = xavg; Scorr[io,1] = xerr
                    else:
                        Scorr[io,0] = zavg; Scorr[io,1] = zerr
                elif mode=='averageXZ':
                    Scorr[io,0] = (2.*xavg+zavg)/3.
                    Scorr[io,1] = (2.*xerr+zerr)/3.
                
    return Scorr

def Get_Saf_orb(data, mode, norb):
    '''
    Obtain Saf(pi,pi) between same orbital from local_orb* files
    Note that local_orb* file orb index starts from 1 instead of 0
    '''
    headlines = 3*(norb+1)+2
    footlines = norb*norb+10
    a = genfromtxt(data, skip_header=headlines, skip_footer=footlines)
    Nlines = norb*norb
    assert(len(a)==Nlines)
    
    Saf = zeros((norb, 2))

    for i in range(0,Nlines):
        if a[i,0]==a[i,1]:
            io = int(a[i,0])-1  # note -1 because DQMC labels from 1 (not 0 !!!)

            if mode=='minerrbar':
                if a[i,3]<a[i,5]:
                    Saf[io,0] = a[i,2]; Saf[io,1] = a[i,3]
                else:
                    Saf[io,0] = a[i,4]; Saf[io,1] = a[i,5]
            elif mode=='averageXZ':
                Saf[io,0] = (2.*a[i,2]+a[i,4])/3.
                Saf[io,1] = (2.*a[i,3]+a[i,5])/3.
                                                   
    return Saf

def Get_Scdw_orb(data, norb):
    '''
    Obtain Scdw(pi,pi) between same orbital from local_orb* files
    Note that local_orb* file orb index starts from 1 instead of 0
    '''
    if os.path.isfile('./data/tmp'):
        os.remove('./data/tmp')
        
    f = open(data, 'r')
    phrase = 'nnprod '
    find_lines(phrase, f, 1, norb+1, './data/tmp')   # +1 accounts for total density Scdw
    a = loadtxt('./data/tmp', skiprows=0)
    assert(len(a)==norb+1)
    
    if os.path.isfile('./data/tmp'):
        os.remove('./data/tmp')
        
    f = open(data, 'r')
    phrase = 'nnsum '
    find_lines(phrase, f, 1, norb+1, './data/tmp')   # +1 accounts for total density Scdw
    b = loadtxt('./data/tmp', skiprows=0)
    assert(len(b)==norb+1)
    
    Scdw = zeros((norb+1, 2))
    
    # normally n+n's (pi,pi) modulation is neglible
    if (max(abs(b[:,1]))<1.e-4):
        # first Scdw for each orbital/layer
        for i in range(norb):
            Scdw[i,0] = a[i+1,1]; Scdw[i,1] = a[i+1,2]
            
        # total Scdw
        Scdw[norb,0] = a[0,1]; Scdw[norb,1] = a[0,2]
                                
    return Scdw

def Get_spin_orb(data, Nline, pairsites, direction):
    '''
    Obtain spin correlation between site/orbital from .out files
    Only get the value for specific orbitals and distance stored in pairsites
    '''
    Npair = len(pairsites)
    spincorre = zeros((Npair,2))
    
    if direction=='XX':
        phrase = ' XX Spin correlation function'
    elif direction=='ZZ':
        phrase = ' ZZ Spin correlation function'
    
    if os.path.isfile('./data/tmp'):
        os.remove('./data/tmp')
            
    f = open(data, 'r')
    find_lines(phrase, f, 1, Nline, './data/tmp') 

    for p in range(Npair):
        site1 = pairsites[p][0]
        site2 = pairsites[p][1]
        x = pairsites[p][2]
        y = pairsites[p][3]
        a = loadtxt('./data/tmp', skiprows=0)   
        
        for i in range(len(a)):   
            if a[i,0]==site1 and a[i,1]==site2 and abs(a[i,3]-x)<1.e-4 and abs(a[i,4]-y)<1.e-4:
                spincorre[p,0] = a[i,6]
                spincorre[p,1] = a[i,7]

    return spincorre

def Get_Pair_orb(data, Nline, pairsites):
    '''
    Obtain pair correlation between site/orbital from .out files
    Only get the value for specific orbitals and distance stored in pairsites
    '''
    Npair = len(pairsites)
    paircorre = zeros((Npair,2))
    
    phrase = ' Pairing correlation function'
    
    if os.path.isfile('./data/tmp'):
        os.remove('./data/tmp')
            
    f = open(data, 'r')
    find_lines(phrase, f, 1, Nline, './data/tmp') 

    for p in range(Npair):
        site1 = pairsites[p][0]
        site2 = pairsites[p][1]
        x = pairsites[p][2]
        y = pairsites[p][3]
        a = loadtxt('./data/tmp', skiprows=0)   
        
        for i in range(len(a)):   
            if a[i,0]==site1 and a[i,1]==site2 and abs(a[i,3]-x)<1.e-4 and abs(a[i,4]-y)<1.e-4:
                paircorre[p,0] = a[i,6]
                paircorre[p,1] = a[i,7]

    return paircorre

def Get_Pair_corre(data, Nline):
    '''
    Obtain all pair correlation between site/orbital from .out files
    '''
    phrase = ' Pairing correlation function'
    
    if os.path.isfile('./data/tmp'):
        os.remove('./data/tmp')
            
    f = open(data, 'r')
    find_lines(phrase, f, 1, Nline, './data/tmp') 
    a = loadtxt('./data/tmp', skiprows=0)   
        
    return a

##########################################################
# Below for tdm measurements
##########################################################
def Get_beta_Gbeta2(data, norb, beta, Ntau):
    '''
    Obtain beta*G(beta/2) from Gr0* files
    '''
    G = zeros((norb,2))
    
    for i in range(0,norb):   
        if os.path.isfile('./data/tmp'):
            os.remove('./data/tmp')
        
        phrase = ' Gfun      '+str(i)+'   '+str(i)+' '
        f = open(data, 'r')
        find_lines(phrase, f, 2, Ntau, './data/tmp') 
        
        a = loadtxt('./data/tmp', skiprows=0)
        assert(len(a)==Ntau and abs(a[0,0])<1.e-5)
        assert(abs(a[Ntau/2,0]-beta/2.)<1.e-5)
        
        G[i,0] = a[Ntau/2,1]
        G[i,1] = a[Ntau/2,2]

    return G

def Get_Gtau_orb(data, norb, beta, Ntau, Gtaudata_name):
    '''
    Obtain local G(tau) from Gr0* files for maxent usage
    and store them in Gtaudata_name file
    '''
    G = zeros((Ntau, norb))
    
    for i in range(0,norb):   
        dataname = Gtaudata_name+'_orb'+str(i)
        if os.path.isfile(dataname):
            os.remove(dataname)
        
        f = open(dataname,'a') 
        f.write(str(Ntau)+'\n')
    
        phrase = ' Gfun      '+str(i)+'   '+str(i)+' '
        f = open(data, 'r')
        find_lines(phrase, f, 2, Ntau, dataname) 
        
        # replace zero errobar to artificial small number
        file = open(dataname, "r")
        text = file.read()
        file.close()
    
        text = text.replace("0.00000000E+00" , "0.10000000E-04")
        file = open(dataname, "w")
        file.write(text)
        file.close()

def Get_chi_r0_orb(data, mode, norb):
    '''
    Obtain chi(r=0) between same orbital from chi_r0_iw0* files
    chi(norb,norb,2): avg and err
    Note: need to divide by 2 between two different orb since DQMC code treats (i,j)=(j,i)
    '''
    a = loadtxt(data, skiprows=1)
    Nlines = norb*(norb+1)/2
    assert(len(a)==Nlines)
    
    chi_r0 = zeros((norb,norb,2))

    for i in range(0,Nlines):
        o1 = int(a[i,0])
        o2 = int(a[i,1])
            
        if o1==o2:
            if mode=='minerrbar':
                if a[i,3]<a[i,5]:
                    chi_r0[o1,o2,0] = a[i,2]; chi_r0[o1,o2,1] = a[i,3]
                else:
                    chi_r0[o1,o2,0] = a[i,4]; chi_r0[o1,o2,1] = a[i,5]
            elif mode=='averageXZ':
                chi_r0[o1,o2,0] = (2.*a[i,2]+a[i,4])/3.
                chi_r0[o1,o2,1] = (2.*a[i,3]+a[i,5])/3.
        else:
            if mode=='minerrbar':
                if a[i,3]<a[i,5]:
                    chi_r0[o1,o2,0] = a[i,2]/2.; chi_r0[o1,o2,1] = a[i,3]/2.
                    chi_r0[o2,o1,0] = a[i,2]/2.; chi_r0[o2,o1,1] = a[i,3]/2.
                else:
                    chi_r0[o1,o2,0] = a[i,4]/2.; chi_r0[o1,o2,1] = a[i,5]/2.
                    chi_r0[o2,o1,0] = a[i,4]/2.; chi_r0[o2,o1,1] = a[i,5]/2.
            elif mode=='averageXZ':
                chi_r0[o1,o2,0] = (2.*a[i,2]+a[i,4])/3./2.
                chi_r0[o1,o2,1] = (2.*a[i,3]+a[i,5])/3./2.
                chi_r0[o2,o1,0] = (2.*a[i,2]+a[i,4])/3./2.
                chi_r0[o2,o1,1] = (2.*a[i,3]+a[i,5])/3./2.
                              
    return chi_r0

def Get_chi_q_orb(data, mode, norb):
    '''
    Obtain chi(q) between same/different orbital at q=(0,0) and (pi,pi) from chi_q_iw0* files
    chi(norb, norb, 2): 2 stores avg and err
    '''
    a = genfromtxt(data, skip_header=2, skip_footer=2)
    norb2 = norb*norb
    Nlines = 2*(norb2+1)
    assert(len(a)==Nlines)
    
    chi_q0  = zeros((norb, norb, 2))
    chi_qpi = zeros((norb, norb, 2))
    chi_q0_total  = zeros(2)
    chi_qpi_total = zeros(2)

    # q=0
    for i in range(0,Nlines/2):
        if a[i,0]!=100:
            o1 = int(a[i,0])
            o2 = int(a[i,1])

            if mode=='minerrbar':
                if a[i,3]<a[i,5]:
                    chi_q0[o1,o2,0] = a[i,2]; chi_q0[o1,o2,1] = a[i,3]
                else:
                    chi_q0[o1,o2,0] = a[i,4]; chi_q0[o1,o2,1] = a[i,5]
            elif mode=='averageXZ':
                chi_q0[o1,o2,0] = (2.*a[i,2]+a[i,4])/3.
                chi_q0[o1,o2,1] = (2.*a[i,3]+a[i,5])/3.
                
        elif a[i,0]==100:
            if mode=='minerrbar':
                if a[i,3]<a[i,5]:
                    chi_q0_total[0] = a[i,2]; chi_q0_total[1] = a[i,3]
                else:
                    chi_q0_total[0] = a[i,4]; chi_q0_total[1] = a[i,5]
            elif mode=='averageXZ':
                chi_q0_total[0] = (2.*a[i,2]+a[i,4])/3.
                chi_q0_total[1] = (2.*a[i,3]+a[i,5])/3.
            
    # q=pi
    for i in range(Nlines/2,Nlines):
        if a[i,0]!=100:
            o1 = int(a[i,0])
            o2 = int(a[i,1])

            if mode=='minerrbar':
                if a[i,3]<a[i,5]:
                    chi_qpi[o1,o2,0] = a[i,2]; chi_qpi[o1,o2,1] = a[i,3]
                else:
                    chi_qpi[o1,o2,0] = a[i,4]; chi_qpi[o1,o2,1] = a[i,5]
            elif mode=='averageXZ':
                chi_qpi[o1,o2,0] = (2.*a[i,2]+a[i,4])/3.
                chi_qpi[o1,o2,1] = (2.*a[i,3]+a[i,5])/3.
                
        elif a[i,0]==100:
            if mode=='minerrbar':
                if a[i,3]<a[i,5]:
                    chi_qpi_total[0] = a[i,2]; chi_qpi_total[1] = a[i,3]
                else:
                    chi_qpi_total[0] = a[i,4]; chi_qpi_total[1] = a[i,5]
            elif mode=='averageXZ':
                chi_qpi_total[0] = (2.*a[i,2]+a[i,4])/3.
                chi_qpi_total[1] = (2.*a[i,3]+a[i,5])/3.
                                                   
    return chi_q0, chi_qpi, chi_q0_total, chi_qpi_total

def Get_chi_q_orb_stackedPAM(data, mode, norb):
    '''
    Obtain chi(q) between same/different orbital at q=(0,0) and (pi,pi) from chi_q_iw0* files
    chi(norb, norb, 2): 2 stores avg and err
    '''
    a = genfromtxt(data, skip_header=2, skip_footer=2)
    norb2 = norb*norb
    Nlines = 2*(norb2+3)  # there are 2 more components for individual chi for f1 and f2
    assert(len(a)==Nlines)
    
    chi_q0  = zeros((norb, norb, 2))
    chi_qpi = zeros((norb, norb, 2))
    chi_q0_total  = zeros((3,2))
    chi_qpi_total = zeros((3,2))

    # q=0
    for i in range(0,Nlines/2-3):
        if a[i,0]!=100:
            o1 = int(a[i,0])
            o2 = int(a[i,1])

            if mode=='minerrbar':
                if a[i,3]<a[i,5]:
                    chi_q0[o1,o2,0] = a[i,2]; chi_q0[o1,o2,1] = a[i,3]
                else:
                    chi_q0[o1,o2,0] = a[i,4]; chi_q0[o1,o2,1] = a[i,5]
            elif mode=='averageXZ':
                chi_q0[o1,o2,0] = (2.*a[i,2]+a[i,4])/3.
                chi_q0[o1,o2,1] = (2.*a[i,3]+a[i,5])/3.
                
    if mode=='minerrbar':
        i = Nlines/2-3
        if a[i,3]<a[i,5]:
            chi_q0_total[0,0] = a[i,2]; chi_q0_total[0,1] = a[i,3]
        else:
            chi_q0_total[0,0] = a[i,4]; chi_q0_total[0,1] = a[i,5]
            
        i = Nlines/2-2
        if a[i,3]<a[i,5]:
            chi_q0_total[1,0] = a[i,2]; chi_q0_total[1,1] = a[i,3]
        else:
            chi_q0_total[1,0] = a[i,4]; chi_q0_total[1,1] = a[i,5]
            
        i = Nlines/2-1
        if a[i,3]<a[i,5]:
            chi_q0_total[2,0] = a[i,2]; chi_q0_total[2,1] = a[i,3]
        else:
            chi_q0_total[2,0] = a[i,4]; chi_q0_total[2,1] = a[i,5]
            
    elif mode=='averageXZ':
        i = Nlines/2-3
        chi_q0_total[0,0] = (2.*a[i,2]+a[i,4])/3.
        chi_q0_total[0,1] = (2.*a[i,3]+a[i,5])/3.
        i = Nlines/2-2
        chi_q0_total[1,0] = (2.*a[i,2]+a[i,4])/3.
        chi_q0_total[1,1] = (2.*a[i,3]+a[i,5])/3.
        i = Nlines/2-1
        chi_q0_total[2,0] = (2.*a[i,2]+a[i,4])/3.
        chi_q0_total[2,1] = (2.*a[i,3]+a[i,5])/3.
            
    # q=pi
    for i in range(Nlines/2,Nlines-3):
        if a[i,0]!=100:
            o1 = int(a[i,0])
            o2 = int(a[i,1])

            if mode=='minerrbar':
                if a[i,3]<a[i,5]:
                    chi_qpi[o1,o2,0] = a[i,2]; chi_qpi[o1,o2,1] = a[i,3]
                else:
                    chi_qpi[o1,o2,0] = a[i,4]; chi_qpi[o1,o2,1] = a[i,5]
            elif mode=='averageXZ':
                chi_qpi[o1,o2,0] = (2.*a[i,2]+a[i,4])/3.
                chi_qpi[o1,o2,1] = (2.*a[i,3]+a[i,5])/3.
                
    if mode=='minerrbar':
        i = Nlines-3
        if a[i,3]<a[i,5]:
            chi_qpi_total[0,0] = a[i,2]; chi_qpi_total[0,1] = a[i,3]
        else:
            chi_qpi_total[0,0] = a[i,4]; chi_qpi_total[0,1] = a[i,5]
            
        i = Nlines-2
        if a[i,3]<a[i,5]:
            chi_qpi_total[1,0] = a[i,2]; chi_qpi_total[1,1] = a[i,3]
        else:
            chi_qpi_total[1,0] = a[i,4]; chi_qpi_total[1,1] = a[i,5]
            
        i = Nlines-1
        if a[i,3]<a[i,5]:
            chi_qpi_total[2,0] = a[i,2]; chi_qpi_total[2,1] = a[i,3]
        else:
            chi_qpi_total[2,0] = a[i,4]; chi_qpi_total[2,1] = a[i,5]
            
    elif mode=='averageXZ':
        i = Nlines-3
        chi_qpi_total[0,0] = (2.*a[i,2]+a[i,4])/3.
        chi_qpi_total[0,1] = (2.*a[i,3]+a[i,5])/3.
        i = Nlines-2
        chi_qpi_total[1,0] = (2.*a[i,2]+a[i,4])/3.
        chi_qpi_total[1,1] = (2.*a[i,3]+a[i,5])/3.
        i = Nlines-1
        chi_qpi_total[2,0] = (2.*a[i,2]+a[i,4])/3.
        chi_qpi_total[2,1] = (2.*a[i,3]+a[i,5])/3.
                                                   
    return chi_q0, chi_qpi, chi_q0_total, chi_qpi_total


def Get_Ps_orb(data, norb):
    '''
    Obtain Ps, Gammas*Ps0 for each orbital from .tdm.out files
    Note that some orb, for example, U=0 orb, does not have Ps value !!
    '''
    Ps  = zeros((norb, 2))
    Ps0 = zeros((norb, 2))
    aGsPs0 = zeros((norb, 2))
    cGsPs0 = 0.0
    
    for i in range(norb):
        Ps[i,0], Ps[i,1], Ps0[i,0], Ps0[i,1] = Get_dataline(data, '     Ps and Ps0 = ', './data/tmp')
        aGsPs0[i,0], aGsPs0[i,1] = Get_dataline(data, 'accumulated Gs*Ps0 = ', './data/tmp')
    cGsPs0 = Get_dataline(data, 'calculated Gd*Ps0 = ', './data/tmp')

    return Ps, Ps0, aGsPs0, cGsPs0


def Get_Pd_orb(data, norb):
    '''
    Obtain Pd, Gammad*Pd0 for each orbital from .tdm.out files
    Note that some orb, for example, U=0 orb, does not have Pd value !!
    '''
    Pd  = zeros((norb, 2))
    Pd0 = zeros((norb, 2))
    
    for i in range(1,norb):
        Pd[i,0], Pd[i,1], Pd0[i,0], Pd0[i,1] = Get_dataline(data, ' Pd and Pd0 =   '+str(i), './data/tmp')

    return Pd, Pd0


def Get_Sigma_k(data, o1, o2, Kind, Ntau):
    '''
    Obtain intra- or inter-orbital (labeled by o1,o2) self-energy 
    at particular k (indexed by Kind) from .tdm.out files
    Ntau = Nw (No. of Matsubara frequency)
    '''
    ws = zeros(Ntau)
    Sigma = zeros((Ntau,2))    # 2 denotes Re and Im parts
    
    if Kind<10:
        phrase = 'Gfun SelfEn k=  '+str(Kind)+'   cell_site_pair=  '+str(o1)+',  '+str(o2)
    else:
        phrase = 'Gfun SelfEn k= ' +str(Kind)+'   cell_site_pair=  '+str(o1)+',  '+str(o2)
    
    if os.path.isfile('./data/tmp'):
        os.remove('./data/tmp')
            
    f = open(data, 'r')
    
    find_lines(phrase, f, 1, Ntau, './data/tmp') 
 
    a = loadtxt('./data/tmp', skiprows=0)   
    
    ws = a[:,0]
    Sigma[:,0] = a[:,1]
    Sigma[:,1] = a[:,3]

    return ws, Sigma


##########################################################
# Below for some useful examples
##########################################################
'''
tpd_nn_hop_fac = {('dx2y2','L','px'):   tpd,\
                  ('dx2y2','R','px'):  -tpd,\
                  ('dx2y2','U','py'):   tpd,\
                  ('dx2y2','D','py'):  -tpd,\
                  # below just inverse dir of the above one by one
                  ('px','R','dx2y2'):   tpd,\
                  ('px','L','dx2y2'):  -tpd,\
                  ('py','D','dx2y2'):   tpd,\
                  ('py','U','dx2y2'):  -tpd}
                  
plot(tps,S01,'--', marker=Ms[iV], color=colors[iV], markersize=8, fillstyle='left', markerfacecolor='None'
'''

