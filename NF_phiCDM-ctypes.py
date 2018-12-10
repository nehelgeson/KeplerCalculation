#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 18:35:10 2018

@author: narayankhadka
"""
from scipy.integrate import quad
import scipy.integrate as integrate
from scipy.integrate import odeint #ODE solver
#import matplotlib.pyplot as plt
import numpy
import math
import time
import scipy
from scipy import LowLevelCallable
import smtplib
import ctypes
from numpy import savetxt, loadtxt, arange, trapz
from numpy.linalg import inv
start_time = time.time()

_integrate = ctypes.CDLL('./kepler.so')
_integrate.SetGlobals.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_double, ctypes.c_double)
_integrate.IntegrateFunc.argtypes = (ctypes.c_double, ctypes.c_void_p)
_integrate.IntegrateFunc.restype = ctypes.c_double


for case in range(1, 3, 1):
    DM_obs = numpy.array([1518, 1977, 2283])
    H_obs = numpy.array([81.5, 90.4, 97.3])
    
    cHobs = 9.0 #Font-Ribera
    dobs = 0.336 #Beutler (Farooq)
    DA_obs = 10.8 #Font-Ribera
    B_obs = 13.94 #Bautista
    
    rfid = 147.60 #Planck
    
    z_obs, Hz_obs, sigHobs = loadtxt('H(z)data.dat',unpack = True)
    
    z = numpy.array([0.38, 0.51, 0.61, 0.106, 0.15, 1.52, 2.33, 2.36])
    
    s = numpy.array([22, 1.9, 27, 1.9, 32, 2.1])
    
    Cov = (1/10000)*(numpy.matrix([[10000,2280,4970,1117,1991,520],
                                   [2280,10000,1536,4873,984,2307],
                                   [4970,1536,10000,2326,5120,1211],
                                   [1117,4873,2326,10000,1571,5449],
                                   [1991,984,5120,1571,10000,2408],
                                   [520,2307,1211,5449,2408,10000]]))
    p1 = numpy.outer(s, s)
    p2 = numpy.multiply(p1, Cov)
    Cinv = inv(p2)
    
    O1 = []
    a1 = []
    Ok1 = []
    Like = []
    
    if case == 1:
        s = 2.8 
        H0_av = 68 
    if case == 2:
        s = 1.74 
        H0_av = 73.24 
    
    
    
    c = 299792458/1000
    Th = 2.7255/2.7 
    gamma = (H0_av**2)/(s**2) + sum((Hz_obs**2)/(sigHobs**2))
    pi = math.pi
    m = 1
    n = 5
    
    def phiCDM(w, t, zz):
        
        p, v, a = w 
        al, k, m, K = zz
    
        f = [v,
             -3*v*(((4/(9*a**3))) + (1/12)*(v**2 + (k*m)/(p**(al))) - K/(a**2))**(1/2)
             + ((k*al*m)/2)/(p**(al + 1)), 
            (((4/9)/a) + (((a**2)/12))*(v**2 + (k*m)*(p**(-al))) - K)**(1/2)]
        return f
    
    def rs(H0, O):
        h = H0/100
        Ob = 0.02221/(h**2)
        b2 = 0.238*((O*(h**2))**(0.223))
        b1 = 0.313*((O*(h**2))**(-0.419))*(1 + (0.607*(O*(h**2))**(0.674)))
        zd = 1291*(((O*(h**2))**(0.251))/(1 
                   + (0.659*(O*(h**2))**(0.828))))*(1 + ((b1*(Ob*(h**2)))**(b2)))
        Rd = 31.5*(Ob*(h**2))*(Th**(-4))*(1000/zd)
        zeq = 25000*(O*(h**2))*(Th**(-4))
        Req = 31.5*(Ob*(h**2))*(Th**(-4))*(1000/zeq)
        keq = 0.0746*(O*(h**2))*(Th**(-2))
        A = (1 + Rd)**(1/2)
        B = (Rd + Req)**(1/2)
        C = 1 + (Req)**(1/2)
        return (2/(3*keq))*((6/Req)**(1/2))*(numpy.log((A + B)/C))
    
    def E(O, red, Ok_0, Ophiz):
        g1 = (O*((1 + red)**3) + Ok_0*((1 + red)**2) + Ophiz)**(1/2)
        return g1
    
    def D_M(H0, Ok_0, x):
        dH = c/H0
        if Ok_0 < 0:            
            y = (1/((-(Ok_0))**(1/2)))*(math.sin(((-(Ok_0))**(1/2))*h0*afin*x))
        if Ok_0 > 0:            
            y = (1/(((Ok_0))**(1/2)))*(math.sinh((((Ok_0))**(1/2))*h0*afin*x))
        if Ok_0 == 0:
            y = h0*afin*x
        return y*dH
    
    def chi_sq(H0, O, Ok_0):
        dH = c/H0
        r = rs(H0, O)/rfid
        DM_th = []
        H_th = []
        for q in range( 0, 8, 1):
            z1 = z[q]
            H1 = H0*E(O, z1, Ok_0, O_phi_z[q])
            y = D_M(H0, Ok_0, rr[q])/dH
            DM = D_M(H0, Ok_0, rr[q]) 
            if q <= 2:
                DM_th.append(DM)
                H_th.append(H1)
            if z1 == 0.15:
                DV = dH*(((y**2)*z1)/(E(O, z1, Ok_0, O_phi_z[q])))**(1/3)
                DVobs1 = 664*r
                unc1 = 25*r
                chi2DV1 = ((DV - DVobs1)**2)/(unc1**2)
            if z1 == 0.106:
                DV = (c/H0)*(((y**2)*(z1))/(E(O, z1, Ok_0, O_phi_z[q])))**(1/3)
                dth = (rs(H0, O))/DV
                chi2d = ((dth - dobs)**2)/(0.015**2)
            if z1 == 1.52:
                DV = dH*(((y**2)*z1)/(E(O, z1, Ok_0, O_phi_z[q])))**(1/3)
                DVobs2 = 3855*r #Ata
                unc2 = 170*r
                chi2DV2 = ((DV - DVobs2)**2)/(unc2**2)
            if z1 == 2.33:
                DH = c/H1
                F = DH**(0.7)
                G = DM**(0.3)
                B_th = (F*G)/(rs(H0, O))
                chi2B = ((B_th - B_obs)**2)/(0.35**2)
            if z1 == 2.36:
                
                cHth = (c/(rs(H0, O)))*(1/H1)
                chi2H = ((cHth - cHobs)**2)/(0.3**2)
        Delta = numpy.array([(DM_th[0]/r - DM_obs[0]), 
                             (r*H_th[0] - H_obs[0]),
                             (DM_th[1]/r - DM_obs[1]),
                             (r*H_th[1] - H_obs[1]),
                             (DM_th[2]/r - DM_obs[2]),
                             (r*H_th[2] - H_obs[2])]) 
        prod = (Cinv.dot(Delta)).T
        chi_sq1 = Delta.dot(prod)
        chi_sq_11 = chi_sq1.item((0,0))
        return chi2B + chi2H + chi2d + chi_sq_11 + chi2DV1 + chi2DV2

    def SetGlobals(): #Insert all the globals we need
        Cinv = inv(p2)
        z_type = ctypes.c_double * len(z)
        O_phi_z_type = ctypes.c_double * len(O_phi_z)
        rr_type = ctypes.c_double * len(rr)
        Cinv_type = ctypes.c_double * (len(Cinv) * len(Cinv))
        DM_obs_type = ctypes.c_double * len(DM_obs)
        H_obs_type = ctypes.c_double * len(H_obs)

        cinv_l = list()
        for row in Cinv.tolist():
            cinv_l.extend(row)
        
        _integrate.SetGlobals(ctypes.c_double(O),
                              ctypes.c_double(Ok_0),
                              ctypes.c_double(c),
                              ctypes.c_double(Th),
                              ctypes.c_double(rfid),
                              z_type(*z),
                              O_phi_z_type(*O_phi_z),
                              rr_type(*rr),
                              ctypes.c_double(h0),
                              ctypes.c_double(afin),
                              ctypes.c_double(dobs),
                              ctypes.c_double(B_obs),
                              ctypes.c_double(cHobs),
                              Cinv_type(*cinv_l),
                              DM_obs_type(*DM_obs),
                              H_obs_type(*H_obs),
                              ctypes.c_double(H0_av),
                              ctypes.c_double(s))


    
    def Li(O, Ok_0):
        SetGlobals()

        test_result = _integrate.IntegrateFunc(ctypes.c_double(2.34), ctypes.c_void_p())
        print("Test result: %f" % (test_result))
        
        new_f = LowLevelCallable(_integrate.IntegrateFunc)
        l, error = integrate.quad(new_f, 
            H0_av - n*s, H0_av + n*s, epsabs=0, epsrel=1.49e-08) #Set error bounds?
        return l/(math.sqrt(2*pi*(s**2)))
        
    def Li_2(alpha0, beta0):
        alpha = 1/(s**2) + alpha0
        beta = H0_av/(s**2) + beta0
        V = (1/2)*((alpha)*(s**2))**(-1/2)
        W = gamma - ((beta)**2)/alpha
        X = beta/(math.sqrt(2*alpha))
        Y = 1 + scipy.special.erf(X)
        Z = math.exp((-1/2)*W)
        return V*Z*Y

    
    def Ofunc(d,K):
        Omegam1 = (4./9.)*(1./(sol[d,2])**3.)
        Omegaphi1 = (1./12.)*(((sol[d,1])**2.) + k/((sol[d,0])**al))
        Omegak1 = -K/((sol[d,2])**2.)
        Omegam = Omegam1/(Omegam1 + Omegaphi1 + Omegak1)
        return Omegam, Omegam1, Omegaphi1, Omegak1
    for K in arange(-0.45,-0.4,0.01):
        print(K)
        for al in arange(0.4, 0.5, 0.01):
           print(al)
           k = (8/3)*((al + 4)/(al + 2))*(((2/3)*(al*(al + 2))))**(al/2) #This is kappa,
        #from eq. (2) of arXiv:1307.7399v1.
        
           t0 = 10**(-4)
           p0 = (((2/3)*(al*(al + 2)))**(1/2))*(t0)**(2/(al + 2)) #Initial value of phi.
           v0 = (((8/3)*al*(1/(al + 2)))**(1/2))/(t0)**((al)/(al + 2)) #Initial value of d(phi)/dt.
           a0 = t0**(2/3) #I assumed a ~ t^(2/3) in the early universe.
        
           t = numpy.arange(0, 150, t0)
        
        #arrays of initial conditions and parameters
           w0 = [p0, v0, a0]
           zz = [al, k, m, K]
        
        #solution array
           sol = odeint(phiCDM, w0, t, args=(zz,))
        #for H0 in arange(50, 85.1, 0.1):
           for O in arange(0.2, 0.4, 0.01):
            #start_time = time.time()
            #print(O)
               for b in range (0, int(150/t0), 1):
                   if (O >= Ofunc(b,K)[0]):
                       break
                
               h0 = (sol[b+1, 2] - sol[b, 2])/(sol[b,2])
               afin = sol[b,2]
               Omegam1, Omegaphi1, Omegak1 = Ofunc(b,K)[1], Ofunc(b,K)[2], Ofunc(b,K)[3]
               Ok_0 = Omegak1/(Omegam1 + Omegaphi1 + Omegak1)
            
               rr = []
               O_phi_z = []
               for q in range(0, 8, 1):
                   z1 = z[q]
                   r = 0
                   qq = 0
                   qt = 0
                   r1 = []
                   t1 = []
                   for d in range (0, b+1, 1):
                       if (sol[d,2]/afin >= 1/(1 + z1)):
                           r += 1/sol[d, 2]
                           qq += 1
                           qt += 1
                           r1.append(1/sol[d, 2])
                           t1.append(d)
                       if qq == 1:
                          Omegaphi_z = (1/12)*(((sol[d,1])**2) + k/((sol[d,0])**al))
                   if qt == 1:
                          rr.append(r)
                   else:
                        rr.append(trapz(r1, t1))
                   O_phi_z.append(Omegaphi_z/(Omegam1 + Omegak1 + Omegaphi1))
            
               alpha0 = 0
               beta0 = 0
            
               for q in range(0, len(z_obs), 1):
                   for g in range (0, int(150/t0), 1):
                       if (sol[g,2]/afin >= 1/(1 + z_obs[q])):
                           Omegaphi2 = (1/12)*(((sol[g,1])**2) + k/((sol[g,0])**al))
                           O_phi_z1 = Omegaphi2/(Omegam1 + Omegak1 + Omegaphi1)
                           break
                   alpha0 += (E(O, z_obs[q], Ok_0, O_phi_z1)**2)/((sigHobs[q])**2)
                   beta0 += (Hz_obs[q])*(E(O, z_obs[q], Ok_0, O_phi_z1))/((sigHobs[q])**2)
               O1.append(O)
               a1.append(al)
               Ok1.append(Ok_0)
               Like.append(Li(O, Ok_0)*Li_2(alpha0, beta0))
    print(Ok1)           
    print(O1)
    print(a1)
    print(Like)
    print(O_phi_z)
    if case == 1:
        savetxt('NF_phiCDM_Omega_m0_68_1.dat', O1)
        savetxt('NF_phiCDM_alpha_68_1.dat', a1)
        savetxt('NF_phiCDM_Like_68_1.dat', Like)
        savetxt('NF_phiCDM_Ok_68_1.dat', Ok1)
    if case == 2:
        savetxt('Nf_phiCDM_Omega_m0_73_1.dat', O1)
        savetxt('NF_phiCDM_alpha_73_1.dat', a1)
        savetxt('NF_phiCDM_Like_73_1.dat', Like)
        savetxt('NF_phiCDM_Ok_73_1.dat', Ok1)
print("--- %s seconds ---" % (time.time() - start_time))
