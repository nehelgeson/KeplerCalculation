#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 18:35:10 2018

@author: narayankhadka
"""

from scipy.integrate import quad
import scipy.integrate as integrate
from scipy.integrate import odeint  # ODE solver
# import matplotlib.pyplot as plt
import numpy
import math
import time
import scipy
import smtplib
from numpy import savetxt, loadtxt, arange, trapz
from numpy.linalg import inv
import multiprocessing

start_time = time.time()
for case in range(1, 3, 1):
    DM_obs = numpy.array([1518, 1977, 2283])
    H_obs = numpy.array([81.5, 90.4, 97.3])

    cHobs = 9.0  # Font-Ribera
    dobs = 0.336  # Beutler (Farooq)
    DA_obs = 10.8  # Font-Ribera
    B_obs = 13.94  # Bautista

    rfid = 147.60  # Planck

    z_obs, Hz_obs, sigHobs = loadtxt('H(z)data.dat', unpack=True)

    z = numpy.array([0.38, 0.51, 0.61, 0.106, 0.15, 1.52, 2.33, 2.36])

    s = numpy.array([22, 1.9, 27, 1.9, 32, 2.1])

    Cov = (1 / 10000) * (numpy.matrix([[10000, 2280, 4970, 1117, 1991, 520],
                                       [2280, 10000, 1536, 4873, 984, 2307],
                                       [4970, 1536, 10000, 2326, 5120, 1211],
                                       [1117, 4873, 2326, 10000, 1571, 5449],
                                       [1991, 984, 5120, 1571, 10000, 2408],
                                       [520, 2307, 1211, 5449, 2408, 10000]]))
    p1 = numpy.outer(s, s)
    p2 = numpy.multiply(p1, Cov)
    Cinv = inv(p2)

    O1 = []
    a1 = []
    Ok1 = []
    Like = []
    O_phi_z = []


    if case == 1:
        s = 2.8
        H0_av = 68
    if case == 2:
        s = 1.74
        H0_av = 73.24

    c = 299792458 / 1000
    Th = 2.7255 / 2.7
    gamma = (H0_av ** 2) / (s ** 2) + sum((Hz_obs ** 2) / (sigHobs ** 2))
    pi = math.pi
    m = 1
    n = 5


    def phiCDM(w, t, zz):

        p, v, a = w
        al, k, m, K = zz

        f = [v,
             -3 * v * ((4 / (9 * a ** 3)) + (1 / 12) * (v ** 2 + (k * m) / (p ** (al))) - K / (a ** 2)) ** (1 / 2)
             + ((k * al * m) / 2) / (p ** (al + 1)),
             (((4 / 9) / a) + ((a ** 2) / 12) * (v ** 2 + (k * m) * (p ** (-al))) - K) ** (1 / 2)]
        return f


    def rs(H0, O):
        h = H0 / 100
        Ob = 0.02221 / (h ** 2)
        b2 = 0.238 * ((O * (h ** 2)) ** 0.223)
        b1 = 0.313 * ((O * (h ** 2)) ** (-0.419)) * (1 + (0.607 * (O * (h ** 2)) ** 0.674))
        zd = 1291 * (((O * (h ** 2)) ** 0.251) / (1
                                                  + (0.659 * (O * (h ** 2)) ** 0.828))) * (
                         1 + ((b1 * (Ob * (h ** 2))) ** b2))
        Rd = 31.5 * (Ob * (h ** 2)) * (Th ** (-4)) * (1000 / zd)
        zeq = 25000 * (O * (h ** 2)) * (Th ** (-4))
        Req = 31.5 * (Ob * (h ** 2)) * (Th ** (-4)) * (1000 / zeq)
        keq = 0.0746 * (O * (h ** 2)) * (Th ** (-2))
        A = (1 + Rd) ** (1 / 2)
        B = (Rd + Req) ** (1 / 2)
        C = 1 + Req ** (1 / 2)
        return (2 / (3 * keq)) * ((6 / Req) ** (1 / 2)) * (numpy.log((A + B) / C))


    def E(O, red, Ok_0, Ophiz):
        g1 = (O * ((1 + red) ** 3) + Ok_0 * ((1 + red) ** 2) + Ophiz) ** (1 / 2)
        return g1


    def D_M(H0, Ok_0, x, h0, afin):
        dH = c / H0
        if Ok_0 < 0:
            y = (1 / ((-Ok_0) ** (1 / 2))) * (math.sin(((-Ok_0) ** (1 / 2)) * h0 * afin * x))
        if Ok_0 > 0:
            y = (1 / ((Ok_0) ** (1 / 2))) * (math.sinh((Ok_0 ** (1 / 2)) * h0 * afin * x))
        if Ok_0 == 0:
            y = h0 * afin * x
        return y * dH


    def chi_sq(H0, O, Ok_0, rr, O_phi_z, h0, afin):
        dH = c / H0
        r = rs(H0, O) / rfid
        DM_th = []
        H_th = []
        for q in range(0, 8, 1):
            z1 = z[q]
            H1 = H0 * E(O, z1, Ok_0, O_phi_z[q])
            y = D_M(H0, Ok_0, rr[q], h0, afin) / dH
            DM = D_M(H0, Ok_0, rr[q], h0, afin)
            if q <= 2:
                DM_th.append(DM)
                H_th.append(H1)
            if z1 == 0.15:
                DV = dH * (((y ** 2) * z1) / (E(O, z1, Ok_0, O_phi_z[q]))) ** (1 / 3)
                DVobs1 = 664 * r
                unc1 = 25 * r
                chi2DV1 = ((DV - DVobs1) ** 2) / (unc1 ** 2)
            if z1 == 0.106:
                DV = (c / H0) * (((y ** 2) * (z1)) / (E(O, z1, Ok_0, O_phi_z[q]))) ** (1 / 3)
                dth = (rs(H0, O)) / DV
                chi2d = ((dth - dobs) ** 2) / (0.015 ** 2)
            if z1 == 1.52:
                DV = dH * (((y ** 2) * z1) / (E(O, z1, Ok_0, O_phi_z[q]))) ** (1 / 3)
                DVobs2 = 3855 * r  # Ata
                unc2 = 170 * r
                chi2DV2 = ((DV - DVobs2) ** 2) / (unc2 ** 2)
            if z1 == 2.33:
                DH = c / H1
                F = DH ** 0.7
                G = DM ** 0.3
                B_th = (F * G) / (rs(H0, O))
                chi2B = ((B_th - B_obs) ** 2) / (0.35 ** 2)
            if z1 == 2.36:
                cHth = (c / (rs(H0, O))) * (1 / H1)
                chi2H = ((cHth - cHobs) ** 2) / (0.3 ** 2)
        Delta = numpy.array([(DM_th[0] / r - DM_obs[0]),
                             (r * H_th[0] - H_obs[0]),
                             (DM_th[1] / r - DM_obs[1]),
                             (r * H_th[1] - H_obs[1]),
                             (DM_th[2] / r - DM_obs[2]),
                             (r * H_th[2] - H_obs[2])])
        prod = (Cinv.dot(Delta)).T
        chi_sq1 = Delta.dot(prod)
        chi_sq_11 = chi_sq1.item((0, 0))
        return chi2B + chi2H + chi2d + chi_sq_11 + chi2DV1 + chi2DV2


    def Li(O, Ok_0, rr, O_phi_z, h0, afin):
        l, error = integrate.quad(lambda H0: math.exp((-1 / 2) * (chi_sq(H0, O, Ok_0, rr, O_phi_z, h0, afin))
                                                      + (-1 / (2 * (s ** 2))) * ((H0 - H0_av) ** 2)),
                                  H0_av - n * s, H0_av + n * s, epsabs=0, epsrel=1.49e-08)  # Set error bounds?
        return l / (math.sqrt(2 * pi * (s ** 2)))


    def Li_2(alpha0, beta0):
        alpha = 1 / (s ** 2) + alpha0
        beta = H0_av / (s ** 2) + beta0
        V = (1 / 2) * ((alpha) * (s ** 2)) ** (-1 / 2)
        W = gamma - ((beta) ** 2) / alpha
        X = beta / (math.sqrt(2 * alpha))
        Y = 1 + scipy.special.erf(X)
        Z = math.exp((-1 / 2) * W)
        return V * Z * Y


    def Ofunc(d, K, k, sol, al):
        Omegam1 = (4. / 9.) * (1. / (sol[d, 2]) ** 3.)
        Omegaphi1 = (1. / 12.) * (((sol[d, 1]) ** 2.) + k / ((sol[d, 0]) ** al))
        Omegak1 = -K / ((sol[d, 2]) ** 2.)
        Omegam = Omegam1 / (Omegam1 + Omegaphi1 + Omegak1)
        return Omegam, Omegam1, Omegaphi1, Omegak1

    def calc_main(K):
        subO1 = []
        suba1 = []
        subOk1 = []
        subLike = []

        K = round(K, 2)
        for al in arange(0.4, 2.0, 0.01):
            al = round(al, 2)
            print(K, al)
            k = (8 / 3) * ((al + 4) / (al + 2)) * ((2 / 3) * (al * (al + 2))) ** (
                    al / 2)  # This is kappa, from eq. (2) of arXiv:1307.7399v1.

            t0 = 10 ** (-4)
            p0 = (((2 / 3) * (al * (al + 2))) ** (1 / 2)) * t0 ** (2 / (al + 2))  # Initial value of phi.
            v0 = (((8 / 3) * al * (1 / (al + 2))) ** (1 / 2)) / t0 ** (al / (al + 2))  # Initial value of d(phi)/dt.
            a0 = t0 ** (2 / 3)  # I assumed a ~ t^(2/3) in the early universe.

            t = numpy.arange(0, 150, t0)

            # arrays of initial conditions and parameters
            w0 = [p0, v0, a0]
            zz = [al, k, m, K]

            # solution array
            sol = odeint(phiCDM, w0, t, args=(zz,))
            sol_col_2 = sol[:int(150 / t0), 2]
            # for H0 in arange(50, 85.1, 0.1):
            for O in arange(0.2, 0.4, 0.01):
                # start_time = time.time()
                # print(O)
                for b in range(0, int(150 / t0), 1):
                    if O >= Ofunc(b, K, k, sol, al)[0]:
                        break

                h0 = (sol[b + 1, 2] - sol[b, 2]) / (sol[b, 2])
                afin = sol[b, 2]
                Omegam1, Omegaphi1, Omegak1 = Ofunc(b, K, k, sol, al)[1], Ofunc(b, K, k, sol, al)[2], Ofunc(b, K, k, sol, al)[3]
                Ok_0 = Omegak1 / (Omegam1 + Omegaphi1 + Omegak1)

                rr = []
                O_phi_z = []
                for q in range(0, 8, 1):
                    z1 = z[q]
                    r = 0
                    qq = 0
                    qt = 0
                    total = 0.0
                    trapezoid_prev_d = 0
                    trapezoid_prev_y = 0
                    trapezoid_cur_y = 0
                    z1_bound = 1 / (1 + z1) * afin  # Optimized loop check
                    first_elem = False
                    for trapezoid_d in range(0, b + 1, 1):  # Manually calculating trapezoidal sum
                        if (sol_col_2[trapezoid_d] >= z1_bound):
                            trapezoid_cur_y = 1 / sol_col_2[trapezoid_d]
                            r += trapezoid_cur_y
                            qq += 1
                            qt += 1
                            if first_elem:
                                total += ((trapezoid_cur_y + trapezoid_prev_y) / 2) * (trapezoid_d - trapezoid_prev_d)
                            else:
                                first_elem = True
                                trapezoid_prev_d = trapezoid_d
                            trapezoid_prev_y = trapezoid_cur_y
                        if qq == 1:
                            Omegaphi_z = (1 / 12) * (((sol[trapezoid_d, 1]) ** 2) + k / ((sol[trapezoid_d, 0]) ** al))
                    if qt == 1:
                        rr.append(r)
                    else:
                        rr.append(total)
                    O_phi_z.append(Omegaphi_z / (Omegam1 + Omegak1 + Omegaphi1))
                alpha0 = 0
                beta0 = 0

                for q in range(0, len(z_obs), 1):
                    afin_bound = afin / (1 + z_obs[q])  # Optimized loop check
                    for g in range(0, int(150 / t0), 1):
                        if (sol_col_2[g] >= afin_bound):
                            Omegaphi2 = (1 / 12) * (((sol[g, 1]) ** 2) + k / ((sol[g, 0]) ** al))
                            O_phi_z1 = Omegaphi2 / (Omegam1 + Omegak1 + Omegaphi1)
                            break
                    alpha0 += (E(O, z_obs[q], Ok_0, O_phi_z1) ** 2) / ((sigHobs[q]) ** 2)
                    beta0 += (Hz_obs[q]) * (E(O, z_obs[q], Ok_0, O_phi_z1)) / ((sigHobs[q]) ** 2)

                subO1.append(O) # each process adds its data to a separate sublist
                suba1.append(al)
                subOk1.append(Ok_0)
                subLike.append(Li(O, Ok_0, rr, O_phi_z, h0, afin) * Li_2(alpha0, beta0))

        return (subO1, suba1, subOk1, subLike) # each process returns a set of sublists


    if __name__ == '__main__':
        pool = multiprocessing.Pool()
        var_range = arange(-0.45, 0.46, 0.01)
        result = pool.map(calc_main, var_range) # the set of sublists, already in order

        for subO1, suba1, subOk1, subLike in result: # appends each process' sublist to the respective full list
            O1.extend(subO1)
            a1.extend(suba1)
            Ok1.extend(subOk1)
            Like.extend(subLike)

        print(Ok1)
        print(O1)
        print(a1)
        print(Like)
        # print(O_phi_z)

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

if __name__ == '__main__':
    print("--- %s seconds ---" % (time.time() - start_time))
