 # -*- coding: utf-8 -*-
"""
@author: t.maybour
"""

import numpy as np
import enh_solvers as f2p
from scipy import optimize

class enh_class():

#################################################################################
#################################################################################
#################################################################################

    def __init__(self):

        # physical constants
        self.c = np.float64(299792458) # speed of light
        self.e0 = np.float64(8.85418782e-12)
        self.u0 = np.float64(1.25663706e-6)
        self.nu = np.sqrt(self.u0 / self.e0)

    ###########################################
    ###########################################
    ###########################################

    def sellmeier(self, l0):

        l0 = l0*1e6

        ## lithium niobate
        a1 = 5.756
        a2 = 0.0983
        a3 = 0.202
        a4 = 189.32
        a5 = 12.52
        a6 = 1.32E-02
        b1 = 2.86E-06
        b2 = 4.70E-08
        b3 = 6.11E-08
        b4 = 1.52E-04

        T = 21
        f = (T - 24.5)*(T + 570.82)

        ne_sqr = a1 + b1*f + (a2 + b2*f) / (l0**2 - (a3 + b3*f)**2) + (a4 + b4*f) / (l0**2 - a5**2) - a6*l0**2
        return np.sqrt(ne_sqr)

    ###########################################
    ###########################################
    ###########################################

    def setup(self):

        ###########################################
        ###########################################
        ###########################################

        self.xv_nopo_input = 1000
        self.res_bare = 10
        self.centre_shift = 1
        self.pig = 1
        self.pout = 0

        ###########################################
        ###########################################
        ###########################################

        self.error_threshold = 1e-8
        self.thres_safety = 1e-2
        self.max_rlevel = 100

        ###########################################
        ###########################################
        ###########################################

        self.phis = 0
        self.phi1 = np.pi/2

        self.lB = self.l0
        self.lQ = self.l0

        ###########################################
        ###########################################
        ###########################################

        self.k0 = np.float64((2 * np.pi) / self.l0)
        self.w0 = np.float64(self.c * self.k0)
        self.f0 = np.float64(self.w0 / (2 * np.pi))

        self.n1 = self.sellmeier(self.l0)
        self.n2 = self.sellmeier(self.l0/2)

        ###########################################
        ###########################################
        ###########################################

        self.kB = np.float64((2 * np.pi) / self.lB)
        self.wB = np.float64(self.c * self.kB)
        self.fB = np.float64(self.wB / (2 * np.pi))

        self.n1B = self.sellmeier(self.lB)
        self.n2B = self.sellmeier(self.lB/2)

        ###########################################
        ###########################################
        ###########################################

        self.kQ = np.float64((2 * np.pi) / self.lQ)
        self.wQ = np.float64(self.c * self.kQ)
        self.fQ = np.float64(self.wQ / (2 * np.pi))

        self.n1Q = self.sellmeier(self.lQ)
        self.n2Q = self.sellmeier(self.lQ/2)

        ###########################################
        ###########################################
        ###########################################

        self.int2amp_fun = np.sqrt((2*self.intensity*1e4)/(self.c*self.n1*self.e0)) # (V/m)
        self.int2amp_har = np.sqrt((2*self.intensity*1e4)/(self.c*self.n2*self.e0)) # (V/m)

        self.ampsqr2int_fun = self.c*self.n1*self.e0/(2*1e4)
        self.ampsqr2int_har = self.c*self.n2*self.e0/(2*1e4)

        ###########################################
        ###########################################
        ###########################################

        self.beta1 = (self.n1 * self.w0) / self.c
        self.beta2 = (2 * self.n2 * self.w0) / self.c
        self.delta_beta = self.beta2-2*self.beta1

        self.scale = (self.x2 * self.w0) / self.c
        self.pm_scale = 2*self.scale/np.pi

        ###########################################
        ###########################################
        ###########################################

        self.dB = np.zeros(2, dtype=np.float64)
        self.dB[0] = self.l0 / (2*self.n1)
        self.dB[1] = self.l0 / (4*self.n2)

        ###########################################
        ###########################################
        ###########################################

        self.d1 = self.lB / (2*self.n1)
        self.d2 = self.lB / (4*self.n2)

        self.detuning = np.zeros(2, dtype=np.float64)
        self.detuning[0] = self.beta1 - np.pi/self.d1
        self.detuning[1] = self.beta2 - np.pi/self.d2

        ###########################################
        ###########################################
        ###########################################

        self.check_save = False
        self.fail = False

    def set_params(self):
        self.kappa = np.zeros(4, dtype=np.complex128)
        self.kappa[0] = np.exp(1j*self.phi1) * (self.dn * self.w0) / (2 * self.c)
        self.kappa[1] = np.exp(-1j*self.phi1) * (self.dn * self.w0) / (2 * self.c)
        self.kappa[2] = np.exp(1j*self.phi2) * (self.dn * self.w0) / self.c
        self.kappa[3] = np.exp(-1j*self.phi2) * (self.dn * self.w0) / self.c
        self.alpha = np.pi /  self.ds

    def xv_setup(self):

        ###########################################
        ############## set postions ###############
        ###########################################

        min_lengths = []
        min_lengths.append(1e-5)
        min_lengths.append(1e-3*self.xlen)

        self.dx_calc_bare = np.min(min_lengths)
        self.l0_res = self.l0 / self.dx_calc_bare

        self.xv_nopo_calc = int(self.xlen / self.dx_calc_bare) + 1
        self.xv_calc = np.linspace(0, self.xlen, self.xv_nopo_calc, dtype=np.float64)
        self.dx_calc = self.xv_calc[1]-self.xv_calc[0]

        self.xv_int = int(self.xv_nopo_calc / self.xv_nopo_input) + 1
        self.xv_nopo = self.xv_calc[::self.xv_int].shape[0]
        self.xv = np.zeros(self.xv_nopo, dtype=np.float64)
        self.xv[:] = self.xv_calc[::self.xv_int]

        self.xv_nopos_calc = np.linspace(0, self.xv_nopo_calc-1, self.xv_nopo_calc, dtype=np.uint32)
        self.xv_nopos = np.zeros(self.xv_nopo, dtype=np.uint32)
        self.xv_nopos[:] = self.xv_nopos_calc[::self.xv_int]

        ###########################################
        ###########################################
        ###########################################

        self.env_moire = np.asfortranarray(np.zeros((2,self.xv_nopo), dtype=np.complex128))
        self.env_x2_moire = np.asfortranarray(np.zeros((4,self.xv_nopo), dtype=np.complex128))

        self.y0 = np.zeros(4, dtype=np.complex128)
        self.y0_2 = np.zeros(2, dtype=np.complex128)

    ###########################################
    ###########################################
    ###########################################

    def norm2param(self, norm, ll, ul):
        return (ul-ll)*norm + ll

##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################

    def run_sim(self):

        int_hold = self.intensity
        loop_count = 50

        loop = True
        recur = False
        count = 0
        int_initial = self.intensity

        str_pass = 'ds '+'{:.3e}'.format(self.ds)+' '

        while loop:

            count += 1
            self.intensity = int_initial
            self.int2amp_fun = np.sqrt((2*int_initial*1e4)/(self.c*self.n1*self.e0))
            self.int2amp_har = np.sqrt((2*int_initial*1e4)/(self.c*self.n2*self.e0)) # (V/m)

            self.run_moire()
            self.x2_moire_cw_fwd([np.real(self.env_moire[1,0]),np.imag(self.env_moire[1,0]),0,0])

            if self.error < self.error_threshold:
                if self.pout:
                    print(str(int_initial) + ', error: ' + str(self.error))
                break
            else:
                recur = True
                int_initial = int_initial / 10
                if self.pout:
                    print(str(int_initial) + ', error: ' + str(self.error))
            if count > loop_count:
                break

        if recur:
            self.fail = False
            self.run_x2_moire_cw_fwd_int(str_pass,self.initial,int_initial,int_hold)

        elif self.pout:
            print(str_pass+'error: '+'{:.3e}'.format(self.error))

        ###########################################
        ###########################################

        self.ratio_fun_fwd = self.ampsqr2int_fun*np.abs(self.env_x2_moire[0,-1])**2/self.intensity
        self.ratio_fun_bwd = self.ampsqr2int_fun*np.abs(self.env_x2_moire[1,0])**2/self.intensity
        self.ratio_har_fwd = self.ampsqr2int_har*np.abs(self.env_x2_moire[2,-1])**2/self.intensity
        self.ratio_har_bwd = self.ampsqr2int_har*np.abs(self.env_x2_moire[3,0])**2/self.intensity

        if self.pout:
            print('------------------------------------')
            print('env x2 moire ratio fund fwd: ' + str(self.ratio_fun_fwd))
            print('env x2 moire ratio fund bwd: ' + str(self.ratio_fun_bwd))
            print('env x2 moire ratio 2har fwd: ' + str(self.ratio_har_fwd))
            print('env x2 moire ratio 2har bwd: ' + str(self.ratio_har_bwd))
            print('env x2 moire ratio sum: ' + str(self.ratio_fun_fwd+self.ratio_fun_bwd+self.ratio_har_fwd+self.ratio_har_bwd))
            print('------------------------------------')

##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################

    def run_moire(self):

        self.y0_2[0] = 1
        self.y0_2[1] = 0

        f2p.lib.rk4_cmt_moire_bwd(kappa=np.array([self.kappa[0],self.kappa[1]]),alpha=self.alpha,detuning=self.detuning[0],phi=self.phis,cs=self.centre_shift,pig=self.pig,ascale=self.ascale,xlen=self.xlen,x0=self.xv_calc[-1],dx_calc=self.dx_calc,y0=self.y0_2,xv_nopo=self.xv_nopo,xv_nopos=self.xv_nopos,env=self.env_moire)

        self.env_moire[:,:] = self.int2amp_fun * self.env_moire[:,:] / self.env_moire[0,0]
        self.env_moire[:,:]

##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################

    def error_func(self):

        #######################
        #######################

        # Normalised Intensities
        nint_fund_fwd_start = self.ampsqr2int_fun*np.abs(self.env_x2_moire[0,0])**2/self.intensity
        nint_fund_fwd_end = self.ampsqr2int_fun*np.abs(self.env_x2_moire[0,-1])**2/self.intensity

        nint_fund_bwd_start = self.ampsqr2int_fun*np.abs(self.env_x2_moire[1,0])**2/self.intensity
        nint_fund_bwd_end = self.ampsqr2int_fun*np.abs(self.env_x2_moire[1,-1])**2/self.intensity

        nint_2har_fwd_start = self.ampsqr2int_har*np.abs(self.env_x2_moire[2,0])**2/self.intensity
        nint_2har_fwd_end = self.ampsqr2int_har*np.abs(self.env_x2_moire[2,-1])**2/self.intensity

        nint_2har_bwd_start = self.ampsqr2int_har*np.abs(self.env_x2_moire[3,0])**2/self.intensity
        nint_2har_bwd_end = self.ampsqr2int_har*np.abs(self.env_x2_moire[3,-1])**2/self.intensity

        #######################
        #######################

        # Input amplitude error
        self.error1 = np.abs(self.int2amp_fun - self.env_x2_moire[0,0])**2

        # Boundary conditions error
        self.error2 = self.ampsqr2int_har*np.abs(self.env_x2_moire[2,0])**2/self.intensity
        self.error3 = self.ampsqr2int_fun*np.abs(self.env_x2_moire[1,-1])**2/self.intensity
        self.error4 = self.ampsqr2int_har*np.abs(self.env_x2_moire[3,-1])**2/self.intensity

        #######################
        #######################

        # total error
        error = self.error1 + self.error2 + self.error3 + self.error4

        #######################
        #######################

        if np.isnan(error):
            return 1e10
        else:
            return error

##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################

    def x2_moire_cw_fwd(self, initial):
        self.initial = self.x2_moire_cw_fwd_error_min(initial)
        self.error = self.x2_moire_cw_fwd_error_func(self.initial)

    ##################################
    ##################################

    def x2_moire_cw_fwd_error_min(self, initial):
        minimum = optimize.fmin(self.x2_moire_cw_fwd_error_func, x0=initial, xtol=1e-25, ftol=1e-25, disp=False)
        return minimum

    ##################################
    ##################################

    def x2_moire_cw_fwd_error_func(self, initial):

        self.y0[0] = self.int2amp_fun
        self.y0[1] = initial[0] + 1j*initial[1]
        self.y0[2] = 0
        self.y0[3] = initial[2] + 1j*initial[3]

        f2p.lib.rk4_cmt_x2_moire_cw_fwd(n1=self.n1,n2=self.n2,pm_scale=self.pm_scale,filter=self.filter,kappa=self.kappa,alpha=self.alpha,phi=self.phis,cs=self.centre_shift,pig=self.pig,ascale=self.ascale,xlen=self.xlen,x0=self.xv_calc[0],dx_calc=self.dx_calc,y0=self.y0,xv_nopo=self.xv_nopo,xv_nopos=self.xv_nopos,env=self.env_x2_moire)

        return self.error_func()

    ##################################
    ##################################

    def run_x2_moire_cw_fwd_int(self, str_pass, initial, ints, inte):

        self.initial = self.x2_moire_cw_fwd_int_recur(ints, inte, initial, False, 0, str_pass+'|')
        self.error = self.x2_moire_cw_fwd_error_func(self.initial)

    ##################################
    ##################################

    def x2_moire_cw_fwd_int_recur(self, ints, inte, initial, fail, level, level_str):

        level_str_extend = level_str+'|'

        first = True
        while True:

            if level > self.max_rlevel:
                self.fail = True

            self.intensity = inte
            self.int2amp_fun = np.sqrt((2*self.intensity*1e4)/(self.c*self.n1*self.e0)) # (V/m)
            self.int2amp_har = np.sqrt((2*self.intensity*1e4)/(self.c*self.n2*self.e0)) # (V/m)

            initial_hold = initial
            initial = self.x2_moire_cw_fwd_error_min(initial_hold)
            error = self.x2_moire_cw_fwd_error_func(initial)

            thres = self.error_threshold
            if first:
                thres = thres*self.thres_safety

            if error > thres:
                level_str_extend = level_str+'0|'
            else:
                level_str_extend = level_str+'1|'

            if self.pout:
                print(str(level_str_extend)+' r'+str(level+1)+', ethres: '+'{:.2e}'.format(self.error_threshold)+', error: '+'{:.2e}'.format(error)+', ints: '+'{:.5e}'.format(ints)+', inte: '+'{:.5e}'.format(inte)+', epsilon: '+'{:.5e}'.format(inte-ints))

            if self.fail:
                break

            if error > thres:
                initial = self.x2_moire_cw_fwd_int_recur(ints, (ints+inte)/2, initial_hold, fail, level + 1, level_str_extend)
                ints = (ints+inte)/2
                first = False
            else:
                break

        self.intensity = inte
        self.int2amp_fun = np.sqrt((2*self.intensity*1e4)/(self.c*self.n1*self.e0)) # (V/m)
        self.int2amp_har = np.sqrt((2*self.intensity*1e4)/(self.c*self.n2*self.e0)) # (V/m)

        return initial

##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
