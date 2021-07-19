# -*- coding: utf-8 -*-
"""
@author: t.maybour
"""
import os
import numpy as np
import enh_class as enh_class
import random

def main():

    enh = enh_class.enh_class()

##########################################
############### OPTO PARAMS ##############
##########################################

    # ascale_norm = 0.183744347
    # ds_norm = 0.828538613
    # phi2_norm = 0.524100845
    # filter_norm = 0.98831637

    ascale_norm = 0.9679022420584498
    ds_norm = 0.7588612239142497
    phi2_norm = 0.5276869342547482
    filter_norm = 0.34131973631201884

##########################################
################# setup ##################
##########################################

    #######################
    enh.xlen = 4e-2
    enh.intensity = 1e3 # Intensity (W/cm^2)
    enh.l0 = 1.064E-06
    enh.x2 = 25e-12
    enh.dn = 1e-3
    #######################

    ####
    ascale_ll = 8
    ascale_ul = 24
    ####
    ds_ll = 1e-3
    ds_ul = 7e-3
    ####
    phi2_ll = 0
    phi2_ul = np.pi
    ####
    filter_ll = 0
    filter_ul = 1.5e-3
    ####

###########################################
###########################################
###########################################

    enh.setup()
    enh.xv_setup()

###########################################
###########################################
###########################################

    enh.ascale = enh.norm2param(ascale_norm, ascale_ll, ascale_ul)
    enh.ds = enh.norm2param(ds_norm, ds_ll, ds_ul)
    enh.phi2 = enh.norm2param(phi2_norm, phi2_ll, phi2_ul)
    enh.filter = enh.norm2param(filter_norm, filter_ll, filter_ul)
    enh.set_params()

    print('ascale:',enh.ascale)
    print('ds:',enh.ds)
    print('phi2:',enh.phi2)
    print('filter:',enh.filter)

###########################################
###########################################
###########################################

    enh.pout = 0
    enh.run_sim()

    print('------------------------------------')
    print('ratio fund fwd:',enh.ratio_fun_fwd)
    print('ratio fund bwd:',enh.ratio_fun_bwd)
    print('ratio 2har fwd:',enh.ratio_har_fwd)
    print('ratio 2har bwd:',enh.ratio_har_bwd)
    print('norm power:',enh.ratio_fun_fwd+enh.ratio_fun_bwd+enh.ratio_har_fwd+enh.ratio_har_bwd)
    print('------------------------------------')

#################################################################################
#################################################################################
#################################################################################

if __name__ == '__main__':
    main()
