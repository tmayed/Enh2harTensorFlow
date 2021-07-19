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
################# setup ##################
##########################################

    #######################
    enh.xlen = 4e-2
    enh.intensity = 1e3 # Intensity (W/cm^2)
    enh.l0 = 1.064E-06
    enh.x2 = 25e-12
    enh.dn = 1e-3
    #######################

    #####
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

    datafile = 'sim_data_work.csv'
    if os.path.isfile(datafile):
        os.remove(datafile)
    open(datafile, 'w').close()

    with open(datafile, 'a') as file:
        file.write('ascale_norm,ds_norm,phi2_norm,filter_norm,har_fwd,power_norm'+'\n')

###########################################
###########################################
###########################################

    enh.setup()
    enh.xv_setup()

###########################################
###########################################
###########################################

    for ii in range(250):

    ###########################################
    ###########################################
    ###########################################

        print("sim no. :",ii)
        ascale_norm = random.uniform(0, 1)
        ds_norm = random.uniform(0, 1)
        phi2_norm = random.uniform(0, 1)
        filter_norm = random.uniform(0, 1)

        enh.ascale = enh.norm2param(ascale_norm, ascale_ll, ascale_ul)
        enh.ds = enh.norm2param(ds_norm, ds_ll, ds_ul)
        enh.phi2 = enh.norm2param(phi2_norm, phi2_ll, phi2_ul)
        enh.filter = enh.norm2param(filter_norm, filter_ll, filter_ul)

        enh.set_params()
        enh.run_sim()

        with open(datafile, 'a') as file:
            file.write(str(ascale_norm)+','+str(ds_norm)+','+str(phi2_norm)+','+str(filter_norm)+','+\
            str(enh.ratio_har_fwd)+','+str(enh.ratio_fun_fwd+enh.ratio_fun_bwd+enh.ratio_har_fwd+enh.ratio_har_bwd)+'\n')

#################################################################################
#################################################################################
#################################################################################

if __name__ == '__main__':
    main()
