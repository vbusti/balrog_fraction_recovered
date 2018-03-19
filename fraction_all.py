#!/usr/bin/env python
# encoding: UTF8
#

from __future__ import print_function
import numpy as np
import pyfits as pf
import matplotlib.pyplot as plt
from matplotlib import cm
import os

realizations = ['0','1','2']
bands = ['g','r','i','z']
tiles = ['DES0347-5540','DES2329-5622','DES2357-6456']

for real in realizations:

    for tile in tiles:

        for band in bands:

            # check if there is a folder for the results, if not create it

            if not os.path.isdir('../results/'+real):
                os.mkdir('../results/'+real)

            if not os.path.isdir('../results/'+real+'/'+tile):
                os.mkdir('../results/'+real+'/'+tile)

            if not os.path.isdir('../results/'+real+'/'+tile+'/'+band):
                os.mkdir('../results/'+real+'/'+tile+'/'+band)


            # data plus injections

            cat_inj = pf.open("/home/vinicius/Documents/balrog/fraction_recovered/data/"+real+"/"+tile+"/"+tile+"_"+band+"_cat.fits")
            cat_inj.info()
            tab_inj = cat_inj[1].data
            my_format_inj = tab_inj.formats
            print(tab_inj.names)
            print(my_format_inj)

            # truth catalogs of injections

            cat_true = pf.open("/home/vinicius/Documents/balrog/fraction_recovered/data/"+real+"/"+tile+"/"+tile+"_"+real+"_balrog_truth_cat.fits")
            cat_true.info()
            tab_true = cat_true[1].data
            my_format_true = tab_true.formats
            print(tab_true.names)
            print(my_format_true)

            ra_i  = tab_inj['ALPHAWIN_J2000']
            dec_i = tab_inj['DELTAWIN_J2000']
            mag_i = tab_inj['MAG_AUTO']
            mag_i_err = tab_inj['MAGERR_AUTO']

            ra_t  = tab_true['ra']
            dec_t = tab_true['dec']
            mag_t = tab_true['cm_mag']

            if(band == 'g'):
                ind_band = 0
            elif(band == 'r'):
                ind_band = 1
            elif(band == 'i'):
                ind_band = 2
            else:
                ind_band = 3       

            mag_t_band = np.array([mag_t[i][ind_band] for i in range(len(mag_t))])

            print(np.min(mag_t_band),np.max(mag_t_band))

            min_mag = np.min(mag_t_band)
            max_mag = 26.#np.max(mag_t_band)

            mag_bins = np.linspace(min_mag-0.01,max_mag+0.01,np.rint((max_mag-min_mag+0.02)/0.25)+1)
            print(mag_bins)


            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            plt.hist(mag_t_band)
            plt.title('Magnitude Distribution of the Injected Galaxies',size=16)
            plt.xlabel(band+' [mag]',size=16)
            plt.ylabel('# of galaxies',size=16)
            plt.savefig('../results/'+real+'/'+tile+'/'+band+'/mag_'+band+'_dist_injected.png')
            plt.close()


            flag_match = np.zeros(len(ra_t))
            flag_blend = np.zeros(len(ra_t))
            matched_elements_all = [[] for _ in range(len(ra_t))]
            matched_elements_min = [[] for _ in range(len(ra_t))]

            for i in range(len(ra_t)):
                mini_mask = (np.abs(ra_i-ra_t[i]) < 0.01)*(np.abs(dec_i-dec_t[i])<0.01)
                ra_temp = ra_i[mini_mask]
                dec_temp = dec_i[mini_mask]
                mag_temp = mag_i[mini_mask]
                mag_err_temp = mag_i_err[mini_mask]
  
 
                for j in range(len(ra_temp)):
                    cos_theta = np.sin(np.pi*dec_t[i]/180.)*np.sin(np.pi*dec_temp[j]/180.) + np.cos(np.pi*dec_t[i]/180.)*np.cos(np.pi*dec_temp[j]/180.)*np.cos(np.pi*(ra_t[i]-ra_temp[j])/180.) 
                    theta = 206300.*np.arccos(cos_theta)
                    if(theta < 1.):
                        flag_match[i] = 1
                        matched_elements_all[i].append([mag_temp[j],mag_err_temp[j],theta])
                if(flag_match[i] == 0):
                    matched_elements_min[i].append([[-999,-999,-999]])
                elif(flag_match[i] == 1 and len(np.array(matched_elements_all[i]) < 3)):
                    matched_elements_min[i].append(matched_elements_all[i])
                else:
                    print('hi ',i)
                    print(matched_elements_all[i])
                    flag_blend[i] = 1
                    z = np.array([matched_elements_all[i][k][2] for k in range(len(matched_elements_all[i]))])
                    ind_min, = np.where(z == np.min(z))  
                    matched_elements_min[i].append(matched_elements_all[i][ind_min[0]]) 
                    print(matched_elements_min[i])            
     
            mean_mag_bin  = np.zeros(len(mag_bins)-1)
            frac_bin      = np.zeros(len(mag_bins)-1)
            mask_zero      = np.zeros(len(mag_bins)-1)

            for i in range(len(mag_bins)-1):
                mask_mag         = (mag_t_band >= mag_bins[i])*(mag_t_band < mag_bins[i+1])               
                mag_temp         = mag_t_band[mask_mag]
                flag_temp        = flag_match[mask_mag]
                mean_mag_bin[i]  = np.mean(mag_temp)
                if(len(flag_temp) > 0):
                    frac_bin[i] = np.sum(flag_temp)/np.float(len(flag_temp))
                    mask_zero[i] = 1  
                else:
                    frac_bin[i] = 2

            mask = (mask_zero == 1)

            mean_mag_bin = mean_mag_bin[mask]
            frac_bin     = frac_bin[mask]

            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            plt.scatter(mean_mag_bin,frac_bin)
            plt.title('Magnitude versus the Fraction of Recovered Galaxies',size=16)
            plt.xlabel(band+' [mag]',size=16)
            plt.ylabel('Fraction',size=16)
            plt.savefig('../results/'+real+'/'+tile+'/'+band+'/mag_'+band+'_vs_fraction.png')
            plt.close()

            np.savetxt('../results/'+real+'/'+tile+'/'+band+'/mag_fraction.txt',np.array([mean_mag_bin,frac_bin]).T)

            mag_t_band = np.array(mag_t_band[flag_match == 1])
            mag_inj    = np.array([matched_elements_min[i][0][0][0] for i in range(len(matched_elements_min))])
            mag_inj    = mag_inj[flag_match == 1]

            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            plt.scatter(mag_t_band,mag_t_band - mag_inj,s=0.2)
            plt.title('Magnitude versus Offset of Recovered Magnitudes',size=16)
            plt.xlabel(band+' [mag]',size=16)
            plt.ylabel('$\Delta m$',size=16)
            plt.ylim(-2,2)
            plt.savefig('../results/'+real+'/'+tile+'/'+band+'/mag_'+band+'_vs_deltamag.png')
            plt.close()

            np.savetxt('../results/'+real+'/'+tile+'/'+band+'/mag_deltamag.txt',np.array([mag_t_band,mag_t_band - mag_inj]).T)

            delta_m = mag_t_band - mag_inj
            mask_d  = (delta_m > -2.)*(delta_m < 2)
            delta_m = delta_m[mask_d]
            mean_dm = np.mean(delta_m)
            std_dm  = np.std(delta_m)
            print(mean_dm,std_dm)

            np.savetxt('../results/'+real+'/'+tile+'/'+band+'/mag_deltamag.txt',np.array([mean_dm,std_dm]))







