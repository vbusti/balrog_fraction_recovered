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
INPUT_DIR = 'blank_test'
OUT_DIR = 'blank_test/no_flags/galaxies2'
OBJECT = 'galaxies' # 'galaxies' or 'stars'
RADIUS_MATCH = 1.0

REMOVE_FLAGGED = True





for real in realizations:

    for tile in tiles:

        # truth catalogs of injections of galaxies

        cat_true = pf.open("/home/vinicius/Documents/balrog/fraction_recovered/data/"+INPUT_DIR+"/"+real+"/"+tile+"/"+tile+"_"+real+"_balrog_truth_cat_gals.fits")
        cat_true.info()
        tab_true = cat_true[1].data
        my_format_true = tab_true.formats

        ra_t  = tab_true['ra']
        dec_t = tab_true['dec']
        mag_t = tab_true['cm_mag']

        # truth catalogs of injections of stars

        cat_trues = pf.open("/home/vinicius/Documents/balrog/fraction_recovered/data/"+INPUT_DIR+"/"+real+"/"+tile+"/"+tile+"_"+real+"_balrog_truth_cat_stars.fits")
        cat_trues.info()
        tab_trues = cat_trues[1].data
        my_format_trues = tab_trues.formats

        ra_ts  = tab_trues['RA_new']
        dec_ts = tab_trues['DEC_new']
        mag_gs = tab_trues['g_Corr']
        mag_rs = tab_trues['g_Corr'] - tab_trues['gr_Corr']
        mag_is = tab_trues['g_Corr'] - tab_trues['gr_Corr'] - tab_trues['ri_Corr']
        mag_zs = tab_trues['g_Corr'] - tab_trues['gr_Corr'] - tab_trues['ri_Corr'] - tab_trues['iz_Corr']
 

        flag_match_all = [[] for _ in range(len(bands))]
        mag_all        = [[] for _ in range(len(bands))]

        for band in bands:

            # check if there is a folder for the results, if not create it

            path = '../results/'+OUT_DIR+'/'+real+'/'+tile+'/'+band
            if not os.path.isdir(path):
                os.makedirs(path)
            
            
            # data plus injections

            cat_inj = pf.open("/home/vinicius/Documents/balrog/fraction_recovered/data/"+INPUT_DIR+"/"+real+"/"+tile+"/"+tile+"_"+band+"_cat.fits")
            cat_inj.info()
            tab_inj = cat_inj[1].data
            my_format_inj = tab_inj.formats
                        
            ra_i  = tab_inj['ALPHAWIN_J2000']
            dec_i = tab_inj['DELTAWIN_J2000']
            mag_i = tab_inj['MAG_AUTO']
            mag_i_err = tab_inj['MAGERR_AUTO']
            flags = tab_inj['FLAGS']

            if(REMOVE_FLAGGED == True):
                mask = (flags == 0)    
                ra_i  = ra_i[mask]
                dec_i = dec_i[mask]
                mag_i = mag_i[mask]
                mag_i_err = mag_i_err[mask]
  

            if(band == 'g'):
                ind_band = 0
                mag_star = mag_gs
            elif(band == 'r'):
                ind_band = 1
                mag_star = mag_rs 
            elif(band == 'i'):
                ind_band = 2
                mag_star = mag_is
            else:
                ind_band = 3     
                mag_star = mag_zs  

            mag_t_band = np.array([mag_t[i][ind_band] for i in range(len(mag_t))])

            if(OBJECT == 'stars'):
                ra_t  = ra_ts
                dec_t = dec_ts
                mag_t_band = mag_star  
            

            min_mag = np.min(mag_t_band)
            max_mag = np.max(mag_t_band)

            mag_bins = np.linspace(min_mag-0.01,max_mag+0.01,np.rint((max_mag-min_mag+0.02)/0.25)+1)


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
                    if(theta < RADIUS_MATCH):
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
                 
            
            
            flag_match_all[ind_band].append(flag_match)
                 
            mean_mag_bin   = np.zeros(len(mag_bins)-1)
            frac_bin       = np.zeros(len(mag_bins)-1)
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
            plt.gcf().subplots_adjust(left=0.20)
            plt.scatter(mean_mag_bin,frac_bin)
            plt.title('Magnitude versus the fraction of recovered '+OBJECT,size=16)
            plt.xlabel(band+' [mag]',size=16)
            plt.ylabel('Fraction',size=16)
            plt.savefig('../results/'+OUT_DIR+'/'+real+'/'+tile+'/'+band+'/mag_'+band+'_vs_fraction.png')
            plt.close()

            np.savetxt('../results/'+OUT_DIR+'/'+real+'/'+tile+'/'+band+'/mag_fraction.txt',np.array([mean_mag_bin,frac_bin]).T)

            mag_t_band = np.array(mag_t_band[flag_match == 1])
            mag_inj    = np.array([matched_elements_min[i][0][0][0] for i in range(len(matched_elements_min))])
            mag_all[ind_band].append(mag_inj)
            print('mag_all')
            print(mag_all[ind_band])   
            mag_inj    = mag_inj[flag_match == 1]

            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            plt.gcf().subplots_adjust(left=0.20) 
            plt.scatter(mag_t_band,mag_t_band - mag_inj,s=0.2)
            plt.title('Magnitude versus Offset of Recovered Magnitudes',size=16)
            plt.xlabel(band+' [mag]',size=16)
            plt.ylabel('$\Delta m$',size=16)
            #plt.ylim(-0.1,0.1)
            plt.savefig('../results/'+OUT_DIR+'/'+real+'/'+tile+'/'+band+'/mag_'+band+'_vs_deltamag.png')
            plt.close()

            np.savetxt('../results/'+OUT_DIR+'/'+real+'/'+tile+'/'+band+'/mag_deltamag.txt',np.array([mag_t_band,mag_t_band - mag_inj]).T)

            delta_m = mag_t_band - mag_inj
            mask_d  = (delta_m > -2.)*(delta_m < 2)
            delta_m = delta_m[mask_d]
            mean_dm = np.mean(delta_m)
            std_dm  = np.std(delta_m)
            print(mean_dm,std_dm)

            np.savetxt('../results/'+OUT_DIR+'/'+real+'/'+tile+'/'+band+'/mag_deltamag.txt',np.array([mean_dm,std_dm]))
        
        if(OBJECT == 'galaxies'):
            # g-r 
            color_det = np.zeros(len(ra_t))
            gmr       = np.zeros(len(ra_t))
            gmr_t     = np.zeros(len(ra_t)) 
            for i in range(len(ra_t)):
                if(flag_match_all[0][0][i] == 1 and flag_match_all[1][0][i] == 1):
                    color_det[i] = 1
                    gmr[i]    = mag_all[0][0][i] - mag_all[1][0][i]
                    gmr_t[i]  = mag_t[i][0] - mag_t[i][1] 

            mask_c = (color_det == 1)
            gmr    = gmr[mask_c]
            gmr_t  = gmr_t[mask_c]

            # check if there is a folder for the results, if not create it

            path = '../results/'+OUT_DIR+'/'+real+'/'+tile+'/colors'
            if not os.path.isdir(path):
                os.makedirs(path)
            

            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            plt.gcf().subplots_adjust(left=0.20) 
            plt.scatter(gmr_t,gmr_t-gmr,s=0.4)
            plt.title('g-r color',size=16)
            plt.xlabel('g-r true',size=16)
            plt.ylabel('$\Delta_{g-r}$',size=16)
            #plt.xlim(-3,3) 
            #plt.ylim(-3,3)
            plt.savefig('../results/'+OUT_DIR+'/'+real+'/'+tile+'/colors/gmr.png')
            plt.close() 

            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            plt.gcf().subplots_adjust(left=0.20)
            plt.scatter(gmr_t,gmr_t-gmr,s=0.4)
            plt.title('g-r color',size=16)
            plt.xlabel('g-r true',size=16)
            plt.ylabel('$\Delta_{g-r}$',size=16)
            plt.xlim(-0.5,2) 
            plt.ylim(-0.1,0.1)
            plt.savefig('../results/'+OUT_DIR+'/'+real+'/'+tile+'/colors/gmr_zoom.png')
            plt.close() 



            # r-i 
            color_det = np.zeros(len(ra_t))
            gmr       = np.zeros(len(ra_t))
            gmr_t     = np.zeros(len(ra_t)) 
            for i in range(len(ra_t)):
                if(flag_match_all[1][0][i] == 1 and flag_match_all[2][0][i] == 1):
                    color_det[i] = 1
                    gmr[i]    = mag_all[1][0][i] - mag_all[2][0][i]
                    gmr_t[i]  = mag_t[i][1] - mag_t[i][2] 

            mask_c = (color_det == 1)
            gmr    = gmr[mask_c]
            gmr_t  = gmr_t[mask_c]

            # check if there is a folder for the results, if not create it

            #path = '../results/'+OUT_DIR+'/'+real+'/'+tile+'/colors'
            #if not os.path.isdir(path):
            #os.makedirs(path)

            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            plt.gcf().subplots_adjust(left=0.20)  
            plt.scatter(gmr_t,gmr_t-gmr,s=0.4)
            plt.title('r-i color',size=16)
            plt.xlabel('r-i true',size=16)
            plt.ylabel('$\Delta_{r-i}$',size=16)
            #plt.xlim(-3,3) 
            #plt.ylim(-3,3)
            plt.savefig('../results/'+OUT_DIR+'/'+real+'/'+tile+'/colors/rmi.png')
            plt.close() 

            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            plt.scatter(gmr_t,gmr_t-gmr,s=0.4)
            plt.gcf().subplots_adjust(left=0.20)
            plt.title('r-i color',size=16)
            plt.xlabel('r-i true',size=16)
            plt.ylabel('$\Delta_{r-i}$',size=16)
            plt.xlim(-0.5,1.5)
            plt.ylim(-0.1,0.1)
            plt.savefig('../results/'+OUT_DIR+'/'+real+'/'+tile+'/colors/rmi_zoom.png')
            plt.close() 

            # i-z 
            color_det = np.zeros(len(ra_t))
            gmr       = np.zeros(len(ra_t))
            gmr_t     = np.zeros(len(ra_t)) 
            for i in range(len(ra_t)):
                if(flag_match_all[2][0][i] == 1 and flag_match_all[3][0][i] == 1):
                    color_det[i] = 1
                    gmr[i]    = mag_all[2][0][i] - mag_all[3][0][i]
                    gmr_t[i]  = mag_t[i][2] - mag_t[i][3] 

            mask_c = (color_det == 1)
            gmr    = gmr[mask_c]
            gmr_t  = gmr_t[mask_c]

            # check if there is a folder for the results, if not create it

            #path = '../results/'+OUT_DIR+'/'+real+'/'+tile+'/colors'
            #if not os.path.isdir(path):
            #os.makedirs(path)
            

            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            plt.gcf().subplots_adjust(left=0.20)
            plt.scatter(gmr_t,gmr_t-gmr,s=0.4)
            plt.title('i-z color',size=16)
            plt.xlabel('i-z true',size=16)
            plt.ylabel('$\Delta_{r-i}$',size=16)
            #plt.xlim(-3,3) 
            #plt.ylim(-3,3)
            plt.savefig('../results/'+OUT_DIR+'/'+real+'/'+tile+'/colors/imz.png')
            plt.close() 

            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            plt.gcf().subplots_adjust(left=0.20)
            plt.scatter(gmr_t,gmr_t-gmr,s=0.4)
            plt.title('i-z color',size=16)
            plt.xlabel('i-z true',size=16)
            plt.ylabel('$\Delta_{i-z}$',size=16)
            plt.xlim(-0.5,1.5) 
            plt.ylim(-0.1,0.1)
            plt.savefig('../results/'+OUT_DIR+'/'+real+'/'+tile+'/colors/imz_zoom.png')
            plt.close() 


        else:
            # g-r 
            color_det = np.zeros(len(ra_t))
            gmr       = np.zeros(len(ra_t))
            gmr_t     = np.zeros(len(ra_t)) 
            for i in range(len(ra_t)):
                if(flag_match_all[0][0][i] == 1 and flag_match_all[1][0][i] == 1):
                    color_det[i] = 1
                    gmr[i]    = mag_all[0][0][i] - mag_all[1][0][i]
                    gmr_t[i]  = mag_gs[i] - mag_rs[i] 

            mask_c = (color_det == 1)
            gmr    = gmr[mask_c]
            gmr_t  = gmr_t[mask_c]

            # check if there is a folder for the results, if not create it

            path = '../results/'+OUT_DIR+'/'+real+'/'+tile+'/colors'
            if not os.path.isdir(path):
                os.makedirs(path)
            

            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            plt.scatter(gmr_t,gmr_t-gmr,s=0.4)
            plt.gcf().subplots_adjust(left=0.20)  
            plt.title('g-r color',size=16)
            plt.xlabel('g-r true',size=16)
            plt.ylabel('$\Delta_{g-r}$',size=16)
            #plt.xlim(-3,3) 
            #plt.ylim(-3,3)
            plt.savefig('../results/'+OUT_DIR+'/'+real+'/'+tile+'/colors/gmr.png')
            plt.close() 

            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            plt.gcf().subplots_adjust(left=0.20)  
            plt.scatter(gmr_t,gmr_t-gmr,s=0.4)
            plt.title('g-r color',size=16)
            plt.xlabel('g-r true',size=16)
            plt.ylabel('$\Delta_{g-r}$',size=16)
            plt.xlim(-0.5,2) 
            plt.ylim(-0.2,0.2)
            plt.savefig('../results/'+OUT_DIR+'/'+real+'/'+tile+'/colors/gmr_zoom.png')
            plt.close() 



            # r-i 
            color_det = np.zeros(len(ra_t))
            gmr       = np.zeros(len(ra_t))
            gmr_t     = np.zeros(len(ra_t)) 
            for i in range(len(ra_t)):
                if(flag_match_all[1][0][i] == 1 and flag_match_all[2][0][i] == 1):
                    color_det[i] = 1
                    gmr[i]    = mag_all[1][0][i] - mag_all[2][0][i]
                    gmr_t[i]  = mag_rs[i] - mag_is[i] 

            mask_c = (color_det == 1)
            gmr    = gmr[mask_c]
            gmr_t  = gmr_t[mask_c]

            

            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            plt.gcf().subplots_adjust(left=0.20) 
            plt.scatter(gmr_t,gmr_t-gmr,s=0.4)
            plt.title('r-i color',size=16)
            plt.xlabel('r-i true',size=16)
            plt.ylabel('$\Delta_{r-i}$',size=16)
            #plt.xlim(-3,3) 
            #plt.ylim(-3,3)
            plt.savefig('../results/'+OUT_DIR+'/'+real+'/'+tile+'/colors/rmi.png')
            plt.close() 

            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            plt.gcf().subplots_adjust(left=0.20)
            plt.scatter(gmr_t,gmr_t-gmr,s=0.4)
            plt.title('r-i color',size=16)
            plt.xlabel('r-i true',size=16)
            plt.ylabel('$\Delta_{r-i}$',size=16)
            plt.xlim(-0.5,2) 
            plt.ylim(-0.1,0.1)
            plt.savefig('../results/'+OUT_DIR+'/'+real+'/'+tile+'/colors/rmi_zoom.png')
            plt.close() 

            # i-z 
            color_det = np.zeros(len(ra_t))
            gmr       = np.zeros(len(ra_t))
            gmr_t     = np.zeros(len(ra_t)) 
            for i in range(len(ra_t)):
                if(flag_match_all[2][0][i] == 1 and flag_match_all[3][0][i] == 1):
                    color_det[i] = 1
                    gmr[i]    = mag_all[2][0][i] - mag_all[3][0][i]
                    gmr_t[i]  = mag_is[i] - mag_zs[i] 

            mask_c = (color_det == 1)
            gmr    = gmr[mask_c]
            gmr_t  = gmr_t[mask_c]

            

            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            plt.gcf().subplots_adjust(left=0.20)
            plt.scatter(gmr_t,gmr_t-gmr,s=0.4)
            plt.title('i-z color',size=16)
            plt.xlabel('i-z true',size=16)
            plt.ylabel('$\Delta_{r-i}$',size=16)
            #plt.xlim(-3,3) 
            #plt.ylim(-3,3)
            plt.savefig('../results/'+OUT_DIR+'/'+real+'/'+tile+'/colors/imz.png')
            plt.close() 

            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            plt.gcf().subplots_adjust(left=0.20) 
            plt.scatter(gmr_t,gmr_t-gmr,s=0.4)
            plt.title('i-z color',size=16)
            plt.xlabel('i-z true',size=16)
            plt.ylabel('$\Delta_{i-z}$',size=16)
            plt.xlim(-0.5,1) 
            plt.ylim(-0.1,0.1)
            plt.savefig('../results/'+OUT_DIR+'/'+real+'/'+tile+'/colors/imz_zoom.png')
            plt.close() 



    








