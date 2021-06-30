import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
from scipy.io.idl import readsav
from datetime import datetime
import matplotlib.dates as mdates
import os
from astropy.time import Time
import pickle
import scipy.io
import math

eventMassMeasured = [['20100203_A', 1.45*1e15, 0.10*1e15, float('Nan'), float('Nan'), float('Nan')],
                    ['20100203_B', 1.45*1e15, 0.10*1e15, float('Nan'), float('Nan'), float('Nan')],
                    ['20100319_A', 1.50*1e15, 0.15*1e15, float('Nan'), float('Nan'), float('Nan')],
                    ['20100319_B', 2.20*1e15, 0.15*1e15, float('Nan'), float('Nan'), float('Nan')],
                    ['20100403_A', 3.90*1e15, 0.10*1e15, float('Nan'), float('Nan'), 6.04*1e15],
                    ['20100403_B', 3.30*1e15, 0.15*1e15, float('Nan'), float('Nan'), 6.04*1e15],
                    ['20100408_A', 4.10*1e15, 0.10*1e15, float('Nan'), float('Nan'), 7.52*1e15],
                    ['20100408_B', 3.50*1e15, 0.10*1e15, float('Nan'), float('Nan'), 7.52*1e15],
                    ['20100523_A', 3.10*1e15, 0.10*1e15, float('Nan'), float('Nan'), float('Nan')],
                    ['20100523_B', 3.20*1e15, 0.10*1e15, float('Nan'), float('Nan'), float('Nan')],
                    ['20101026_A', 5.35*1e15, 0.10*1e15, float('Nan'), float('Nan'), 7.14*1e15],
                    ['20101026_B', 4.70*1e15, 0.20*1e15, float('Nan'), float('Nan'), 7.14*1e15],
                    ['20110130_A', 2.00*1e15, 0.15*1e15, float('Nan'), float('Nan'), 5.60*1e15],
                    ['20110130_B', 2.00*1e15, 0.15*1e15, float('Nan'), float('Nan'), 5.60*1e15],
                    ['20110214_A', 6.85*1e15, 0.05*1e15, float('Nan'), float('Nan'), 6.83*1e15],
                    ['20110214_B', 6.20*1e15, 0.10*1e15, float('Nan'), float('Nan'), 6.83*1e15],
                    ['20110906_A', 5.00*1e15, 0.10*1e15, float('Nan'), float('Nan'), float('Nan')],
                    ['20110906_B', 2.50*1e15, 0.10*1e15, float('Nan'), float('Nan'), float('Nan')],
                    ['20120123_A', 1.20*1e16, 0.10*1e16, 1.80*1e16, 0.20*1e16, float('Nan')],
                    ['20120123_B', 6.50*1e15, 0.20*1e15, 6.50*1e15, 0.20*1e15, float('Nan')],
                    ['20120712_A', 9.60*1e15, 0.10*1e15, 1.65*1e16, 0.15*1e16, 1.84*1e16],
                    ['20120712_B', 2.70*1e15, 0.10*1e15, 1.15*1e16, 0.15*1e16, 1.84*1e16]]

eventMassMeasured = np.array(eventMassMeasured)

def getcat(filename):
#    print('reading CAT ' + filename)
    print(filename)
    cat = scipy.io.readsav(filename)  # , verbose='false')
#    print('done reading CAT')
    return cat

# ############################################################################
# ################################ main program ##############################
# ############################################################################
# ##################################### CONTROLS

# directory of current event
# Animation of ensemble simulations for ElEvoHI

# Author: T. Amerstorfer, J. Hinterreiter IWF Graz, Austria
# November 2018
# This work is published under the MIT LICENSE (see bottom)
# parameter:
#            scriptPath: Path to the ELEvoHI runs

def main(scriptPath):

#    rcParams['font.size'] = 15
    # run='CMEs/Feb2019BGinsitu'
    #run = 'HI_animate/events/Paper_BGSW/'
    run = scriptPath
    
    BGSWmodel = run.split('/')[len(run.split('/'))-2]

    # path = 'VisualizeCMEPropagation/' + run + '/'
    path = run + '/'
    eventdate = next(os.walk(path))[1]
    eventdate.sort()

    outfolder = path + '/plots/'
    if os.path.exists(outfolder) is False:
            os.makedirs(path + '/plots/')

    scs = next(os.walk(path + eventdate[0] + '/results/'))[1]
    scs.sort()
    
    sc = 'Earth'
    if True:
        meanDT = []
        events = []
        predmeds = []
        for j in range(0, len(eventdate)):
            sc_folder = path + eventdate[j] + '/results/'+sc+'/'
            
#            print(path + eventdate[j] + '/all_apex_variables.p')
            pickleFile = path + eventdate[j] + '/all_apex_variables.p'
            if os.path.exists(pickleFile) is True:
                [all_apex_t, all_apex_r, all_apex_lat, all_apex_lon, all_apex_f,
                 all_apex_w, all_apex_s, all_apex_run, all_apex_flag,
                 CME_start_time, dur_days, startcutFit, endcutFit, frontDataKins] = pickle.load(
                 open(pickleFile, "rb"))
            else:
                frontFile = path + eventdate[j] + '/results/frontDataAll.sav'
                if os.path.exists(frontFile) is True:
                    frontData = getcat(frontFile)
                    for i in np.arange(0, len(frontData.frontkins.arrtimeearth)):
                        if frontData.frontkins.arrtimeearth[i] != b'      -1':
                            frontData.frontkins.arrtimeearth[i] = mdates.date2num(
                            Time.strptime(frontData.frontkins.arrtimeearth[i], '%Y-%m-%dT%H:%M:%S.%f').datetime)
                        else:
                            frontData.frontkins.arrtimeearth[i] = float('Nan')
        
                    frontDataKins = frontData.frontkins

            if os.path.exists(sc_folder + '/prediction.sav') is True:
                elevohi_results = readsav(path + eventdate[j] +
                                          '/eELEvoHI_results.sav', verbose=0)


                data = readsav(sc_folder + '/gamma.sav', verbose=0)
                data1 = readsav(sc_folder + '/rinit.sav', verbose=0)
                data2 = readsav(sc_folder + '/sw.sav', verbose=0)
                data3 = readsav(sc_folder + '/vinit.sav', verbose=0)
                tt = readsav(sc_folder + '/transittimes.sav', verbose=0)
                arrival = readsav(sc_folder + '/prediction.sav', verbose=0)
                arr = readsav(sc_folder + '/arrivaltimes.sav', verbose=0)
                plotarr = readsav(sc_folder + '/plottimes.sav', verbose=0)
                tt_colors = readsav(sc_folder + '/labels.sav', verbose=0)
                issc = sc
                if arr.insitu_sc.decode() == 'B':
                    issc = 'STEREO-B'
                elif arr.insitu_sc.decode() == 'A':
                    issc = 'STEREO-A'
                elif arr.insitu_sc.decode() == 'Earth':
                    issc = 'Earth'
                elif arr.insitu_sc.decode() == 'MES':
                    issc = 'MESSENGER'
                elif arr.insitu_sc.decode() == 'VEX':
                    issc = 'Venus Express'
                else:
                    print('No in situ spacecraft defined!')

                count_nan = 0
                arrtime = []
                for i in range(len(frontDataKins.arrtimeearth)):
                    if frontDataKins.arrtimeearth[i] != float('Nan'):
                        arrtime.append(frontDataKins.arrtimeearth[i])
                    else:
                        count_nan = count_nan + 1
  
                count_nan = 0
                arrtime = []
                for at in frontDataKins.arrtimeearth:
                    if np.isnan(at):
                        count_nan = count_nan + 1
                    else:
                        arrtime.append(at)

#                print("Number of predicted misses: ", count_nan, "from", 
#                      len(frontDataKins.arrtimeearth), " total predictions.")
    
#                print(len(arrtime))
#                count_nan = 0
#                arrtime = []
                # arrtime=np.zeros(np.size(0))

 #               for i in range(len(arr.arrivaltimes)):
 #                   if arr.arrivaltimes[i].decode() != "NaN":
 #                       arrtime.append(Time.strptime(
 #                           arr.arrivaltimes[i].decode(), '%Y-%m-%dT%H:%M').datetime)
 #                   else:
 #                       count_nan = count_nan + 1

 #               print("Number of predicted misses: ", count_nan, "from", len(
 #                   arr.arrivaltimes), " total predictions.")
                
                # convert bytes to strings

 #               arrtime_ = np.zeros(np.size(arrtime))

#                for i in range(len(arrtime)):
#                    arrtime_[i] = mdates.date2num(arrtime[i])

                colormap = plt.cm.viridis
                col = np.arange(8)
                
                arrtime_ = arrtime

                for i in range(8):
                    col[i] = i * 35

                # print(arr.insitu_sc.decode())

                # four panel figure with frequency distributions of output parameter

                fig_fourpanels, ax = plt.subplots(
                    2, 2, sharex=False, sharey=False, figsize=(10, 8))
                #fig_fourpanels, ax = plt.subplots(
                #    2, 3, sharex=False, sharey=False, figsize=(16, 10))

                # plt.rc('text', usetex=True)
                # plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
                rcParams['font.family'] = 'sans-serif'
                rcParams['font.sans-serif'] = ['Tahoma']

                labels = tt_colors.labels.astype('U13')
                # ['$54 - 58$ h','$58 - 64$ h','$64 - 69$ h','$69 - 73$ h','$74 -
                # 78$ h','$78 - 83$ h','$83 - 88$ h','$88 - 90$ h']
                
                binBoundaries1 = np.linspace(min(data.gamma), max(data.gamma), 8)
                dat = [data.gamma1, data.gamma2, data.gamma3, data.gamma4,
                       data.gamma5, data.gamma6, data.gamma7, data.gamma8]
                
                ax[0, 0].hist(dat[7], bins=binBoundaries1, color='lightblue', ec='gray', zorder=2)
                ax[0, 0].set_ylabel("Number of runs", fontsize=14)
                ax[0, 0].set_xlabel(
                    "($\pm$) Drag Parameter [10$^{-7}$ km$^{-1}$]", fontsize=12)
                
                #######

                binBoundaries2 = np.linspace(min(data1.rinit), max(data1.rinit), 8)
                dat = [data1.rinit1, data1.rinit2, data1.rinit3, data1.rinit4,
                       data1.rinit5, data1.rinit6, data1.rinit7, data1.rinit8]
                ax[0, 1].hist(dat[7], bins=binBoundaries2, color='lightblue', ec='gray', zorder=2)
                ax[0, 1].set_ylabel("Number of runs", fontsize=14)
                ax[0, 1].set_xlabel("Initial Distance [R$_\odot$]", fontsize=12)
            
                print('len: ', len(dat[7]))
            
                #######

                binBoundaries3 = np.linspace(min(data3.vinit), max(data3.vinit), 8)
                dat = [data3.vinit1, data3.vinit2, data3.vinit3, data3.vinit4,
                       data3.vinit5, data3.vinit6, data3.vinit7, data3.vinit8]
                
                ax[1, 1].hist(dat[7], bins=binBoundaries3, color='lightblue', ec='gray', zorder=2)
                ax[1, 1].set_ylabel("Number of runs", fontsize=14)
                ax[1, 1].set_xlabel("Initial Speed [km s$^{-1}$]", fontsize=12)
                
                
                #######
                
                binBoundaries4 = np.linspace(min(data2.sw), max(data2.sw), 8)
                dat = [data2.sw1, data2.sw2, data2.sw3, data2.sw4, data2.sw5,
                       data2.sw6, data2.sw7, data2.sw8]

                ax[1, 0].hist(dat[7], bins=binBoundaries4, color='lightblue', ec='gray', zorder=2)
                ax[1, 0].set_ylabel("Number of runs", fontsize=14)
                ax[1, 0].set_xlabel(
                    "Ambient Solar Wind Speed [km s$^{-1}$]", fontsize=12)

                #######
                

                fig_fourpanels.subplots_adjust(hspace=0.2)
                fig_fourpanels.suptitle(
                    'ELEvoHI Ensemble Prediction for ' + issc, fontsize=16)

                # fig_fourpanels.show()

                fig_fourpanels.savefig(outfolder + '/' + eventdate[j] + '_' + sc + '_results.pdf')
                fig_fourpanels.clf()
                
                plt.close('all')
                
                
                
                
                
                # ############################deformable front plot################
#                fig_fourpanels, ax = plt.subplots(
#                    2, 3, sharex=False, sharey=False, figsize=(16, 10))
                fig_fourpanels, ax = plt.subplots(
                    2, 3, sharex=False, sharey=False, figsize=(13, 7))
#                plt.xticks(fontsize=14)
#                plt.yticks(fontsize=14)
                
                dragsOri = frontDataKins.dragparameter
#                drags = np.log10(dragsOri)
                drags = dragsOri
                binBoundaries1 = np.linspace(min(drags), max(drags), 8)
                ax[0, 0].hist(drags, bins=binBoundaries1, color='lightblue', ec='gray', zorder=2)
                ax[0, 0].hist(dat[7], bins=binBoundaries1, color='lightblue', ec='gray', zorder=2)
                ax[0, 0].axvline(x = dragsOri[0], color='r')
                ax[0, 0].set_ylabel("Number of runs")#, fontsize=16)
                ax[0, 0].set_xlabel(
                    #"Initial Drag parameter [km$^{-1}$]", fontsize=12)
                    "Drag Parameter: ($\pm$)$\gamma$ [km$^{-1}$]")#, fontsize=16)
                
                print('len(drags): ', len(drags))
                
                distMass = frontDataKins.distmass
                binBoundaries2 = np.linspace(min(distMass), max(distMass), 8)
                ax[0, 1].hist(distMass, bins=binBoundaries2, color='lightblue', ec='gray', zorder=2)
                ax[0, 1].axvline(x = distMass[0], color='r')
                ax[0, 1].set_ylabel("Number of runs")#, fontsize=16)
                ax[0, 1].set_xlabel("Distance [R$_\odot$]")#, fontsize=16)
            
            
            
                defSpeed = []
                for allVs in frontDataKins.varr:
                    defSpeed.append(allVs[int(len(allVs)/2)][0])
                binBoundaries3 = np.linspace(min(defSpeed), max(defSpeed), 8)

                ax[0, 2].hist(defSpeed, bins=binBoundaries3, color='lightblue', ec='gray', zorder=2)
                ax[0, 2].axvline(x = defSpeed[0], color='r')
                ax[0, 2].set_ylabel("Number of runs")#, fontsize=16)
                ax[0, 2].set_xlabel(
                     "CME speed [km s$^{-1}$]")#, fontsize=16)

                   
                defDens = frontDataKins.rho
                binBoundaries3 = np.linspace(min(defDens), max(defDens), 8)

                ax[1, 0].hist(defDens, bins=binBoundaries3, color='lightblue', ec='gray', zorder=2)
                ax[1, 0].axvline(x = defDens[0], color='r')
                ax[1, 0].set_ylabel("Number of runs")#, fontsize=16)
                ax[1, 0].set_xlabel(
                     "Ambient Solar Wind Density [g km$^{-3}$]")#, fontsize=16)
                
                
                areasOri = frontDataKins.area
                areas = np.log10(areasOri)
#                areas = areasOri
                binBoundaries6 = np.linspace(min(areas), max(areas), 8)
                ax[1, 1].hist(areas, bins=binBoundaries6, color='lightblue', ec='gray', zorder=2)
                ax[1, 1].axvline(x = areas[0], color='r')
                ax[1, 1].set_ylabel("Number of runs")#, fontsize=16)
                ax[1, 1].set_xlabel(
                    "Cross-sectional Area: log(A [km$^2$])")#, fontsize=16)
                
                
                indEvent = np.where(eventMassMeasured[:,0] == eventdate[j])[0][0]
                
                measMass1 = float(eventMassMeasured[indEvent][1])
                errMass1 = float(eventMassMeasured[indEvent][2])
                measMass2 = float(eventMassMeasured[indEvent][3])
                errMass2 = float(eventMassMeasured[indEvent][4])
                measMass3 = float(eventMassMeasured[indEvent][5])
                
                               
                massesOri = frontDataKins.mass
                masses = np.log10(massesOri)
#                masses = massesOri
                binBoundaries5 = np.linspace(min(masses), max(masses), 8)
                ax[1, 2].hist(masses, bins=binBoundaries5, color='lightblue', ec='gray', zorder=2)
                ax[1, 2].axvline(x = masses[0], color='r')
                ax[1, 2].set_ylabel("Number of runs")#, fontsize=16)
                ax[1, 2].set_xlabel(
                    #"CME mass [g]", fontsize=12)
                    "CME mass: log(m [g])")#, fontsize=16)
                
                
                shadelow1 = np.log10(measMass1-errMass1)
                shadehigh1 = np.log10(measMass1+errMass1)
                measMass1 = np.log10(measMass1)
                                
                ylim1 = ax[1,2].get_ylim()[0]
                ylim2 = ax[1,2].get_ylim()[1]
 
                ax[1, 2].axvline(x = measMass1, color='blue')
                ax[1, 2].fill_between([shadelow1, shadehigh1],
                                [ax[1,2].get_ylim()[0], ax[1,2].get_ylim()[0]],
                                [ax[1,2].get_ylim()[1], ax[1,2].get_ylim()[1]],
                                facecolor='lightblue', zorder=4, alpha=0.3,
                                edgecolor='none')
                
                if not np.isnan(measMass2):
                    shadelow2 = np.log10(measMass2-errMass2)
                    shadehigh2 = np.log10(measMass2+errMass2)
                    measMass2 = np.log10(measMass2)
 
                    ax[1, 2].plot([measMass2, measMass2], [ylim1,
                                           ylim2], color='green')
                    ax[1, 2].fill_between([shadelow2, shadehigh2],
                                    [ylim1,ylim1],
                                    [ylim2,ylim2],
                                    facecolor='lightgreen', zorder=3, alpha=0.3,
                                    edgecolor='none')
                if not np.isnan(measMass3):
                    measMass3 = np.log10(measMass3)

                    ax[1, 2].plot([measMass3, measMass3], [ylim1,
                                           ylim2], color='red')
                
                    
                ax[1,2].set_ylim([ylim1, ylim2])
                fig_fourpanels.subplots_adjust(hspace=0.2)
                fig_fourpanels.suptitle(
                    'ELEvoHI/'+BGSWmodel+'\nInput Parameters for Deformable Front', fontsize=16)

                # fig_fourpanels.show()
                fig_fourpanels.savefig(outfolder + '/' + eventdate[j] + '_' + sc + '_defFront_inputs.pdf')
                fig_fourpanels.clf()
                
                plt.close('all')
                
              

                
                # ############################prediction plot######################


                predmedian = plotarr.arrplotmedian
                predmean = plotarr.arrplotmean
                prederr = plotarr.arrerr
                

                insitu = ''

                if plotarr.insituarr != '':
                    insitu = plotarr.insituarr.decode()
                if insitu != '':
                    num_insitu = mdates.date2num(Time.strptime(
                        plotarr.insituarr.decode(), '%Y-%m-%d %H:%M').datetime)

                arri = np.median(arrtime)
                arri2 = np.mean(arrtime)
                shadelow = arri - np.std(arrtime)
                shadehigh = arri + np.std(arrtime)
                
                predmedian = mdates.num2date(arri).strftime("%Y-%m-%d %H:%M")
                predmean = mdates.num2date(arri2).strftime("%Y-%m-%d %H:%M")
                prederr = str(np.round(np.std(arrtime)*24, 1))
     
 
                fig, ax = plt.subplots(1, 1, figsize=(8, 5))
                ax.hist(arrtime_, bins=10, color='lightblue', ec='gray', zorder=2)

                ymin = 0
                ymax = ax.get_ylim()[1] + 1
                ax.set_ylim(ymin, ymax)

                ax.plot([arri, arri], [ax.get_ylim()[0],
                                       ax.get_ylim()[1]], color='blue')
                ax.fill_between([shadelow, shadehigh],
                                [ax.get_ylim()[0], ax.get_ylim()[0]],
                                [ax.get_ylim()[1], ax.get_ylim()[1]],
                                facecolor='slategrey', zorder=3, alpha=0.3,
                                edgecolor='none')
                ax.plot([arri2, arri2], [ax.get_ylim()[0], ax.get_ylim()[1]],
                        color='green')
                if insitu != '':
                    ax.plot([num_insitu, num_insitu], [ax.get_ylim()[0],
                                                       ax.get_ylim()[1]], color='red')

                # ax.xaxis_date()

                ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=12))
                ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%b-%d\n%H:%M"))
                # ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
                xmin = ax.get_xlim()[0] - 0.05
                xmax = ax.get_xlim()[1] + 0.05
                ax.set_xlim(xmin, xmax)

                ax.xaxis_date()

                ax.set_ylabel("Number of Runs", fontsize=14)
                ax.set_xlabel(datetime.strftime(mdates.num2date(arri), '%Y %b %d') +
                              "\nPredicted Shock Arrival Time [UT]", fontsize=12)
                # ax.set_xlabel("Predicted Shock Arrival Time [UT]", fontsize=12)

                plt.grid(color='white', zorder=1)
                ax.set_facecolor('whitesmoke')

                plt.text(1.15, 0.8, 'Arrival Time [UT]:', horizontalalignment='center',
                         verticalalignment='center', transform=ax.transAxes)
                plt.text(1.15, 0.7, 'Median \n' + predmedian,
                         horizontalalignment='center', verticalalignment='center',
                         transform=ax.transAxes, color='blue')
                plt.text(1.15, 0.6, 'Mean \n' + predmean,
                         horizontalalignment='center', verticalalignment='center',
                         transform=ax.transAxes, color='green')
                plt.text(1.15, 0.45, 'Standard Deviation\n$\pm$' +
                         prederr + ' hours', horizontalalignment='center',
                         verticalalignment='center', transform=ax.transAxes,
                         color='black', bbox={'facecolor': 'slategray', 'alpha': 0.2,
                                              'pad': 10})

                if insitu != '':
                    plt.text(1.15, 0.2, 'Insitu Arrival \n' + insitu,
                             horizontalalignment='center', verticalalignment='center',
                             transform=ax.transAxes, color='red')

                SCHit = 0
                NoHit = 0
                if sc == 'Earth':
                    arrtime_sc = elevohi_results.eelevohi.arrtime_earth[0]
                if sc == 'MES':
                    arrtime_sc = elevohi_results.eelevohi.arrtime_mes[0]
                if sc == 'VEX':
                    arrtime_sc = elevohi_results.eelevohi.arrtime_vex[0]
                if sc == 'STA':
                    arrtime_sc = elevohi_results.eelevohi.arrtime_sta[0]
                if sc == 'STB':
                    arrtime_sc = elevohi_results.eelevohi.arrtime_stb[0]
                if sc == 'SOLO':
                    arrtime_sc = elevohi_results.eelevohi.arrtime_solo[0] 
                if sc == 'PSP':
                    arrtime_sc = elevohi_results.eelevohi.arrtime_psp[0] 
                NrRuns = np.size(arrtime_sc)
                for i in range(0, NrRuns):
                    if 'NaN' in str(arrtime_sc[i]):
                        NoHit = NoHit + 1
                    else:
                        SCHit = SCHit + 1

                plt.text(1.15, 0.1, str(SCHit) + '/' + str(NrRuns) +
                         ' runs hit ' + issc, horizontalalignment='center',
                         verticalalignment='center', transform=ax.transAxes,
                         color='black')

                plt.title('ELEvoHI Ensemble Prediction for ' + issc, fontsize=16)

                fig.subplots_adjust(hspace=0.2)

                plt.xticks(rotation=45)

                # fig.show()

                fig.savefig(outfolder + '/' + eventdate[j] + '_' + sc + '_ELEvoHI_prediction.pdf',
                            bbox_inches="tight")
                fig.clf()

                plt.close('all')
                
                print(eventdate[j],':', predmedian, '   ', insitu, '=> ', (arri-num_insitu)*24, 'hrs')
                meanDT.append((arri-num_insitu)*24)
                
                events.append(eventdate[j])
                predmeds.append(predmedian)
                # input()

        print(meanDT)
        print('mean: ', np.mean(meanDT))
        print('MAE: ', np.mean(np.abs(meanDT)))
        print('median: ', np.median(meanDT))
        print('std: ', np.std(meanDT))
        print('RMSE: ', np.sqrt(np.sum(np.power(meanDT,2))/np.size(meanDT)) )
        
        pickle.dump((events, predmeds, meanDT), open(path + "/preds.p", "wb"))
        
        
    catPath = 'HI_animate/cats/ShockNormal/'
    [CMEDates, CMEArrTimes, shockAngles] = pickle.load(
            open(catPath + 'shockAngles.p', "rb"))

    path = 'HI_animate/events/DefFront/shortTracks/'

    eventdate = next(os.walk(path))[1]

    eventdate.sort()

    for evs in eventdate:
        eDate =  evs[0:len(evs)-3]
        sc = evs[len(evs)-1]

        file = path + evs + '/tpWind_AB.p'
        if os.path.exists(file):

            [tpWind_a, tpWind_b, et_time_num_interp_a, et_time_num_interp_b, tangAngle_a, 
                             tangAngle_b, tangentPoints_a, tangentPoints_b, shockAngles_a, shockAngles_b,
                             diffMeasure_a, diffMeasure_b] = pickle.load(
                    open(file, "rb"))


            shockAngleMeas = float('Nan')
            i = -1
            for CMED in CMEDates:
                i = i+1
                if CMED == eDate:
                    shockAngleMeas = shockAngles[i]


            fig, ax = plt.subplots(1, 2, figsize=(12, 6))
            fig.suptitle('CME Shock Normal')
            ax[0].set_title('STEREO-A')
            ax[0].hist(shockAngles_a, bins=10, color='lightblue', ec='gray', zorder=2)
            ax[0].set_xlabel('Shock Angle [°]', fontsize=14)
            ax[0].set_ylabel("Number of runs", fontsize=14)
            ax[1].set_title('STEREO-B')
            ax[1].hist(shockAngles_b, bins=10, color='lightblue', ec='gray', zorder=2)
            ax[1].set_xlabel('Shock Angle [°]', fontsize=14)
            ax[1].set_ylabel("Number of runs", fontsize=14)

            ylim1 = ax[0].get_ylim()[0]
            ylim2 = ax[0].get_ylim()[1]
            ax[0].plot([shockAngleMeas, shockAngleMeas], [ylim1, ylim2], color='red')
            ylim1 = ax[1].get_ylim()[0]
            ylim2 = ax[1].get_ylim()[1]
            ax[1].plot([shockAngleMeas, shockAngleMeas], [ylim1, ylim2], color='red')

            fig.savefig(path + '/plots/' + eDate + '_shockAngles.png',
                bbox_inches="tight")
            fig.clf()

            plt.close('all')
            

if __name__ == '__main__':
    # eventslist = ['20090623']
    main('/nas/helio/ELEvoHI_plotting/runs/')

# ########################## MIT license


# Copyright 2018 Dr. Tanja Amerstorfer
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions: The above copyright
# notice and this permission notice shall be included in all copies or
# substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS",
# WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED
# TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
# FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
# THE USE OR OTHER DEALINGS IN THE SOFTWARE.