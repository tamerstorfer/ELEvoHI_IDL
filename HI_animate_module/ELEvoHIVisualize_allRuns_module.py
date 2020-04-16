import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
from scipy.io.idl import readsav
from datetime import datetime
#from sunpy.time import parse_time
#import sunpy
import matplotlib.dates as mdates
import os
from astropy.time import Time


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

    # run='CMEs/Feb2019BGinsitu'
    #run = 'HI_animate/events/Paper_BGSW/'
    run = scriptPath

    # path = 'VisualizeCMEPropagation/' + run + '/'
    path = run + '/'
    eventdate = next(os.walk(path))[1]
    eventdate.sort()

    outfolder = path + '/plots/'
    if os.path.exists(outfolder) is False:
            os.makedirs(path + '/plots/')

    # path2 = path + eventdate[0] + '/results/Earth/'
    # predictions = next(os.walk(path2))[1]

    sc_folder = path + eventdate[0] + '/results/Earth/'

    for j in range(0, len(eventdate)):
        sc_folder = path + eventdate[j] + '/results/Earth/'

        if os.path.exists(sc_folder + '/prediction.sav') is True:
            elevohi_results = readsav(path + eventdate[j] +
                                      '/eELEvoHI_results.sav', verbose=1)

            data = readsav(sc_folder + '/gamma.sav', verbose=1)
            data1 = readsav(sc_folder + '/rinit.sav', verbose=1)
            data2 = readsav(sc_folder + '/sw.sav', verbose=1)
            data3 = readsav(sc_folder + '/vinit.sav', verbose=1)
            tt = readsav(sc_folder + '/transittimes.sav', verbose=1)
            arrival = readsav(sc_folder + '/prediction.sav', verbose=1)
            arr = readsav(sc_folder + '/arrivaltimes.sav', verbose=1)
            plotarr = readsav(sc_folder + '/plottimes.sav', verbose=1)
            tt_colors = readsav(sc_folder + '/labels.sav', verbose=1)
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

            print('Prediction for in situ s/c:')
            print(arr.insitu_sc.decode())

            count_nan = 0
            arrtime = []
            # arrtime=np.zeros(np.size(0))

            for i in range(len(arr.arrivaltimes)):
                if arr.arrivaltimes[i].decode() != "NaN":
                    arrtime.append(Time.strptime(
                        arr.arrivaltimes[i].decode(), '%Y-%m-%dT%H:%M').datetime)
                else:
                    count_nan = count_nan + 1

            print("Number of predicted misses: ", count_nan, "from", len(
                arr.arrivaltimes), " total predictions.")

            # convert bytes to strings

            arrtime_ = np.zeros(np.size(arrtime))

            for i in range(len(arrtime)):
                arrtime_[i] = mdates.date2num(arrtime[i])

            print('mean transit time:')
            print((np.mean(tt.tt)))

            print('standard deviation transit time:')
            print((np.std(tt.tt)))

            colormap = plt.cm.viridis
            col = np.arange(8)

            for i in range(8):
                col[i] = i * 35

            # print(arr.insitu_sc.decode())

            # four panel figure with frequency distributions of output parameter

            fig_fourpanels, ax = plt.subplots(
                2, 2, sharex=False, sharey=False, figsize=(10, 8))

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

            ax[0, 0].hist(dat[7], bins=binBoundaries1, histtype='bar', align='mid',
                          color=colormap.colors[col[0]],
                          edgecolor=colormap.colors[col[4]], alpha=1,
                          label=labels[7])
            ax[0, 0].hist(dat[6], bins=binBoundaries1, histtype='bar', align='mid',
                          color=colormap.colors[col[1]],
                          edgecolor=colormap.colors[col[4]], alpha=1,
                          label=labels[6])
            ax[0, 0].hist(dat[5], bins=binBoundaries1, histtype='bar', align='mid',
                          color=colormap.colors[col[2]],
                          edgecolor=colormap.colors[col[4]], alpha=1,
                          label=labels[5])
            ax[0, 0].hist(dat[4], bins=binBoundaries1, histtype='bar', align='mid',
                          color=colormap.colors[col[3]],
                          edgecolor=colormap.colors[col[4]], alpha=1,
                          label=labels[4])
            ax[0, 0].hist(dat[3], bins=binBoundaries1, histtype='bar', align='mid',
                          color=colormap.colors[col[4]],
                          edgecolor=colormap.colors[col[4]], alpha=1,
                          label=labels[3])
            ax[0, 0].hist(dat[2], bins=binBoundaries1, histtype='bar', align='mid',
                          color=colormap.colors[col[5]],
                          edgecolor=colormap.colors[col[4]], alpha=1,
                          label=labels[2])
            ax[0, 0].hist(dat[1], bins=binBoundaries1, histtype='bar', align='mid',
                          color=colormap.colors[col[6]],
                          edgecolor=colormap.colors[col[4]], alpha=1,
                          label=labels[1])
            ax[0, 0].hist(dat[0], bins=binBoundaries1, histtype='bar', align='mid',
                          color=colormap.colors[col[7]],
                          edgecolor=colormap.colors[col[4]], alpha=1,
                          label=labels[0])

            ax[0, 0].legend(loc=2, title='Transit Time')
            ax[0, 0].set_ylabel("Number of runs", fontsize=14)
            ax[0, 0].set_xlabel(
                "Drag parameter [10$^{-7}$ km$^{-1}$]", fontsize=12)

            #######

            binBoundaries2 = np.linspace(min(data1.rinit), max(data1.rinit), 8)
            dat = [data1.rinit1, data1.rinit2, data1.rinit3, data1.rinit4,
                   data1.rinit5, data1.rinit6, data1.rinit7, data1.rinit8]

            ax[0, 1].hist(dat[7], bins=binBoundaries2, histtype='bar',
                          color=colormap.colors[col[0]],
                          edgecolor=colormap.colors[col[4]], alpha=1,
                          label=labels[7])
            ax[0, 1].hist(dat[6], bins=binBoundaries2, histtype='bar',
                          color=colormap.colors[col[1]],
                          edgecolor=colormap.colors[col[4]], alpha=1,
                          label=labels[6])
            ax[0, 1].hist(dat[5], bins=binBoundaries2, histtype='bar',
                          color=colormap.colors[col[2]],
                          edgecolor=colormap.colors[col[4]], alpha=1,
                          label=labels[5])
            ax[0, 1].hist(dat[4], bins=binBoundaries2, histtype='bar',
                          color=colormap.colors[col[3]],
                          edgecolor=colormap.colors[col[4]], alpha=1,
                          label=labels[4])
            ax[0, 1].hist(dat[3], bins=binBoundaries2, histtype='bar',
                          color=colormap.colors[col[4]],
                          edgecolor=colormap.colors[col[4]], alpha=1,
                          label=labels[3])
            ax[0, 1].hist(dat[2], bins=binBoundaries2, histtype='bar',
                          color=colormap.colors[col[5]],
                          edgecolor=colormap.colors[col[4]], alpha=1,
                          label=labels[2])
            ax[0, 1].hist(dat[1], bins=binBoundaries2, histtype='bar',
                          color=colormap.colors[col[6]],
                          edgecolor=colormap.colors[col[4]], alpha=1,
                          label=labels[1])
            ax[0, 1].hist(dat[0], bins=binBoundaries2, histtype='bar',
                          color=colormap.colors[col[7]],
                          edgecolor=colormap.colors[col[4]], alpha=1,
                          label=labels[0])

            ax[0, 1].set_ylabel("Number of runs", fontsize=14)
            ax[0, 1].set_xlabel("Initial Distance [R$_\odot$]", fontsize=12)

            #######

            binBoundaries3 = np.linspace(min(data3.vinit), max(data3.vinit), 8)
            dat = [data3.vinit1, data3.vinit2, data3.vinit3, data3.vinit4,
                   data3.vinit5, data3.vinit6, data3.vinit7, data3.vinit8]

            ax[1, 1].hist(dat[7], bins=binBoundaries3, histtype='bar',
                          color=colormap.colors[col[0]],
                          edgecolor=colormap.colors[col[4]],
                          alpha=1, label=labels[7])
            ax[1, 1].hist(dat[6], bins=binBoundaries3, histtype='bar',
                          color=colormap.colors[col[1]],
                          edgecolor=colormap.colors[col[4]], alpha=1,
                          label=labels[6])
            ax[1, 1].hist(dat[5], bins=binBoundaries3, histtype='bar',
                          color=colormap.colors[col[2]],
                          edgecolor=colormap.colors[col[4]], alpha=1,
                          label=labels[5])
            ax[1, 1].hist(dat[4], bins=binBoundaries3, histtype='bar',
                          color=colormap.colors[col[3]],
                          edgecolor=colormap.colors[col[4]], alpha=1,
                          label=labels[4])
            ax[1, 1].hist(dat[3], bins=binBoundaries3, histtype='bar',
                          color=colormap.colors[col[4]],
                          edgecolor=colormap.colors[col[4]], alpha=1,
                          label=labels[3])
            ax[1, 1].hist(dat[2], bins=binBoundaries3, histtype='bar',
                          color=colormap.colors[col[5]],
                          edgecolor=colormap.colors[col[4]], alpha=1,
                          label=labels[2])
            ax[1, 1].hist(dat[1], bins=binBoundaries3, histtype='bar',
                          color=colormap.colors[col[6]],
                          edgecolor=colormap.colors[col[4]], alpha=1,
                          label=labels[1])
            ax[1, 1].hist(dat[0], bins=binBoundaries3, histtype='bar',
                          color=colormap.colors[col[7]],
                          edgecolor=colormap.colors[col[4]], alpha=1,
                          label=labels[0])

            ax[1, 1].set_ylabel("Number of runs", fontsize=14)
            ax[1, 1].set_xlabel("Initial Speed [km s$^{-1}$]", fontsize=12)

            #######

            binBoundaries4 = np.linspace(min(data2.sw), max(data2.sw), 8)
            dat = [data2.sw1, data2.sw2, data2.sw3, data2.sw4, data2.sw5,
                   data2.sw6, data2.sw7, data2.sw8]

            ax[1, 0].hist(dat[7], bins=binBoundaries4, histtype='bar',
                          color=colormap.colors[col[0]],
                          edgecolor=colormap.colors[col[4]], alpha=1,
                          label=labels[7])
            ax[1, 0].hist(dat[6], bins=binBoundaries4, histtype='bar',
                          color=colormap.colors[col[1]],
                          edgecolor=colormap.colors[col[4]], alpha=1,
                          label=labels[6])
            ax[1, 0].hist(dat[5], bins=binBoundaries4, histtype='bar',
                          color=colormap.colors[col[2]],
                          edgecolor=colormap.colors[col[4]], alpha=1,
                          label=labels[5])
            ax[1, 0].hist(dat[4], bins=binBoundaries4, histtype='bar',
                          color=colormap.colors[col[3]],
                          edgecolor=colormap.colors[col[4]], alpha=1,
                          label=labels[4])
            ax[1, 0].hist(dat[3], bins=binBoundaries4, histtype='bar',
                          color=colormap.colors[col[4]],
                          edgecolor=colormap.colors[col[4]], alpha=1,
                          label=labels[3])
            ax[1, 0].hist(dat[2], bins=binBoundaries4, histtype='bar',
                          color=colormap.colors[col[5]],
                          edgecolor=colormap.colors[col[4]], alpha=1,
                          label=labels[2])
            ax[1, 0].hist(dat[1], bins=binBoundaries4, histtype='bar',
                          color=colormap.colors[col[6]],
                          edgecolor=colormap.colors[col[4]], alpha=1,
                          label=labels[1])
            ax[1, 0].hist(dat[0], bins=binBoundaries4, histtype='bar',
                          color=colormap.colors[col[7]],
                          edgecolor=colormap.colors[col[4]], alpha=1,
                          label=labels[0])

            # ax[1,0].legend(loc=2, title='Transit Time')

            ax[1, 0].set_ylabel("Number of runs", fontsize=14)
            ax[1, 0].set_xlabel(
                "Background Solar Wind Speed [km s$^{-1}$]", fontsize=12)

            #######

            fig_fourpanels.subplots_adjust(hspace=0.2)
            fig_fourpanels.suptitle(
                'ELEvoHI Ensemble Prediction for ' + issc, fontsize=16)

            # fig_fourpanels.show()

            fig_fourpanels.savefig(outfolder + '/' + eventdate[j] + '_results.png')
            fig_fourpanels.clf()
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

            arri = mdates.date2num(Time.strptime(predmedian.decode(), '%Y-%m-%d %H:%M').datetime)

            arri2 = mdates.date2num(Time.strptime(predmean.decode(), '%Y-%m-%d %H:%M').datetime)

            shadelow = mdates.date2num(Time.strptime(arrival.shadelimlow.decode(), '%Y-%m-%dT%H:%M:%S.%f').datetime)
            shadehigh = mdates.date2num(Time.strptime(arrival.shadelimhigh.decode(), '%Y-%m-%dT%H:%M:%S.%f').datetime)

            # date_arrtime=mdates.num2date(arrtime)

            # figure displaying ELEvoHI prediction

            # start=min(arrtime_)
            # end=max(arrtime_)

            # dt=(end-start)/10

            # binBoundaries = []
            # binBoundaries.append(start)

            # for i in range(1,10):
            #   binBoundaries.append(binBoundaries[i-1] + dt)

            # print(binBoundaries)

            # , color=['lightblue']*len(binBoundaries)

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
            plt.text(1.15, 0.7, 'Median \n' + predmedian.decode(),
                     horizontalalignment='center', verticalalignment='center',
                     transform=ax.transAxes, color='blue')
            plt.text(1.15, 0.6, 'Mean \n' + predmean.decode(),
                     horizontalalignment='center', verticalalignment='center',
                     transform=ax.transAxes, color='green')
            plt.text(1.15, 0.45, 'Standard Deviation\n$\pm$' +
                     prederr.decode() + ' hours', horizontalalignment='center',
                     verticalalignment='center', transform=ax.transAxes,
                     color='black', bbox={'facecolor': 'slategray', 'alpha': 0.2,
                                          'pad': 10})

            if insitu != '':
                plt.text(1.15, 0.2, 'Insitu Arrival \n' + insitu,
                         horizontalalignment='center', verticalalignment='center',
                         transform=ax.transAxes, color='red')

            EarthHit = 0
            NoHit = 0
            arrtime_earth = elevohi_results.eelevohi.arrtime_earth[0]
            NrRuns = np.size(arrtime_earth)
            for i in range(0, NrRuns):
                if 'NaN' in str(arrtime_earth[i]):
                    NoHit = NoHit + 1
                else:
                    EarthHit = EarthHit + 1

            plt.text(1.15, 0.1, str(EarthHit) + '/' + str(NrRuns) +
                     ' runs hit Earth', horizontalalignment='center',
                     verticalalignment='center', transform=ax.transAxes,
                     color='black')

            plt.title('ELEvoHI Ensemble Prediction for ' + issc, fontsize=16)

            fig.subplots_adjust(hspace=0.2)

            plt.xticks(rotation=45)

            # fig.show()

            fig.savefig(outfolder + '/' + eventdate[j] + '_ELEvoHI_prediction.png',
                        bbox_inches="tight")
            fig.clf()

            plt.close('all')
            # input()


if __name__ == '__main__':
    # eventslist = ['20090623']
    main('/nas/helio/ELEvoHI_plotting/runs/')