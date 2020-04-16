import matplotlib
matplotlib.use('Agg')
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import scipy.io
import sunpy.time
import pickle
import seaborn as sns
import shutil
import datetime
from astropy.time import Time


# ###################################################### functions

# for reading catalogues
def getcat(filename):
    print('reading CAT ' + filename)
    cat = scipy.io.readsav(filename)  # , verbose='false')
    print('done reading CAT')
    return cat


def decode_array(bytearrin):
    # for decoding the strings from the IDL .sav file to a list of python
    # strings, not bytes make list of python lists with arbitrary length
    bytearrout = ['' for x in range(len(bytearrin))]
    for i in range(0, len(bytearrin) - 1):
        bytearrout[i] = bytearrin[i].decode()
    # has to be np array so to be used with numpy "where"
    bytearrout = np.array(bytearrout)
    return bytearrout


def time_to_num_cat(time_in):
    # for time conversion from catalogue .sav to numerical time
    # this for 1-minute data or lower time resolution
    # for all catalogues
    # time_in is the time in format: 2007-11-17T07:20:00 or 2007-11-17T07:20Z
    # for times help see:
    # http://docs.sunpy.org/en/latest/guide/time.html
    # http://matplotlib.org/examples/pylab_examples/date_demo2.html

    j = 0
    # time_str=np.empty(np.size(time_in),dtype='S19')
    time_str = ['' for x in range(len(time_in))]
    # =np.chararray(np.size(time_in),itemsize=19)
    time_num = np.zeros(np.size(time_in))

    for i in time_in:
        # convert from bytes (output of scipy.readsav) to string
        time_str[j] = time_in[j][0:16].decode() + ':00'
        year = int(time_str[j][0:4])
        time_str[j]
        # convert time to sunpy friendly time and to matplotlibdatetime
        # only for valid times so 9999 in year is not converted
        # pdb.set_trace()
        if year < 2100:
            time_num[j] = mdates.date2num(Time.strptime(time_str[j], '%Y-%m-%dT%H:%M:%S').datetime)
        j = j + 1

    # the date format in matplotlib is e.g. 735202.67569444
    # this is time in days since 0001-01-01 UTC, plus 1.
    # return time_num which is already an array and convert the list of strings
    # to an array
    return time_num, np.array(time_str)


def roundTime(dt=None, roundTo=60):
    # Round a datetime object to any time lapse in seconds
    # dt : datetime.datetime object, default now.
    # roundTo : Closest number of seconds to round to, default 1 minute.
    # Author: Thierry Husson 2012 - Use it as you want but don't blame me.

    if dt is None:
        dt = datetime.datetime.now()
    seconds = (dt.replace(tzinfo=None) - dt.min).seconds
    rounding = (seconds + roundTo / 2) // roundTo * roundTo
    return dt + datetime.timedelta(0, rounding - seconds, -dt.microsecond)


def plot_ellipse(ax, dayjump, pos, timeind, cmeind, k, all_apex_f, all_apex_w,
                 all_apex_r, all_apex_lon, all_apex_s, all_apex_flag,
                 frame_time_num, et_time_num_interp, et_elon_interp,
                 et_time_num, startcutFit, endcutFit):
    # Plot all the different ellipses (different runs for each time step) with
    # the elongation profile

    slope = np.NaN
    intercept = np.NaN
    # ############################## plot all active CME ellipses
    if np.size(cmeind) > 0:
        for p in range(0, np.size(cmeind)):

            # print('CME active ',p)

            # derive values for ellipse
            theta = np.arctan((all_apex_f[cmeind[0][p]] ** 2) * np.tan(
                all_apex_w[cmeind[0][p]]))
            omega = np.sqrt((np.cos(theta) ** 2) * (
                all_apex_f[cmeind[0][p]] ** 2 - 1) + 1)

            # ellipse values, depending on R and lamda and f, from Moestl et
            # al. 2015
            # Nat. Comm.
            b = (all_apex_r[cmeind[0][p]] * omega * np.sin(
                all_apex_w[cmeind[0][p]])) / (np.cos(
                    all_apex_w[cmeind[0][p]] - theta) + omega * np.sin(
                    all_apex_w[cmeind[0][p]]))
            a = b / all_apex_f[cmeind[0][p]]
            c = all_apex_r[cmeind[0][p]] - b  # center distance of ellipse

            # ellipse apex and center
            # [xapex,yapex]=np.array([np.cos(all_apex_lon[cmeind[0][p]]*np.pi/180),np.sin(all_apex_lon[cmeind[0][p]]*np.pi/180)])*all_apex_r[cmeind[0][p]]
            [xc, yc] = np.array([np.cos(
                all_apex_lon[cmeind[0][p]] * np.pi / 180),
                np.sin(all_apex_lon[cmeind[0][p]] * np.pi / 180)]) * c

            # convert only apex to show
            # now convert to polar coordinates
            # rapex=np.sqrt(xapex**2+yapex**2)
            # longapex=np.arctan2(yapex,xapex)
            # print(rapex,longapex*180/np.pi)
            # ax.scatter(longapex,rapex,c='k',s=20)

            # rc=np.sqrt(xc**2+yc**2)
            # lc=np.arctan2(yc,xc)
            # print(rc,lc*180/np.pi)
            # ax.scatter(lc,rc,c='r',s=20)
            # point at x=1 y=1
            # r1=np.sqrt(0**2+1**2)
            # l1=np.arctan2(0,1)
            # ax.scatter(l1,r1,c='b',s=50)

            # make points on ellipse
            circ_ang = ((np.arange(111) * 2 - 110) * np.pi / 180)

            xe = b * np.cos(circ_ang)   # Parameterized equation of ellipse
            ye = a * np.sin(circ_ang)

            # rotation angle
            cosang = np.cos(all_apex_lon[cmeind[0][p]] * np.pi / 180)
            # -np.deg2rad(90))
            sinang = np.sin(all_apex_lon[cmeind[0][p]] * np.pi / 180)
            # -np.deg2rad(90))

            xell = xc + xe * cosang - ye * sinang    # Rotate to desired
            # position angle
            yell = yc + xe * sinang + ye * cosang

            # now convert to polar coordinates
            rell = np.sqrt(xell ** 2 + yell ** 2)
            longell = np.arctan2(yell, xell)

            # plot in correct color
            if all_apex_s[cmeind[0][p]] == 'A':
                # make alpha dependent on distance to solar equatorial plane
                # ax.plot(longell,rell, c='grey', alpha=1-abs(all_apex_lat[
                # cmeind[0][p]]/50),
                # lw=1.5)
                ax.plot(longell, rell, c='red', alpha=0.02, lw=1.5)
            if all_apex_s[cmeind[0][p]] == 'B':
                # ax.plot(longell,rell, c='royalblue', alpha=1-abs(
                # all_apex_lat[cmeind[0][p]]/50), lw=1.5)

                # alpha should depend on colorflag
                if all_apex_flag[cmeind[0][p]] == 0:
                    # ax.plot(longell,rell, c='grey', alpha=0.6, lw=1,zorder=1)
                    ax.plot(
                        longell, rell, c='silver', alpha=0.6, lw=1, zorder=1)
                    # if all_apex_flag[cmeind[0][p]] ==1:
                    # ax.plot(longell,rell, c='silver', alpha=0.8,
                    # lw=1,zorder=2)

                if all_apex_flag[cmeind[0][p]] == 1:
                    # ax.plot(longell,rell, c='silver', alpha=0.6, lw=1,
                    # zorder=1)
                    ax.plot(
                        longell, rell, c='blue', alpha=0.02, lw=1, zorder=1)

                if all_apex_flag[cmeind[0][p]] == 2:
                    ax.plot(
                        longell, rell, c='black', alpha=1, lw=1, zorder=3)

    # ##############################plot elongation
        # difference array of current frame time frame_time_num+k to
        # position time
        # frame_time_num
        elondt = frame_time_num + k - et_time_num_interp
        # get indices where difference is less than half the time resolution
        elonind = np.where(abs(elondt) < dayjump / 2.0)

        # print( 'elonind', cmeind)

        if np.size(elonind) > 0:
            # todo: only for STB also do STA
            # check how far out tracked
            # end of fit at AU 0.4557, this is h_time_num[26]
            # = 734081.4166666657
            # before this time plot elongation as straight line
            # if frame_time_num+k <  h_time_num[26]:
            tangent_size = 1.2  # AU
            if all_apex_s[cmeind[0][p]] == 'B':
                # for ElEvoHI2 paper Amerstorfer et al. 2017
                # ######## add tangent from STEREO-B to ellipseusing the time
                # elongation profile
                # this is the currently active epsilon for the active CME
                angletox = np.deg2rad(180 - et_elon_interp[elonind[0]] - abs(
                    np.rad2deg(pos.stb[1, timeind])))  # +np.pi/2
                # make x y coordinates of tangent vector from 0/0
                vecx1 = tangent_size * np.cos(angletox)
                vecy1 = tangent_size * np.sin(angletox)
                stx = pos.stb[0, timeind] * np.cos(pos.stb[1, timeind])
                sty = pos.stb[0, timeind] * np.sin(pos.stb[1, timeind])
                elonx1 = stx + vecx1
                elony1 = sty + vecy1
                if sty > 0:
                    elony1 = sty - vecy1
                elonr = np.sqrt(elonx1 ** 2 + elony1 ** 2)
                elonlong = np.arctan2(elony1, elonx1)

                if (frame_time_num + k > et_time_num[startcutFit] and
                        frame_time_num + k < et_time_num[endcutFit]):
                    ax.plot(
                        [pos.stb[1, timeind], elonlong],
                        [pos.stb[0, timeind], elonr], c='navy', alpha=1, lw=1)
                else:
                    ax.plot(
                        [pos.stb[1, timeind], elonlong],
                        [pos.stb[0, timeind], elonr], c='navy', alpha=1,
                        lw=1, ls='--')
            if all_apex_s[cmeind[0][p]] == 'A':

                # Original
                angletox = np.deg2rad(90 - et_elon_interp[elonind[0]] - abs(
                    np.rad2deg(pos.sta[1, timeind])))

                stx = pos.sta[0, timeind] * np.cos(pos.sta[1, timeind])
                sty = pos.sta[0, timeind] * np.sin(pos.sta[1, timeind])
                vecx1 = np.sin(angletox) * tangent_size
                vecy1 = np.cos(angletox) * tangent_size
                elonx1 = stx - vecx1
                elony1 = sty - vecy1
                if sty < 0:
                    elony1 = sty + vecy1

                elonr = np.sqrt(elonx1 ** 2 + elony1 ** 2)
                elonlong = np.arctan2(elony1, elonx1)

                if (frame_time_num + k > et_time_num[startcutFit] and
                        frame_time_num + k < et_time_num[endcutFit]):
                    ax.plot(
                        [pos.sta[1, timeind], elonlong],
                        [pos.sta[0, timeind], elonr], c='darkred',
                        alpha=1, lw=1)
                else:
                    ax.plot(
                        [pos.sta[1, timeind], elonlong],
                        [pos.sta[0, timeind], elonr], c='darkred', alpha=1,
                        lw=1, ls='--')

            slope = (elony1 - sty) / (elonx1 - stx)
            intercept = sty - slope * stx

    return [slope, intercept]


def read_CME_data(read_data, dayjump, current_event_dir, ensemble_results,
                  d_days, cme_start_date_time, tracksav):
    # read all the data needed (CME parameters and elongation profiles)

    print('Start reading CME data')
    # ############ read file with ensemble results, dump as pickle to use later
    if read_data == 1:

        h = getcat(current_event_dir + ensemble_results)
        all_apex_t = h.elevo_kin.all_apex_t[0]
        startcutFit = int(h.startcut)
        endcutFit = int(h.endcut)
        [all_apex_t_num_non_interp,
            all_apex_t_num_non_interp_str] = time_to_num_cat(all_apex_t)
        # get all parameters

        all_apex_r_non_interp = h.elevo_kin.all_apex_r[0]
        all_apex_lat_non_interp = h.elevo_kin.all_apex_lat[0]  # degree
        all_apex_lon_non_interp = h.elevo_kin.all_apex_lon[0]  # degree
        # f
        all_apex_f_non_interp = h.elevo_kin.all_apex_f[0]
        # width
        all_apex_w_non_interp = np.deg2rad(h.elevo_kin.all_apex_w[0])

        # constants
        all_apex_s_non_interp = decode_array(h.elevo_kin.all_apex_s[0])
        all_apex_run_non_interp = h.elevo_kin.runnumber[0]
        all_apex_flag_non_interp = h.elevo_kin.colorflag[0]

        if d_days == 0:
            time_gt1AU = all_apex_t_num_non_interp[np.where(
                all_apex_r_non_interp > 1.5)][0]
            d_days = time_gt1AU + 0.5 - all_apex_t_num_non_interp[0]
            d_days = np.round(d_days / dayjump) * dayjump

        dur_days = d_days
        if cme_start_date_time == '':
            CME_start_time = mdates.date2num(roundTime(
                Time.strptime(all_apex_t_num_non_interp_str[0], '%Y-%m-%dT%H:%M:%S').datetime,
                roundTo=60 * 60))
            # define cme frame times
            h_time_num = np.arange(
                CME_start_time, CME_start_time + d_days, dayjump)

        # go through each run and interpolate data for each run
        # final array size -> time array of CME frames * run numbers
        finarrs = np.size(h_time_num) * np.size(
            np.unique(all_apex_run_non_interp))

        eventsize = np.size(h_time_num)

        # initialise arrays
        all_apex_t = np.zeros(finarrs)
        all_apex_r = np.zeros(finarrs)
        all_apex_lat = np.zeros(finarrs)
        all_apex_lon = np.zeros(finarrs)
        all_apex_f = np.zeros(finarrs)
        all_apex_w = np.zeros(finarrs)
        all_apex_s = [''] * finarrs
        all_apex_run = np.zeros(finarrs)
        all_apex_flag = np.zeros(finarrs)

        print('start interpolation')
    #    for q in np.arange(0, np.max(all_apex_run_non_interp)):
        for q in np.arange(0, np.size(np.unique(all_apex_run_non_interp))):

            # print(q)
            # get indices of kinematic data for this run
            thisrunind = np.where(all_apex_run_non_interp == np.unique(
                all_apex_run_non_interp)[q])

            # if there is data available for this run, interpolate to CME times
            if np.size(thisrunind) > 0:
                # these variables change with time
                # this is time, fill with frame times
                all_apex_t[eventsize * q:eventsize * (q + 1)] = h_time_num

                # fill with interpolation variables
                all_apex_r[eventsize * q:eventsize * (q + 1)] = np.interp(
                    h_time_num, all_apex_t_num_non_interp[thisrunind],
                    all_apex_r_non_interp[thisrunind])
                all_apex_lon[eventsize * q:eventsize * (q + 1)] = np.interp(
                    h_time_num, all_apex_t_num_non_interp[thisrunind],
                    all_apex_lon_non_interp[thisrunind])
                all_apex_lat[eventsize * q:eventsize * (q + 1)] = np.interp(
                    h_time_num, all_apex_t_num_non_interp[thisrunind],
                    all_apex_lat_non_interp[thisrunind])
                all_apex_f[eventsize * q:eventsize * (q + 1)] = np.interp(
                    h_time_num, all_apex_t_num_non_interp[thisrunind],
                    all_apex_f_non_interp[thisrunind])
                all_apex_w[eventsize * q:eventsize * (q + 1)] = np.interp(
                    h_time_num, all_apex_t_num_non_interp[thisrunind],
                    all_apex_w_non_interp[thisrunind])

                # fill with run numbers
                all_apex_run[eventsize * q:eventsize * (
                    q + 1)] = all_apex_run_non_interp[thisrunind][0:eventsize]

                # fill with flag numbers
                all_apex_flag[eventsize * q:eventsize * (
                    q + 1)] = all_apex_flag_non_interp[thisrunind][0:eventsize]

                # fill with observatory string
                all_apex_s[eventsize * q:eventsize * (
                    q + 1)] = all_apex_s_non_interp[thisrunind][0:eventsize]

            else:  # set all to np.nan for this run
                all_apex_t[eventsize * q:eventsize * (q + 1)] = np.nan
                all_apex_r[eventsize * q:eventsize * (q + 1)] = np.nan
                all_apex_lon[eventsize * q:eventsize * (q + 1)] = np.nan
                all_apex_lat[eventsize * q:eventsize * (q + 1)] = np.nan
                all_apex_f[eventsize * q:eventsize * (q + 1)] = np.nan
                all_apex_w[eventsize * q:eventsize * (q + 1)] = np.nan
                all_apex_run[eventsize * q:eventsize * (q + 1)] = np.nan
                all_apex_s[eventsize * q:eventsize * (q + 1)] = ''
                all_apex_flag[eventsize * q:eventsize * (q + 1)] = np.nan

        print('end interpolation')
        pickle.dump((all_apex_t, all_apex_r, all_apex_lat, all_apex_lon,
                     all_apex_f, all_apex_w, all_apex_s, all_apex_run,
                     all_apex_flag, CME_start_time, dur_days, startcutFit,
                     endcutFit),
                    open(current_event_dir + "all_apex_variables.p", "wb"))

    if read_data == 0:
        [all_apex_t, all_apex_r, all_apex_lat, all_apex_lon, all_apex_f,
            all_apex_w, all_apex_s, all_apex_run, all_apex_flag,
            CME_start_time, dur_days, startcutFit, endcutFit] = pickle.load(
                open(current_event_dir + 'all_apex_variables.p', "rb"))

        if d_days == 0:
            d_days = dur_days

        if cme_start_date_time == '':
            # CME_start_time = all_apex_t[np.isfinite(all_apex_t)][0]
            # define cme frame times
            h_time_num = np.arange(
                CME_start_time, CME_start_time + d_days, dayjump)

    # define times
    if cme_start_date_time != '':
        CME_start_time = mdates.date2num(Time.strptime(
            cme_start_date_time, '%Y-%b-%d %H:%M:%S').datetime)
        # define cme frame times
        h_time_num = np.arange(
            CME_start_time, CME_start_time + d_days, dayjump)
        # h_time_str=mdates.num2date(h_time_num)

    # ######## read and interpolate e-t profile to movie frame times - used for
    # making line from spacecraft to front
    # #et_time_num
    # #h_time_num
    # #et_elon
    # #h_et_elon
    #
    #
    # get elongation-time profile from track
    et = getcat(current_event_dir + tracksav)
    et_time = et.track.track_date[0]
    et_time_num = time_to_num_cat(et_time)[0]
    et_elon = et.track['elon'][0]

    # todo automatic interpolate values
    # linearly interpolate to hourly values make automatic later
    # et_start_time=mdates.date2num(sunpy.time.parse_time(cme_start_date_time))
    et_start_time = CME_start_time  # +dayjump
    # et_time_num_interp=np.arange(et_start_time,et_start_time+duration_days,dayjump)
    et_time_num_interp = np.arange(et_start_time, max(
        et_time_num), dayjump)
    et_elon_interp = np.interp(et_time_num_interp, et_time_num, et_elon)

    print('Finished reading CME data')
    return [CME_start_time, d_days, startcutFit, endcutFit, all_apex_t,
            all_apex_r, all_apex_lat, all_apex_lon, all_apex_f, all_apex_w,
            all_apex_s, all_apex_run, all_apex_flag, et_time_num,
            et_time_num_interp, et_elon_interp]


# ############################################################################
# ################################ main program ##############################
# ############################################################################
# ##################################### CONTROLS

# directory of current event
# Animation of ensemble simulations for ElEvoHI

# Author: C. Moestl, J. Hinterreiter, IWF Graz, Austria
# twitter @chrisoutofspace, https://github.com/cmoestl
# November 2018
# This work is published under the MIT LICENSE (see bottom)
# parameter:
#            eventsList: List with the events for which the movies should be
#                        generated
#            spaceCraft: None or 'AB' for A and B, 'A' for A and 'B'' for B
#            readData: set to 1 if you want to create the pickle file
#            coordSys: HEEQ or HEE, None for HEE
#            catPath: Path to the catalogs
#            scriptPath: Path to the ELEvoHI ensemble output
#            outPath: Path where to save the movies
#            plotBGSW: set to 'True' to plot the plot the background solar wind, if available
#            showMag: set to 'True' to plot the magnetic field legend
#            ffmpegPath: Path to ffmpeg


def main(eventsList, spaceCraft=None, readData=None, coordSys=None,
         catPath=None, scriptPath=None, outPath=None, plotBGSW=None, showMag=None, ffmpegPath=None):
    """_Bisschen Text.

    Parameters...
    eventList...

    Returns
    ...
    """

    if catPath is None:
        catPath = 'cats/'
    if showMag is None:
        showMag = False
    if scriptPath is None:
        scriptPath = 'events/'
    if spaceCraft is None or spaceCraft == 'AB':
        spacecraft = 'AB'
    if spaceCraft == 'A':
        spacecraft = 'A'
    if spaceCraft == 'B':
        spacecraft = 'B'
    plotSolarWind = True
    if plotBGSW is None or not plotBGSW:
        plotSolarWind = False
    if ffmpegPath is None:
        ffmpegPath = '/nas/helio/ffmpeg/'

    r_sun = 695700.
    au = 149597870.
    rotSun = 27.27  # days
    startBGSW = 5 # Start of the background solar wind from HUX model

    # eventsList = ['20090623']
    for l in range(0, np.size(eventsList)):
        startTime = datetime.datetime.now()
        current_event = eventsList[l]
        # current_event = '20100408'

        # IDL sav file with ensemble simulation results
        ensemble_results = 'ForMovie/formovie_all_flag.sav'

        # set 1 on the first run to produce .p save files for interpolated
        # variables needed for the movie
        read_data = 1
        if readData == 0:
            read_data = 0

        # how much time is between frames in days
        dayjump = np.double(1 / 24.0)

        # how long the movie takes in days
        # if set to zero then the duration is calculated until the CME hits
        # Earth
        duration_days = 0

        # set these values to use defined start and end time
        # set to blank to start movie 3h before and start cme 1h before
        # movie_start_date_time='2009-Jun-23 18:00:00'
        # cme_start_date_time='2009-Jun-23 20:00:00'
        movie_start_date_time = ''
        cme_start_date_time = ''

        # how long an in situ arrival stays visible in fade mode
        fadedays = 20

        # font size on bottom labels
        labelfontsize = 8

        bscale = 4

        # whether HEEQ or HEE positions are used for plot
        if coordSys is None or coordSys == 'HEE':
            HEE = 1
        if coordSys == 'HEEQ':
            HEE = 0

        coordSysString = 'HEEQ'
        if HEE == 1:
            coordSysString = 'HEE'

        # save file with elongation tracks
        # tracksav='track_B_img_ccsds.sav'
        # tracksav='satplot_track_ccsds.sav'
        # tracksav = current_event + '_ccsds.sav'

        ##########################################

        plt.close('all')

        # current_event_dir = 'events/PropDirFromElon/' + current_event
        current_event_dir = scriptPath + current_event

#        if not os.path.isdir(current_event_dir):
#            os.mkdir(current_event_dir)

#        if not os.path.isdir(current_event_dir + '/frames'):
#            os.mkdir(current_event_dir + '/frames')

        print()
        print('Start ELEvoHI animation program.')
        print()
        print('Current event ', current_event)
        print()

        if plotSolarWind:
            scDir = spaceCraft
            if scDir == 'AB':
                scDir = 'A'
            bgswFile = current_event_dir + '_' + scDir + '/bgsw_data.sav'
            if not os.path.exists(bgswFile):
                plotSolarWind = False
                print('No data for background solar wind found!')
        if plotSolarWind:
            print('Using background solar wind data')
            bgsw = getcat(bgswFile)
            bgswData = bgsw.bgsw_data
            rLen = len(bgswData)
            thetaLen = len(bgswData[1])

        # #########get ICMECAT
        filename_icmecat = catPath + 'HELCATS_ICMECAT_v10_SCEQ.sav'
        i = getcat(filename_icmecat)
        print()

        # get parameters
        # bscale makes circles larger in movie
        bmean = i.icmecat['MO_BMEAN'] * bscale
        # hee long converted to rad
        long = i.icmecat['SC_LONG_HEEQ'] * np.pi / 180
        rdist = i.icmecat['sc_heliodistance']  # AU
        sc = i.icmecat['sc_insitu']  # string
        sc = decode_array(sc)

        # get indices of events in different spacecraft
        vexind = np.where(sc == 'VEX')
        staind = np.where(sc == 'STEREO-A')
        stbind = np.where(sc == 'STEREO-B')
        winind = np.where(sc == 'Wind')
        mesind = np.where(sc == 'MESSENGER')
        ulyind = np.where(sc == 'ULYSSES')

        # make time conversion for all icme_start_time variables
        # save it as string
        icme_start_time_str = i.icmecat['icme_start_time']
        # save it as matplotlib date number
        [icme_start_time_num, icme_start_time_str] = time_to_num_cat(
            icme_start_time_str)

        # for each spacecraft, make a zeros array
        active_icme_vex = np.zeros(np.size(icme_start_time_num))
        active_icme_stb = np.zeros(np.size(icme_start_time_num))
        active_icme_sta = np.zeros(np.size(icme_start_time_num))
        active_icme_win = np.zeros(np.size(icme_start_time_num))
        active_icme_mes = np.zeros(np.size(icme_start_time_num))
        active_icme_uly = np.zeros(np.size(icme_start_time_num))

        # ####get spacecraft and planet positions
        pos = getcat(catPath + 'positions_2007_2023_' + coordSysString +
                     '_6hours.sav')
        pos_time_num = time_to_num_cat(pos.time)[0]
        # positions are available as pos.mercury etc.

        if spacecraft == 'AB' or spacecraft == 'A':
            [CME_start_time_a, duration_days_a, startcutFit_a, endcutFit_a,
             all_apex_t_a, all_apex_r_a, all_apex_lat_a, all_apex_lon_a,
             all_apex_f_a, all_apex_w_a, all_apex_s_a, all_apex_run_a,
             all_apex_flag_a, et_time_num_a, et_time_num_interp_a,
             et_elon_interp_a] = read_CME_data(
                read_data, dayjump, current_event_dir + '_A/',
                ensemble_results, duration_days, cme_start_date_time,
                current_event + '_A_ccsds.sav')

        if spacecraft == 'AB' or spacecraft == 'B':
            [CME_start_time_b, duration_days_b, startcutFit_b, endcutFit_b,
             all_apex_t_b, all_apex_r_b, all_apex_lat_b, all_apex_lon_b,
             all_apex_f_b, all_apex_w_b, all_apex_s_b, all_apex_run_b,
             all_apex_flag_b, et_time_num_b, et_time_num_interp_b,
             et_elon_interp_b] = read_CME_data(
                read_data, dayjump, current_event_dir + '_B/',
                ensemble_results, duration_days, cme_start_date_time,
                current_event + '_B_ccsds.sav')

        current_event_dir = current_event_dir + '_' + spacecraft + '/'
        if not os.path.isdir(current_event_dir):
            os.mkdir(current_event_dir)
        if not os.path.isdir(current_event_dir + 'frames'):
            os.mkdir(current_event_dir + 'frames')
        if spacecraft == 'A':
            CME_start_time = CME_start_time_a
            duration_days = duration_days_a
        if spacecraft == 'B':
            CME_start_time = CME_start_time_b
            duration_days = duration_days_b
        if spacecraft == 'AB':
            CME_start_time = min(CME_start_time_a, CME_start_time_b)
            duration_days = max(duration_days_a, duration_days_b)

        shutil.rmtree(current_event_dir + 'frames/')
        os.mkdir(current_event_dir + '/frames')

        # ################################## MAKE MOVIE FRAMES

        # initiate plot
        plt.figure(1, figsize=(8, 6), dpi=100, facecolor='w', edgecolor='w')

        sns.set_context('talk')
        sns.set_style('whitegrid')
        if not plotSolarWind:
            sns.set_style('darkgrid')

        # set start time of movie
        frame_time_num = CME_start_time - 3 / 24.
        if movie_start_date_time != '':
            frame_time_num = mdates.date2num(Time.strptime(
                movie_start_date_time, '%Y-%b-%d %H:%M:%S').datetime)

        # ##### loop over all movie frames
        for k in np.arange(0, duration_days, dayjump):
            # to current frame time, the days need to be added, so +k is done
            # save frame time as string to write on plot

            framestr = '%04i' % np.round(k * 1.0 / dayjump)
            frame_time_str = str(mdates.num2date(frame_time_num + k))

            print('frame ', framestr, '  ', frame_time_str)
            if spacecraft == 'AB' or spacecraft == 'A':
                # difference array of current frame time frame_time_num+k to
                # position time frame_time_num
                cmedt_a = frame_time_num + k - all_apex_t_a
                # get indices where difference is less than half the time
                # resolution use this to avoid nan in np.where
                cmedt_a[np.isnan(cmedt_a)] = 10000
                cmeind_a = np.where(np.abs(cmedt_a) < dayjump / 2)
                # print( 'cmeind', cmeind)

            if spacecraft == 'AB' or spacecraft == 'B':
                cmedt_b = frame_time_num + k - all_apex_t_b
                # get indices where difference is less than half the time
                # resolution use this to avoid nan in np.where
                cmedt_b[np.isnan(cmedt_b)] = 10000
                cmeind_b = np.where(np.abs(cmedt_b) < dayjump / 2)
                # print( 'cmeind', cmeind)

            # ########################################### make plot

            ax = plt.subplot(111, projection='polar')

            if plotSolarWind:
                rotAngle = (2 * np.pi / rotSun) * k

                angle = np.deg2rad(np.arange(0, 362, 362 / thetaLen)) + np.deg2rad(180) + rotAngle
                radius = np.arange(startBGSW, rLen + startBGSW) / au * r_sun
                
                thetaBGSW, rBGSW = np.meshgrid(angle, radius)
                rBGSW = np.transpose(rBGSW)
                thetaBGSW = np.transpose(thetaBGSW)
                bgswFinal = bgswData
                bgswFinal = np.flip(bgswFinal, axis=1)

                levels = np.arange(np.min(bgswFinal), np.max(bgswFinal), 1)
                cf = ax.contourf(thetaBGSW, rBGSW, bgswFinal.T, levels,
                                 cmap=plt.cm.get_cmap('rainbow'), alpha=0.2,
                                 vmin=np.min(bgswFinal), vmax=np.max(bgswFinal))

            # difference array of current frame time frame_time_num+k to
            # position time frame_time_num
            dct = frame_time_num + k - pos_time_num
            # get index of closest to 0, use this for position
            timeind = np.argmin(abs(dct))

            if spacecraft == 'AB' or spacecraft == 'A':
                [slope_a, intercept_a] = plot_ellipse(ax, dayjump, pos, timeind,
                                                      cmeind_a, k, all_apex_f_a,
                                                      all_apex_w_a, all_apex_r_a,
                                                      all_apex_lon_a, all_apex_s_a,
                                                      all_apex_flag_a, frame_time_num,
                                                      et_time_num_interp_a,
                                                      et_elon_interp_a, et_time_num_a,
                                                      startcutFit_a, endcutFit_a)
            if spacecraft == 'AB' or spacecraft == 'B':
                [slope_b, intercept_b] = plot_ellipse(ax, dayjump, pos, timeind,
                                                      cmeind_b, k, all_apex_f_b,
                                                      all_apex_w_b, all_apex_r_b,
                                                      all_apex_lon_b, all_apex_s_b,
                                                      all_apex_flag_b, frame_time_num,
                                                      et_time_num_interp_b,
                                                      et_elon_interp_b, et_time_num_b,
                                                      startcutFit_b, endcutFit_b)

            if spacecraft == 'AB':
                if (np.isfinite(slope_a) and np.isfinite(slope_b) and
                        np.isfinite(intercept_a) and np.isfinite(intercept_b)):

                    intersect_x = (intercept_b - intercept_a) / (slope_a - slope_b)
                    intersect_y = slope_a * intersect_x + intercept_a

                    elonr = np.sqrt(intersect_x ** 2 + intersect_y ** 2)
                    elonlong = np.arctan2(intersect_y, intersect_x)

                    angToA = pos.sta[1, timeind] - elonlong
                    angToB = abs(pos.stb[1, timeind]) + elonlong
                    # print('Long of apex: ', np.rad2deg(elonlong))
                    # print('Angle to A: ', np.rad2deg(angToA))
                    # print('Angle to B: ', np.rad2deg(angToB))

                    # ax.plot([0, elonlong], [0, elonr], c='green', alpha=1, lw=1)

            # index 1 is longitude, 0 is rdist
            ax.scatter(
                pos.venus[1, timeind], pos.venus[0, timeind],
                s=50, c='orange', alpha=1, lw=0, zorder=3)
            ax.scatter(
                pos.mercury[1, timeind], pos.mercury[0, timeind],
                s=50, c='dimgrey', alpha=1, lw=0, zorder=3)
            ax.scatter(
                pos.messenger[1, timeind], pos.messenger[0, timeind],
                s=25, c='dimgrey', marker='s', alpha=1, lw=0, zorder=3)
            ax.scatter(
                pos.sta[1, timeind], pos.sta[0, timeind],
                s=25, c='red', alpha=1, marker='s', lw=0, zorder=3)
            ax.scatter(
                pos.stb[1, timeind], pos.stb[0, timeind],
                s=25, c='royalblue', alpha=1, marker='s', lw=0, zorder=3)
            ax.scatter(
                pos.earth[1, timeind], pos.earth[0, timeind],
                s=50, c='mediumseagreen', alpha=1, lw=0, zorder=3)
            ax.scatter(
                pos.mars[1, timeind], pos.mars[0, timeind],
                s=50, c='orangered', alpha=1, lw=0, zorder=3)
            ax.scatter(
                pos.msl[1, timeind], pos.msl[0, timeind],
                s=25, c='magenta', marker='s', alpha=1, lw=0, zorder=3)
            ax.scatter(
                pos.maven[1, timeind], pos.maven[0, timeind],
                s=25, c='steelblue', marker='s', alpha=1, lw=0, zorder=3)
            ax.scatter(
                pos.rosetta[1, timeind], pos.rosetta[0, timeind],
                s=25, c='black', marker='s', alpha=1, lw=0, zorder=3)
            ax.scatter(
                pos.ulysses[1, timeind], pos.ulysses[0, timeind],
                s=25, c='darkolivegreen', marker='s', alpha=1, lw=0, zorder=3)

            # ######################  plot ICME detections
            # ####### for each frame time, check active ICMEs looking into ICMECAT:
            for m in range(0, len(icme_start_time_num)):

                # calculate difference in arrival_time_num_time to current frame
                icme_diff_to_frame = (frame_time_num + k) - icme_start_time_num[m]

                # for all arrival_time_num_times that are later than the current frame,
                # make them active for fadedays (fading) or infinite (keeping).
                if icme_diff_to_frame > 0 and icme_diff_to_frame < fadedays:
                    # check if this active icme belongs to a spacecraft
                    # in1d compares to arrays; true or 1 if m is contained in vexind
                    if np.in1d(m, vexind) == 1:
                        active_icme_vex[m] = icme_diff_to_frame
                    # same for the other spacecraft
                    if np.in1d(m, stbind) == 1:
                        active_icme_stb[m] = icme_diff_to_frame
                    if np.in1d(m, staind) == 1:
                        active_icme_sta[m] = icme_diff_to_frame
                    if np.in1d(m, winind) == 1:
                        active_icme_win[m] = icme_diff_to_frame
                    if np.in1d(m, mesind) == 1:
                        active_icme_mes[m] = icme_diff_to_frame
                    if np.in1d(m, ulyind) == 1:
                        active_icme_uly[m] = icme_diff_to_frame
                else:
                    # if no detection, set the index to 0
                    active_icme_vex[m] = 0
                    active_icme_stb[m] = 0
                    active_icme_sta[m] = 0
                    active_icme_win[m] = 0
                    active_icme_mes[m] = 0
                    active_icme_uly[m] = 0

            # look which ICMEs are active
            active_index_vex = np.where(active_icme_vex > 0)
            active_index_stb = np.where(active_icme_stb > 0)
            active_index_sta = np.where(active_icme_sta > 0)
            active_index_win = np.where(active_icme_win > 0)
            active_index_mes = np.where(active_icme_mes > 0)
            active_index_uly = np.where(active_icme_uly > 0)

            # fader style plot alpha dependent on time difference for this loop
            # over each element:

            for y in range(0, np.size(active_index_vex)):
                # access elements in tuple that is produced by where
                z = active_index_vex[0][y]
                # fadedays is maximum difference in time, and alpha from 0 to 1
                fadealpha = 1 - active_icme_vex[z] / (fadedays)
                ax.scatter(
                    long[z], rdist[z], s=bmean[z], c='orange',
                    alpha=fadealpha, zorder=4)

            for y in range(0, np.size(active_index_sta)):
                z = active_index_sta[0][y]
                # 30 days is maximum difference in time, and alpha from 0 to 1
                fadealpha = 1 - active_icme_sta[z] / (fadedays)
                ax.scatter(long[z], rdist[z], s=bmean[z], c='red',
                           alpha=fadealpha, zorder=4)

            for y in range(0, np.size(active_index_stb)):
                z = active_index_stb[0][y]
                # 30 days is maximum difference in time, and alpha from 0 to 1
                fadealpha = 1 - active_icme_stb[z] / (fadedays)
                ax.scatter(
                    long[z], rdist[z], s=bmean[z], c='royalblue',
                    alpha=fadealpha, zorder=4)

            for y in range(0, np.size(active_index_win)):
                z = active_index_win[0][y]
                # 30 days is maximum difference in time, and alpha from 0 to 1
                fadealpha = 1 - active_icme_win[z] / (fadedays)
                ax.scatter(
                    long[z], rdist[z], s=bmean[z], c='mediumseagreen',
                    alpha=fadealpha, zorder=4)

            for y in range(0, np.size(active_index_mes)):
                z = active_index_mes[0][y]
                # 30 days is maximum difference in time, and alpha from 0 to 1
                fadealpha = 1 - active_icme_mes[z] / (fadedays)
                ax.scatter(
                    long[z], rdist[z], s=bmean[z], c='dimgrey',
                    alpha=fadealpha, zorder=4)

            for y in range(0, np.size(active_index_uly)):
                z = active_index_uly[0][y]
                # 30 days is maximum difference in time, and alpha from 0 to 1
                fadealpha = 1 - active_icme_uly[z] / (fadedays)
                ax.scatter(
                    long[z], rdist[z], s=bmean[z], c='darkolivegreen',
                    alpha=fadealpha, zorder=4)

            # ##################### legend and additional text

            plt.suptitle('ELEvoHI ensemble simulation ')

            # Sun
            ax.scatter(0, 0, s=100, c='yellow', edgecolors='yellow')
            #plt.figtext(0.51, 0.5, 'Sun', fontsize=10, ha='center')

            # Earth
        #    plt.figtext(0.51, 0.28, 'Earth', fontsize=10, ha='center')
            plt.figtext(0.648, 0.47, 'Earth', fontsize=10, ha='center')

            plt.figtext(0.55, 0.1, coordSysString + ' longitude',
                        fontsize=10, ha='left')

            plt.figtext(0.1 - 0.02, 0.02, 'Mercury', color='dimgrey',
                        ha='center', fontsize=labelfontsize)
            plt.figtext(0.2 - 0.02, 0.02, 'MESSENGER', color='dimgrey',
                        ha='center', fontsize=labelfontsize)
            plt.figtext(0.3 - 0.02, 0.02, 'Venus', color='orange',
                        ha='center', fontsize=labelfontsize)
            plt.figtext(0.4 - 0.02, 0.02, 'STEREO-A', color='red',
                        ha='center', fontsize=labelfontsize)
            plt.figtext(0.53 - 0.02, 0.02, 'STEREO-B', color='royalblue',
                        ha='center', fontsize=labelfontsize)
            plt.figtext(0.62 - 0.02, 0.02, 'Earth', color='mediumseagreen',
                        ha='center', fontsize=labelfontsize)
            plt.figtext(0.68 - 0.02, 0.02, 'Mars', color='orangered',
                        ha='center', fontsize=labelfontsize)
            plt.figtext(0.78 - 0.02, 0.02, 'Maven', color='steelblue',
                        ha='center', fontsize=labelfontsize)
            plt.figtext(0.73 - 0.02, 0.02, 'MSL', color='magenta',
                        ha='center', fontsize=labelfontsize)
            plt.figtext(0.84 - 0.02, 0.02, 'Rosetta', color='black',
                        ha='center', fontsize=labelfontsize)
            plt.figtext(0.90 - 0.02, 0.02, 'Ulysses', color='darkolivegreen',
                        ha='center', fontsize=labelfontsize)

            
            if showMag:
                # add legend for bmean
                bleg = np.array([10, 50, 100]) * bscale
                blegstr = ['10nT', '50nT', '100nT']

                ablegr = np.zeros(len(bleg)) + [1.61, 1.6, 1.625]
                ablegt = np.radians([172.5, 180, 190])
                ax.scatter(blegt, blegr, s=bleg, c='violet', edgecolor='violet')

                for p in range(0, len(bleg)):
                    ax.annotate(
                        blegstr[p], xy=(blegt[p], blegr[p] - 0.15), ha='left',
                        va='center', fontsize=8)

            # set axes
            plt.thetagrids(
                range(0, 360, 45),
                (u'0\u00b0', u'45\u00b0', u'90\u00b0', u'135\u00b0',
                 u'\u00B1180\u00b0', u'- 135\u00b0', u'- 90\u00b0',
                 u'- 45\u00b0'),
                fmt='%d', fontsize=10)  # , frac = 1.05)
            ax.set_theta_zero_location('E')
            ax.set_ylim(0, 2.0)

            # plot text for date extra so it does not move
            # year
            # plt.figtext(0.52,0.85,frame_time_str[0:13]+':00:00',
            # fontsize=13, ha='center')
            plt.figtext(0.435, 0.85, frame_time_str[0:4], fontsize=13,
                        ha='center')
            # month
            plt.figtext(0.474, 0.85, frame_time_str[5:7], fontsize=13,
                        ha='center')
            # day
            plt.figtext(0.50, 0.85, frame_time_str[8:10], fontsize=13,
                        ha='center')
            # hours
            plt.figtext(0.54, 0.85, frame_time_str[11:13], fontsize=13,
                        ha='center')
            # mins
            plt.figtext(0.554, 0.85, ':', fontsize=13, ha='center')
            plt.figtext(0.571, 0.85, frame_time_str[14:16], fontsize=13,
                        ha='center')
            # secs
            plt.figtext(0.585, 0.85, ':', fontsize=13, ha='center')
            plt.figtext(0.601, 0.85, frame_time_str[17:19], fontsize=13,
                        ha='center')

            ax.grid(True, linestyle='--', linewidth=0.5)
            ax.set_rticks([0.5, 1, 1.5, 2])
            ax.set_rlabel_position(135)
            ax.tick_params(labelsize=labelfontsize)

            # signature
            # plt.figtext(0.95,0.01/2,r'$C. M\ddot{o}stl, T. Amerstorfer$',
            # fontsize=4, ha='center')

            if plotSolarWind:
                cax = plt.axes([0.05, 0.2, 0.02, 0.6])
                cbar = plt.colorbar(cf, cax=cax, ticks=np.arange(200, 800, 50))
                ticklabs = cbar.ax.get_yticklabels()
                cbar.ax.set_yticklabels(ticklabs, fontsize=7)
                cbar.set_label('Solar wind speed [km/s]', fontsize=7)
                
            # ##################### save frame
            plt.savefig(
                current_event_dir + '/frames/elevohi_' + framestr + '.png',
                dpi=300)
            # clears plot window
            plt.clf()

        # ########### end of loop

        if outPath is None:
            outpath = current_event_dir

        if not os.path.isdir(outpath):
            os.mkdir(outpath)

        outpath = outpath + '/'

        # convert to jpg
        # os.system('ffmpeg -i "'+current_event_dir+'frames/elevo_%04d.png" '
        # +current_event_dir+'frames/elevo_%04d.jpg -y -loglevel quiet')
        # make mp4

        
        os.system(ffmpegPath + 'ffmpeg -r 20 -i "' + current_event_dir +
                  'frames/elevohi_%04d.png" -c:v libx264 -vf "fps=25,format=yuv420p" ' + outpath +
                  current_event + '_' + spacecraft +
                  '_ensemble_movie.mp4 -y -loglevel quiet')
        plt.close('all')

        print('Made movie.')
        print('End ElEvoHI animation program.')
        print('The run took: ', datetime.datetime.now() - startTime)


if __name__ == '__main__':

    #eventslist = ['20100203', '20100319', '20100403', '20100408',
    #              '20100523', '20101026', '20110130', '20110214',
    #              '20110906', '20120123', '20120614', '20120712']
    # eventslist = ['20110906', '20120123', '20120614', '20120712']
    eventslist = ['20100523']

    main(eventslist, spaceCraft='AB', scriptPath='/nas/helio/ELEvoHI_plotting/runs/',
         catPath='/nas/helio/ELEvoHI_plotting/HI_animate_module/cats/', readData=0,
         plotBGSW=True)

# ########################## MIT license


# Copyright 2018 Mag. Dr. Christian Moestl
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

