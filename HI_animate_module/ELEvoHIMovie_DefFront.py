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
import pdb
from sympy.solvers import solve
from sympy import Symbol
import multiprocessing
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import glob
import h5py
import astropy.units as u


# ###################################################### functions

def getShockNormalAngle(pos, longell, rell, timeind, frameTime, ArrTime, plotLines):

#    print('TI: ', mdates.num2date(frameTime))
#    print('AT: ', mdates.num2date(ArrTime))

#    print('Tdiff [min]: ', TimeDiff)
    
    LonEarth = pos.earth[1, timeind]# + 0.55
#    print('LonEll: ', longell)
#    print('lonEarth: ', LonEarth)
    minDiffLonEll = min(abs(longell-LonEarth))
    indMinLon = np.where(abs(longell-LonEarth) == minDiffLonEll)[0]
    EarthHit = False
    if indMinLon < np.size(longell)-1 and indMinLon > 0:
        EarthHit = True
        
    TimeDiff = 100
#    if ArrTime != b'      -1':
    if ArrTime != float('Nan'):
        TimeDiff = abs(frameTime - ArrTime)*60*24
        
    if EarthHit and TimeDiff < 30:
        REarth = pos.earth[0, timeind]
#        plt.plot([0, LonEarth], [0, REarth], color='pink', lw=0.8, alpha=1)
 
        #if plotLines:
        #    plt.scatter(longell[indMinLon-1], rell[indMinLon-1], s=2)
        #    plt.scatter(longell[indMinLon+1], rell[indMinLon+1], s=2)

        x = rell[indMinLon]*np.cos(longell[indMinLon])
        y = rell[indMinLon]*np.sin(longell[indMinLon])
        x = REarth*np.cos(LonEarth)
        y = REarth*np.sin(LonEarth)
        x1 = rell[indMinLon-1]*np.cos(longell[indMinLon-1])
        x2 = rell[indMinLon+1]*np.cos(longell[indMinLon+1])
        y1 = rell[indMinLon-1]*np.sin(longell[indMinLon-1])
        y2 = rell[indMinLon+1]*np.sin(longell[indMinLon+1])

        k = (y1-y2)/(x1-x2)
        d = y1-k*x1

        #normale: steigung = -1/k
        fact = 1
        #if x[ind] < 0:
        #    fact = -1

        kNew = -1/k
        dNew = y-kNew*x

        dCent = 0
        kCent = y/x

        alpha = np.arctan(kCent)
 #       print('kCent [°]: ', np.rad2deg(alpha))


        # alpha = arctan(abs((m1-m2)/(1+m1*m2)))
        angleDiff = np.arctan((kNew-kCent)/(1+kNew*kCent))
        angleDiffDeg = np.rad2deg(angleDiff)

        alpha = np.arctan(kNew)
#        print('kNew [°]: ', np.rad2deg(alpha))
        dist = 0.2
        #print('x: ', x)
        #print('y: ', y)
        #print('rell: ', rell[indMinLon])
        #print('longell: ', longell[indMinLon])
        tmpXN = dist*np.cos(alpha) + x
        tmpYN = dist*np.sin(alpha) + y

        rellNew = np.sqrt(tmpXN ** 2 + tmpYN ** 2)
        longellNew = np.arctan2(tmpYN, tmpXN)
        r1 = np.sqrt(x1 ** 2 + y1 ** 2)
        l1 = np.arctan2(y1, x1)
        r2 = np.sqrt(x2 ** 2 + y2 ** 2)
        l2 = np.arctan2(y2, x2)
#        if plotLines:
#            plt.plot([LonEarth, longellNew], [REarth, rellNew], color='black', lw=0.3, alpha=1)

    #    print('angle Diff [°]= ', angleDiffDeg)
        return angleDiffDeg[0]


def plot_bgsw_speed(time, speed, angle, label, vmin, vmax, plotPath):    #arr = np.array(np.size(time_b), max(speed_b) - min(speed_b))
    ysize = np.int(max(speed) - min(speed))
    xsize = np.size(time)
    
    arr = np.zeros(shape=(xsize, ysize))
    for i in np.arange(0, xsize):
        arr[i,:] = speed[i]


    elons = np.zeros(xsize)
    for i in np.arange(0, np.size(elons)):
        elons[i] = i +1

    fig = plt.figure(figsize=(16, 5))
    ax1 = fig.add_subplot(111)
    ax1.grid(b = None, axis='both')

    #cf = ax1.imshow(arr.T, cmap=plt.cm.get_cmap('rainbow'), vmin=vmin, vmax=vmax, aspect = (xsize / ysize), origin='lower')
    cf = ax1.imshow(arr.T, cmap=plt.cm.get_cmap('coolwarm'), vmin=vmin, vmax=vmax, aspect = (xsize / ysize), origin='lower')
    #ax = plt.axes()
    plt.yticks([])
    plt.xticks(np.arange(xsize), time, rotation = 45)
    ax1.xaxis.set_major_locator(plt.MaxNLocator(np.int(xsize/8)))

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    ax2.grid(b = None, axis='both')
    ax2.set_ylabel('Elongation [°]')  # we already handled the x-label with ax1
    ax2.plot(time, np.rad2deg(angle), 'black')
    ax2.yaxis.set_ticks_position('left')
    ax2.yaxis.set_label_position('left')
    ax2.xaxis.set_major_locator(plt.MaxNLocator(np.int(xsize/8)))
    ax2.legend([label], handlelength=0, handletextpad=0, loc='upper left')
    
    
    cax = plt.axes([-0.01, 0.125, 0.02, 0.75])
    cbar = plt.colorbar(cf, cax=cax, ticks=np.arange(vmin, vmax, 50))
    cbar.set_label('Solar wind speed [km/s]')
    
    plt.savefig(plotPath + 'BGSW_' + label + '.png', dpi=300, bbox_inches='tight')
            # clears plot window
    plt.clf()


def plot_BGSW_tangent(path):
    ######################################################
    ######################################################
    # FOR a nicer plot see 'PlotAmbientSolarWinds.ipynb' #
    ######################################################
    ######################################################
    
    
    #path = 'HI_animate/events/test/20100203_AB/'
    [tpWind_a, tpWind_b, et_time_a, et_time_b, angle_a, angle_b, tp_a, tp_b] = pickle.load(
        open(path + 'tpWind_AB.p', "rb"))
    #[tpWind_a, et_time_a] = pickle.load(
    #    open('HI_animate/events/test/20100203_A/tpWind_A.p', "rb"))
    fig = plt.figure(figsize=(16, 8))

    time_a = []
    speed_a = []
    for i in np.arange(0, np.int(np.size(tpWind_a)/2)):
        #print((tpWind_a[i][0])[0:19])
        time_a.append((tpWind_a[i][0])[0:19])
        speed_a.append(tpWind_a[i][1])

    time_b = []
    speed_b = []
    for i in np.arange(0, np.int(np.size(tpWind_b)/2)):
        time_b.append((tpWind_b[i][0])[0:19])
        speed_b.append(tpWind_b[i][1])


    #x = time_a
    x = mdates.date2num(Time.strptime(time_a, '%Y-%m-%d %H:%M:%S').datetime)
    x = x - x.min()
    y = np.arange(0, len(x), 1)
    y = np.array(np.rad2deg(angle_a))

    speeds = np.array(speed_a)
    ymin = 0
    ymax = np.round(np.nanmax([np.rad2deg(angle_a), np.rad2deg(angle_b)]),-1)+10

    # Create a set of line segments so that we can color them individually
    # This creates the points as a N x 1 x 2 array so that we can stack points
    # together easily to get the segments. The segments array for line collection
    # needs to be (numlines) x (points per line) x 2 (for x and y)
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    plt.rcParams.update({'font.size': 21})
    fig, axs = plt.subplots(2, 1, sharex=True, sharey=True, figsize=[16,10])


    # Create a continuous norm to map from data points to colors
    norm = plt.Normalize(vmin, vmax)
    lc = LineCollection(segments, cmap='coolwarm', norm=norm)
    # Set the values used for colormapping
    lc.set_array(speeds)
    lc.set_linewidth(7)
    line = axs[0].add_collection(lc)
    #fig.colorbar(line, ax=axs[0])

    axs[0].set_xlim(x.min(), x.max())
    axs[0].set_ylim(ymin, ymax)
    axs[0].set_ylabel('Elongation [°]')

    #x = time_a
    x = mdates.date2num(Time.strptime(time_b, '%Y-%m-%d %H:%M:%S').datetime)
    x = x - x.min()
    y = np.array(np.rad2deg(angle_b))

    speeds = np.array(speed_b)

    # Create a set of line segments so that we can color them individually
    # This creates the points as a N x 1 x 2 array so that we can stack points
    # together easily to get the segments. The segments array for line collection
    # needs to be (numlines) x (points per line) x 2 (for x and y)
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    norm = plt.Normalize(vmin, vmax)
    lc = LineCollection(segments, cmap='coolwarm', norm=norm)
    # Set the values used for colormapping
    lc.set_array(speeds)
    lc.set_linewidth(7)
    line = axs[1].add_collection(lc)

    axs[1].set_xlim(x.min(), x.max())
    axs[1].set_ylim(ymin, ymax)

    plt.yticks(np.arange(ymin, ymax, 20.0))
    #plt.xticks(np.arange(x.min(), x.max(), 0.083))
    plt.xticks(x[0::12], time_a[0::12])
    axs[1].set_ylabel('Elongation [°]')
    plt.setp(axs[1].xaxis.get_majorticklabels(), rotation=25)

    #fig.text(0.02, 0.5, 'Elongation [°]', ha='center', va='center', rotation='vertical')

    cax = plt.axes([0.92, 0.125, 0.015, 0.755])
    cbar = plt.colorbar(line, cax=cax, ticks=np.arange(vmin, vmax, 40))
    cbar.set_label('Solar wind speed [km/s]')

    axs[0].text(0.2, ymax-5, 'a)', fontsize=28, ha='center', va='top', wrap=True)
    axs[1].text(0.2, ymax-5, 'b)', fontsize=28, ha='center', va='top', wrap=True)

    fig.savefig(path + '/BGSW_elon.png',
                bbox_inches="tight")
    fig.clf()
    plt.close('all')

    print('done')
    

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

def getTangentPoint(a, b, xc, yc, px, py, elon, sc, plot):
    tilt = 90

    pxOri = px
    pyOri = py
    
    px = px - xc
    py = py - yc

    ti = np.deg2rad(elon)

    pxRot = px*np.cos(ti) - py*np.sin(ti)
    pyRot = px*np.sin(ti) + py*np.cos(ti)

    px = pxRot
    py = pyRot

    ellipseResolution = 211
    circ_ang = ((np.arange(ellipseResolution) * 2 - (ellipseResolution-1)) * np.pi / 180)

    xe = b * np.cos(circ_ang)   # Parameterized equation of ellipse
    ye = a * np.sin(circ_ang)

    cosang = np.cos(tilt * np.pi / 180)
    sinang = np.sin(tilt * np.pi / 180)


    xell = xe * cosang - ye * sinang    # Rotate to desired
                # position angle
    yell = xe * sinang + ye * cosang


    if py != 0:
        xSolve = Symbol('xSolve')
        xSol = solve(b**2*xSolve**2 + a**2*((a**2*b**2-b**2*xSolve*px)/(a**2*py))**2-a**2*b**2, xSolve)

        #print(xSol)
        xs = []
        for xst in xSol:
            xs.append(float(xst))
        #print(xs)
        xs =[np.max(xs)]
        ys = []
        ytmp = Symbol('ytmp')
        for xtmp in xs:
            tmp = solve((b**2*xtmp**2 + a**2*ytmp**2 - a**2*b**2))
            ys.append(tmp)

        if sc == 'A':
            if np.max(xell) < px:
                ys = np.min(ys)
            else:
                ys = np.max(ys)
        if sc == 'B':
            if np.max(xell) < px:
                ys = np.max(ys)
            else:
                ys = np.min(ys)

    xt1 = 0
    xt2 = 0
    yt1 = 0
    yt2 = 0
    if plot == 1:
        #d = (py - k * px)
        k = Symbol('k')
        kSol = solve(a**2*k**2+b**2-(py-k*px)**2, k)
        #print('kSol =', kSol)
        k1 = float(kSol[0])
        d1 = (py - k1*px)

        k2 = float(kSol[1])
        d2 = (py - k2*px)

        #y = (k*x + d)
        xtest = np.arange(-1, 1, 0.005)


        plt.gca().set_aspect('equal', adjustable='box')
        plt.xlim(-1,1)
        plt.ylim(-1,1)
        plt.plot(xell, yell)
        plt.scatter(px, py, s=10, c='red')
        for i in np.arange(0, np.size(xtest)):
            plt.scatter(xtest[i], k1*xtest[i] + d1, s=1, c='red')
            plt.scatter(xtest[i], k2*xtest[i] + d2, s=1, c='red')
        for i in np.arange(0, np.size(xs)):
            plt.scatter(xs[i], ys, s=10, c='black')



        plt.xlim(-1,1)
        plt.ylim(-1,1)
        
        #ytest2 = k2*xtest + d2
        #ytest1 = k1*xtest + d1

        #xt1 = xtest*np.cos(ti) - ytest1*np.sin(ti) + xc
        #yt1 = xtest*np.sin(ti) + ytest1*np.cos(ti) + yc
        #xt2 = xtest*np.cos(ti) - ytest2*np.sin(ti) + xc
        #yt2 = xtest*np.sin(ti) + ytest2*np.cos(ti) + yc
        
    points = []
    ti = -ti
    for i in np.arange(0, np.size(xs)):
        xRot = xs[i]*np.cos(ti)-ys*np.sin(ti) + xc
        yRot = xs[i]*np.sin(ti)+ys*np.cos(ti) + yc
        points.append([xRot, yRot])
  
    a = np.array(np.float64(points[0]))
    b = np.array([pxOri, pyOri])
    c = np.array([0,0])
 
    ba = a - b
    bc = c - b

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)

    return points, xt1, yt1, xt2, yt2, angle


def plot_ellipse(ax, dayjump, pos, timeind, cmeind, k, all_apex_f, all_apex_w,
                 all_apex_r, all_apex_lon, all_apex_s, all_apex_flag,
                 frame_time_num, et_time_num_interp, et_elon_interp,
                 et_time_num, startcutFit, endcutFit, tangentPoints, frontData, ELEvoHIResults, shockAngles, diffMeasure):
    # Plot all the different ellipses (different runs for each time step) with
    # the elongation profile

    r_sun = 695700.
    au = 149597870.
    slope = np.NaN
    intercept = np.NaN
    angle = np.NaN
    stopAnimation = False
    # ############################## plot all active CME ellipses
    if np.size(cmeind) > 0:
#        print(np.size(cmeind))
        
        indizes = np.arange(0, np.size(cmeind))
        indizes = np.roll(indizes, -1)

        for p in indizes:
#        for p in range(0, np.size(cmeind)):
#            if p == 0:
            if True:
                # print('CME active ',p)

#                print('f: ', all_apex_f[cmeind[0][p]])
#                print('lon: ', all_apex_lon[cmeind[0][p]])
#                print('w: ', np.rad2deg(all_apex_w[cmeind[0][p]]))
#                print('lon: ', all_apex_s[cmeind[0][p]], all_apex_lon[cmeind[0][p]])
                
                # derive values for ellipse
        
                phi = all_apex_lon[cmeind[0][p]]
                lamb = np.rad2deg(all_apex_w[cmeind[0][p]])
                f = all_apex_f[cmeind[0][p]]
                
#                print('phi', phi)
#                print('lamb', lamb)
#                print('f', f)

               
                LonE = pos.earth[1, timeind]
                LonST = pos.sta[1, timeind]
                if all_apex_s[cmeind[0][p]] == 'B':
                    LonST = pos.stb[1, timeind]
                
  #              print('diff: ', np.rad2deg(LonST - LonE) - phi, lamb, f)
  #              print(ELEvoHIResults[:,p])                
                
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

               
                #ellipseResolution = 211
#                circ_ang = ((np.arange(111) * 2 - 110) * np.pi / 180)
                circ_ang = ((np.arange(101) * 2 - 100) * np.pi / 180)
#                circ_ang1 = ((np.arange(0, 201, 0.5) * 2 - 100) * np.pi / 180)
#                circ_ang1 = circ_ang1[0:201]
#                circ_ang = circ_ang1
                
                #circ_ang = ((np.arange(ellipseResolution) * 2 - (ellipseResolution-1)) * np.pi / 180
                #circ_ang = (np.arange(-180, 180.25, 0.25) * np.pi / 180)
            

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
                #xell = xc + xe
                #yell = yc + ye
                
                rell = np.sqrt(xell ** 2 + yell ** 2)
                longell = np.arctan2(yell, xell)
                
#                print('r0: ', all_apex_r[cmeind[0][p]])
#                print(rell)
                                
                frameTime = frame_time_num + k
                
                #print('f: ', all_apex_f[cmeind[0][p]])
                #print('min lon: ', np.min(longell))
                #print('max lon: ', np.max(longell))

        
                # plot in correct color
                if all_apex_s[cmeind[0][p]] == 'A':
                    # make alpha dependent on distance to solar equatorial plane
                    # ax.plot(longell,rell, c='grey', alpha=1-abs(all_apex_lat[
                    # cmeind[0][p]]/50),
                    # lw=1.5)
                    #if True:
                    #    ax.plot(longell, rell, c='red', alpha=1, lw=0.8)
                    #else:
                    #    ax.plot(longell, rell, c='red', alpha=0.08, lw=1.5)
                        
                    stx = pos.sta[0, timeind] * np.cos(pos.sta[1, timeind])
                    sty = pos.sta[0, timeind] * np.sin(pos.sta[1, timeind])
                    elon = all_apex_lon[cmeind[0][p]]
                    
#                    pos.earth[1, timeind], pos.earth[0, timeind]

#                    shockAngle = getShockNormalAngle(pos, longell, rell, timeind, frameTime, ELEvoHIResults[:,p][3], False)
#                    print('SA: ', shockAngle)
#                    if shockAngle != None:
#                        shockAngles.append(shockAngle)
            
            
#                    frontKins = frontData.frontkins[p]

#                    frontTimeNum = []
#                    for tn in frontKins.timearr:
#                        frontTimeNum.append(mdates.date2num(Time.strptime(tn, '%Y-%m-%dT%H:%M:%S.%f').datetime))

#                    frontTime = np.array(frontTimeNum)
#                    ftArrTime = mdates.date2num(Time.strptime(frontKins.arrtimeearth, '%Y-%m-%dT%H:%M:%S.%f').datetime)

                    frontKins = frontData[p]
                    frontTime = frontKins.timearr
                    ftArrTime = frontKins.arrtimeearth
                    
                    #print('ft: ', mdates.num2date(frontTime))
                    indFront = np.where(abs(frontTime - frameTime) == min(abs(frontTime - frameTime)))[0]
                    indexFront = 0
                    if len(indFront) > 0:
                        indexFront = indFront[0]
                        
                    # get a measure for how deformed the CME front is
                    rEllipse = rell
                    longEllipse = longell 
                    if indexFront != 0:   
                        rell = frontKins.frontarr[:,indexFront]*r_sun/au
                        longell = np.deg2rad(frontKins.longitude*-1)
                        
                        apexR = rell[int(len(rell)/2)]
                        if (apexR > 1.0):
                            dM = np.array(diffMeasure)
                            rDiff = np.abs(rell - rEllipse)
                            if len(dM) == 0:
                                diffMeasure.append([p, np.mean(rDiff), np.std(rDiff)])
                            elif p not in dM[:,0]:
                                diffMeasure.append([p, np.mean(rDiff), np.std(rDiff)])

                                                
                        shockAngle = getShockNormalAngle(pos, longell, rell, timeind, frameTime, ftArrTime, True)
                        if shockAngle != None:
                            print('DefFront ShockNormal A: ', shockAngle)
                            shockAngles.append(shockAngle)
                        
                        if p == 0:
                            ax.plot(longell, rell, c='darkred', alpha=1, lw=0.8)
                        else:
                            ax.plot(longell, rell, c='red', alpha=0.08, lw=0.4)

                            
                    if indexFront == 0:
                        ax.plot(longEllipse, rEllipse, c='red', alpha=0.08, lw=0.4)
                    if p == 0:
                        ax.plot(longEllipse, rEllipse, c='forestgreen', alpha=1, lw=0.8)




                    #if indexFront > 0 and rell[0] == 0:
                        #stopAnimation = True
                        #print('stop animation A')


                    if False: # remove, use next line
                    #if p == 0:
                        points, xtest1, ytest1, xtest2, ytest2, angle = getTangentPoint(b, a, xc, yc, stx, sty, -elon, all_apex_s[cmeind[0][p]], 0)

                        tpR = np.sqrt(float(points[0][0]) ** 2 + float(points[0][1]) ** 2)
                        tpLon = np.arctan2(float(points[0][1]), float(points[0][0]))

                        ax.plot([pos.sta[1, timeind], tpLon], [pos.sta[0, timeind], tpR], c='black', linewidth=0.2, alpha = 1)

                        tangentPoints.append([tpR, tpLon])
                        #ax.scatter(np.deg2rad(all_apex_lon[cmeind[0][p]]), all_apex_r[cmeind[0][p]], s=0.7, c='green')
                        for pt in tangentPoints:
                            tpR = pt[0]
                            tpLon = pt[1]
                            ax.scatter(tpLon, tpR, s=0.5, c='red')
                        

                if all_apex_s[cmeind[0][p]] == 'B':
                    # ax.plot(longell,rell, c='royalblue', alpha=1-abs(
                    # all_apex_lat[cmeind[0][p]]/50), lw=1.5)

                        
                    frontKins = frontData[p]
                    frontTime = frontKins.timearr
                    ftArrTime = frontKins.arrtimeearth
                    frameTime = frame_time_num + k
                    #print('ft: ', mdates.num2date(frontTime))
                    indFront = np.where(abs(frontTime - frameTime) == min(abs(frontTime - frameTime)))[0]
                    indexFront = 0
                    if len(indFront) > 0:
                        indexFront = indFront[0]
                     
                    #print('FT: ', mdates.num2date(frontTime[indexFront]))
                    
                    rEllipse = rell
                    longEllipse = longell 
                    if indexFront != 0:
                        rEllipse = rell
                        longEllipse = longell                            
                            
                        rell = frontKins.frontarr[:,indexFront]*r_sun/au
                        longell = np.deg2rad(frontKins.longitude*-1)
                        
                        apexR = rell[int(len(rell)/2)]
                        if (apexR > 1.0):
                            dM = np.array(diffMeasure)
                            rDiff = np.abs(rell - rEllipse)
                            if len(dM) == 0:
                                diffMeasure.append([p, np.mean(rDiff), np.std(rDiff)])
                            elif p not in dM[:,0]:
                                diffMeasure.append([p, np.mean(rDiff), np.std(rDiff)])
                                
                        shockAngle = getShockNormalAngle(pos, longell, rell, timeind, frameTime, ftArrTime, True)
                        if shockAngle != None:
                            print('DefFront ShockNormal B: ', shockAngle)
                            shockAngles.append(shockAngle)
                        
                        if p == 0:
                            ax.plot(longell, rell, c='navy', alpha=1, lw=0.8)
                        else:
                            ax.plot(longell, rell, c='blue', alpha=0.08, lw=0.4)

                            
                    if indexFront == 0:
                        ax.plot(longEllipse, rEllipse, c='blue', alpha=0.08, lw=0.4)
                    if p == 0:
                        ax.plot(longEllipse, rEllipse, c='forestgreen', alpha=1, lw=0.8)

                    #if indexFront > 0 and rell[0] == 0:
                        #stopAnimation = True
                        #print('Stop animation B')
                        #print('index front: ', indexFront)
                        #print('rell[0]: ', rell[0])
                        #print('p: ', p)
                        #print('FrameTime: ', mdates.num2date(frameTime))
                        
                    if False:
                        if p == 0:
                            stx = pos.stb[0, timeind] * np.cos(pos.stb[1, timeind])
                            sty = pos.stb[0, timeind] * np.sin(pos.stb[1, timeind])
                            elon = all_apex_lon[cmeind[0][p]]
                            points, xtest1, ytest1, xtest2, ytest2, angle = getTangentPoint(b, a, xc, yc, stx, sty, -elon, all_apex_s[cmeind[0][p]], 0)

                            tpR = np.sqrt(float(points[0][0]) ** 2 + float(points[0][1]) ** 2)
                            tpLon = np.arctan2(float(points[0][1]), float(points[0][0]))

                            ax.plot([pos.stb[1, timeind], tpLon], [pos.stb[0, timeind], tpR], c='black', linewidth=0.2, alpha = 1)

                            tangentPoints.append([tpR, tpLon])
                            #ax.scatter(np.deg2rad(all_apex_lon[cmeind[0][p]]), all_apex_r[cmeind[0][p]], s=0.7, c='green')
                            for pt in tangentPoints:
                                tpR = pt[0]
                                tpLon = pt[1]
                                ax.scatter(tpLon, tpR, s=0.5, c='blue')
                                
    # ##############################plot elongation
        # difference array of current frame time frame_time_num+k to
        # position time
        # frame_time_num
        elondt = frame_time_num + k - et_time_num_interp
        # get indices where difference is less than half the time resolution
        elonind = np.where(abs(elondt) < dayjump / 2.0)

        # print( 'elonind', cmeind)

        if np.size(elonind) > 0:
            tangent_size = 1.2  # AU
            tangent_size = np.arange(0, 1.2, 0.0025)
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
                        [pos.stb[1, timeind], elonlong[-1]],
                        [pos.stb[0, timeind], elonr[-1]], c='navy', alpha=1, lw=0.5)
                else:
                    ax.plot(
                        [pos.stb[1, timeind], elonlong[-1]],
                        [pos.stb[0, timeind], elonr[-1]], c='navy', alpha=1,
                        lw=0.5, ls='--')

                    
                tangentPoint = get_tangentPoint(longell, rell, elonlong, elonr)
                tangentPointCart = get_tangentPointCart(xell, yell, elonx1, elony1)
                tpR = np.sqrt(tangentPointCart[0][0] ** 2 + tangentPointCart[0][1] ** 2)
                tpLon = np.arctan2(tangentPointCart[0][1], tangentPointCart[0][0])
                

                #ax.scatter(tangentPoint[0][0], tangentPoint[0][1], s=0.5, c='black')
                #ax.scatter(tangentPoint[1][0], tangentPoint[1][1], s=0.5, c='blue')
                #ax.scatter(tpLon, tpR, s=0.5, c='red')
                
                
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
                        [pos.sta[1, timeind], elonlong[-1]],
                        [pos.sta[0, timeind], elonr[-1]], c='darkred',
                        alpha=1, lw=0.5)
                else:
                    ax.plot(
                        [pos.sta[1, timeind], elonlong[-1]],
                        [pos.sta[0, timeind], elonr[-1]], c='darkred', alpha=1,
                        lw=0.5, ls='--')
            
            if False:
                slope = (elony1 - sty) / (elonx1 - stx)
                intercept = sty - slope * stx
                slope = slope[-1]
                intercept = intercept[-1]

    return [slope, intercept, tangentPoints, angle, stopAnimation, shockAngles, diffMeasure]


def read_CME_data(read_data, dayjump, current_event_dir, ensemble_results,
                  d_days, cme_start_date_time, tracksav):
    # read all the data needed (CME parameters and elongation profiles)

    print('Start reading CME data')
    # ############ read file with ensemble results, dump as pickle to use later
    if read_data == 1:
        print('start transforming front data')
        frontFile = current_event_dir + '/results/frontDataAll.sav'
        frontData = getcat(frontFile)
        for i in np.arange(0, len(frontData.frontkins.timearr)):
            for j in np.arange(0, len(frontData.frontkins.timearr[i])):
                if frontData.frontkins.timearr[i][j] != b'      -1':
                    frontData.frontkins.timearr[i][j] = mdates.date2num(
                        Time.strptime(frontData.frontkins.timearr[i][j], '%Y-%m-%dT%H:%M:%S.%f').datetime)
                else:
                    frontData.frontkins.timearr[i][j] = float('Nan')
        for i in np.arange(0, len(frontData.frontkins.arrtimeearth)):
            if frontData.frontkins.arrtimeearth[i] != b'      -1':
                frontData.frontkins.arrtimeearth[i] = mdates.date2num(
                Time.strptime(frontData.frontkins.arrtimeearth[i], '%Y-%m-%dT%H:%M:%S.%f').datetime)
            else:
                frontData.frontkins.arrtimeearth[i] = float('Nan')
        

        frontDataKins = frontData.frontkins
        
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
                all_apex_r_non_interp > 2.0)][0]
            d_days = time_gt1AU + 1.0 - all_apex_t_num_non_interp[0]
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
                     endcutFit, frontDataKins), open(current_event_dir + "all_apex_variables.p", "wb"))

    if read_data == 0:
        [all_apex_t, all_apex_r, all_apex_lat, all_apex_lon, all_apex_f,
            all_apex_w, all_apex_s, all_apex_run, all_apex_flag,
            CME_start_time, dur_days, startcutFit, endcutFit, frontDataKins] = pickle.load(
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
            et_time_num_interp, et_elon_interp, frontDataKins]


def get_tangentPointCart(ell_x, ell_y, tangent_x, tangent_y):
    minDist = 10000
    indEll = -1
    indTang = -1
    for i in np.arange(0, np.size(ell_x)):
        for j in np.arange(0, np.size(tangent_x)):
            dist = (ell_x[i] - tangent_x[j]) ** 2 + (ell_y[i] - tangent_y[j]) ** 2
            if dist < minDist:
                minDist = dist
                indEll = i
                indTang = j
    #print(minDist)
    #print(ell_elon[indEll], ell_r[indEll])
    #print(tangent_elon[indTang], tangent_r[indTang])
    return [[ell_x[indEll], ell_y[indEll]], [tangent_x[indTang], tangent_y[indTang]], indEll, indTang]



def get_tangentPoint(ell_elon, ell_r, tangent_elon, tangent_r):
    minDist = 100
    indEll = -1
    indTang = -1
    for i in np.arange(0, np.size(ell_elon)):
        for j in np.arange(0, np.size(tangent_elon)):
            dist = tangent_r[j]**2 + ell_r[i]**2 - 2*tangent_r[j]*ell_r[i]*np.cos(tangent_elon[j] - ell_elon[i])
            if dist < minDist:
                minDist = dist
                indEll = i
                indTang = j
    #print(minDist)
    #print(ell_elon[indEll], ell_r[indEll])
    #print(tangent_elon[indTang], tangent_r[indTang])
    return [[ell_elon[indEll], ell_r[indEll]], [tangent_elon[indTang], tangent_r[indTang]], indEll, indTang]


# ############################################################################
# ################################ main program ##############################
# ############################################################################
# ##################################### CONTROLS

# directory of current event
# Animation of ensemble simulations for ElEvoHI

# Author: C. Moestl, IWF Graz, Austria
# twitter @chrisoutofspace, https://github.com/cmoestl
# November 2018
# This work is published under the MIT LICENSE (see bottom)
# parameter:
#            eventsList: List with the events for which the movies should be
#                        generated
#            spaceCraft: None for A and B, A for A and B for B
#            readData: set to 1 if you want to create the pickle file
#            coordSys: HEEQ or HEE, None for HEE
#            catPath: Path to the catalogs
#            scriptPath: Path to the ELEvoHI ensemble output
#            outPath: Path where to save the movies


def main(eventsList, spaceCraft=None, readData=None, coordSys=None,
         catPath=None, scriptPath=None, outPath=None, plotBGSW=None, showMag=None, bgswModel=None, huxtPath=None, ffmpegPath=None):
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
    if bgswModel is None or not bgswModel:
        plotSolarWind = False
    if huxtPath is None:
        huxtPath = '/nas/helio/data/bgsw_HUXt/'
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
        #movie_start_date_time='2012-Jun-14 18:00:00'
        cme_start_date_time = ''

        # how long an in situ arrival stays visible in fade mode
        fadedays = 20

        # font size on bottom labels
        labelfontsize = 13

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
            if bgswModel == 'HUX':
                scDir = spaceCraft
                if scDir == 'AB':
                    scDir = 'A'
                bgswFile = current_event_dir + '_' + scDir + '/results/bgsw_Data.sav'
                bgsw = getcat(bgswFile)
                bgswDat = bgsw.bgsw_data
                bgswStartTime = (bgsw.bgswStartTime[0]).decode()
                bgswStartTime_num = mdates.date2num(Time.strptime(bgswStartTime, '%Y-%m-%d %H:%M').datetime)
                bgswData = bgswDat.copy()
                bgswDataNew = bgswDat.copy()
    #            bgswData[208:212, 30] = 600
                rLen = len(bgswData)
                thetaLen = len(bgswData[1])
            
            if bgswModel == 'HUXt':
                data_path = huxtPath + current_event + '.hdf5'
                data = h5py.File(data_path, 'r')
                time = data['time_out'][()] # model time in seconds
                r = data['r'][()]
                r_grid = data['r_grid'][()]
                lon = data['lon'][()]
                lon_grid = data['lon_grid'][()]
                dlon = data['dlon'][()] * u.Unit(data['dlon'].attrs['unit'])
                v_grid = data['v_grid'][()]
                print(np.min(v_grid), np.max(v_grid))
                angle = lon
                radius = r*r_sun/au
                fmt = '%Y-%m-%d %H:%M'
#                startDate = '2010-02-03 22:00'
                f = open(huxtPath + current_event + '_tinit.txt', 'r')
                startDate = f.read()
                startDate = mdates.date2num(Time.strptime(startDate, fmt).datetime)
                bgswTime = startDate + time/24/60/60

                bgswStartTime_num = bgswTime[0]
                
            if bgswModel == 'EUHFORIA':
                scDir = spaceCraft
                if scDir == 'AB':
                    scDir = 'A'
                bgswFile = current_event_dir + '_' + scDir + '/results/bgswData_EUHFORIA.sav'
                bgsw = getcat(bgswFile)
                bgsw = bgsw.bgswData
                bgswData = bgsw.varr[0]
                bgswStartTime = bgsw.tinitstr[0]
                bgswStartTime_num = mdates.date2num(Time.strptime(bgswStartTime, '%Y-%m-%d %H:%M').datetime)
                bgswDataNew = bgswData.copy()
    #            bgswData[208:212, 30] = 600
                bgswRadii = bgsw.r[0]
                bgswLons = bgsw.lon[0]
            

        # #########get ICMECAT
        filename_icmecat = catPath + 'HELCATS_ICMECAT_v10_SCEQ.sav'
        i = getcat(filename_icmecat)

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

        if True:
            print('read_data = ', read_data)
            print('dayjump = ', dayjump)
            print('current_event_dir = ', current_event_dir)
            print('ensemble_results = ', ensemble_results)
            print('cme_start_date_time = ', cme_start_date_time)
            print('current_event = ', current_event)
            
        if spacecraft == 'AB' or spacecraft == 'A':
            [CME_start_time_a, duration_days_a, startcutFit_a, endcutFit_a,
             all_apex_t_a, all_apex_r_a, all_apex_lat_a, all_apex_lon_a,
             all_apex_f_a, all_apex_w_a, all_apex_s_a, all_apex_run_a,
             all_apex_flag_a, et_time_num_a, et_time_num_interp_a,
             et_elon_interp_a, frontData_a] = read_CME_data(
                read_data, dayjump, current_event_dir + '_A/',
                ensemble_results, duration_days, cme_start_date_time,
                current_event + '_A_ccsds.sav')
            
            CME_start_time = CME_start_time_a
#            duration_days = duration_days_a
            duration_days_a = max(frontData_a.arrtimeearth)-CME_start_time + 1.2
            duration_days = duration_days_a

        if spacecraft == 'AB' or spacecraft == 'B':
            [CME_start_time_b, duration_days_b, startcutFit_b, endcutFit_b,
             all_apex_t_b, all_apex_r_b, all_apex_lat_b, all_apex_lon_b,
             all_apex_f_b, all_apex_w_b, all_apex_s_b, all_apex_run_b,
             all_apex_flag_b, et_time_num_b, et_time_num_interp_b,
             et_elon_interp_b, frontData_b] = read_CME_data(
                read_data, dayjump, current_event_dir + '_B/',
                ensemble_results, duration_days, cme_start_date_time,
                current_event + '_B_ccsds.sav')
            
            CME_start_time = CME_start_time_b
            duration_days_b = max(frontData_b.arrtimeearth)-CME_start_time + 1.2
            duration_days = duration_days_b

        current_event_dir = current_event_dir + '_' + spacecraft + '/'
        if not os.path.isdir(current_event_dir):
            os.mkdir(current_event_dir)
        if not os.path.isdir(current_event_dir + 'frames'):
            os.mkdir(current_event_dir + 'frames')

        if spacecraft == 'AB':
            CME_start_time = min(CME_start_time_a, CME_start_time_b)
            duration_days = max(duration_days_a, duration_days_b) + abs(CME_start_time_a-CME_start_time_b)
            
        print('duration_days = ', duration_days)
        print('StartTime: ', mdates.num2date(CME_start_time))
        print('EndTime: ', mdates.num2date(CME_start_time + duration_days))
        
        ####################
        ### front data:
        ####################
#        print('CED: ', current_event_dir)
#        for file in glob.glob(current_event_dir + '/results/frontDataAll.sav'):
#            frontFile = file
            
#        front=frontFile
#        frontData = getcat(front)

        
        #print(frontTimeNum)

        shutil.rmtree(current_event_dir + 'frames/')
        os.mkdir(current_event_dir + '/frames/')
        
        ####################
        ### ELEvoHI restuls:
        ####################
        ELEvoHIResults = None
#        Resfile = 'eELEvoHI_results.sav'
#        h = getcat(current_event_dir + Resfile)
#        arrTimeNum = []
#        for at in h.eelevohi.arrtime_earth[0]:
#            if at.decode() != 'NaN':
#                arrTimeNum.append(mdates.date2num(Time.strptime(at.decode(), '%Y-%m-%dT%H:%M').datetime))
#            else:
#                arrTimeNum.append(float('Nan'))
#        ELEvoHIResults = np.zeros((4,np.size(h.eelevohi.arrtime_earth[0])))
#        ELEvoHIResults[0] = h.eelevohi.phi[0]
#        ELEvoHIResults[1] = h.eelevohi.LAMBDA[0]
#        ELEvoHIResults[2] = h.eelevohi.f[0]
#        ELEvoHIResults[3] = arrTimeNum


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

        tangentPoints_a = []
        tangentPoints_b = []
        tangAngle_a = []
        tangAngle_b = []
        tpWind_a = []
        tpWind_b = []
        shockAngles_a = []
        shockAngles_b = []
        diffMeasure_a = []
        diffMeasure_b = []
        rn = 0
        # ##### loop over all movie frames
        for k in np.arange(0, duration_days, dayjump):
#        for k in np.arange(3.3, duration_days, dayjump):
#        for k in np.arange(-0.587, duration_days, dayjump):
#        for k in np.arange(0.0793, duration_days, dayjump):
            # to current frame time, the days need to be added, so +k is done
            # save frame time as string to write on plot
            rn = rn+1

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
                bgswMin = 250
                bgswMax = 700
                if bgswModel == 'HUX':
                    diffFrameTime =  (frame_time_num - bgswStartTime_num) * (2 * np.pi / rotSun)
                    rotAngle = (2 * np.pi / rotSun * k)# + (2 * np.pi / rotSun)*diffFrameTime
                    rotAngle = np.deg2rad(60) + rotAngle + diffFrameTime# 60° comes from the BGSW output
                    angle = np.deg2rad(np.arange(0, 362, 362 / thetaLen)) + rotAngle
                    radius = np.arange(startBGSW, rLen + startBGSW) / au * r_sun
                    thetaBGSW, rBGSW = np.meshgrid(angle, radius)
                    rBGSW = np.transpose(rBGSW)
                    thetaBGSW = np.transpose(thetaBGSW)
                    bgswFinal = bgswData
                    bgswFinal = np.flip(bgswFinal, axis=1)
                    levels = np.arange(bgswMin, bgswMax, 1)
                    cf = ax.contourf(thetaBGSW, rBGSW, bgswFinal.T, levels,
                                     cmap=plt.cm.get_cmap('coolwarm'), alpha=0.2,
                                     vmin=bgswMin, vmax=bgswMax)


                    
                if bgswModel == 'HUXt':
                    ft = frame_time_num + k
                    indMin = np.argmin(np.abs(bgswTime-ft))
                    v_sub = v_grid[indMin]
                    bgswData = v_sub

                    pad = lon_grid[:, 0].reshape((lon_grid.shape[0], 1)) + 2*np.pi
                    lonTmp = np.concatenate((lon_grid, pad), axis=1)
                    pad = r_grid[:, 0].reshape((r_grid.shape[0], 1))
                    radTmp = np.concatenate((r_grid, pad), axis=1)*r_sun/au
                    pad = v_sub[:, 0].reshape((v_sub.shape[0], 1))
                    v = np.concatenate((v_sub, pad), axis=1)
                    rLen = len(bgswData)
                    thetaLen = len(bgswData[1])
                    levels = np.arange(bgswMin, bgswMax, 1)
                    cf = ax.contourf(lonTmp, radTmp, v, levels=levels, cmap=plt.cm.get_cmap('coolwarm'), extend='both', alpha=0.2)
            
                if bgswModel == 'EUHFORIA':
                    diffFrameTime =  (frame_time_num - bgswStartTime_num) * (2 * np.pi / rotSun)
                    rotAngle = (2 * np.pi / rotSun * k)# + (2 * np.pi / rotSun)*diffFrameTime
                    rotAngle = rotAngle + diffFrameTime# 60° comes from the BGSW output
                    angle = np.deg2rad(bgswLons) + rotAngle
                    radius = bgswRadii*r_sun/au
                    thetaBGSW, rBGSW = np.meshgrid(angle, radius)
                    rBGSW = np.transpose(rBGSW)
                    thetaBGSW = np.transpose(thetaBGSW)
                    bgswFinal = bgswData
                    levels = np.arange(bgswMin, bgswMax, 1)
                    cf = ax.contourf(thetaBGSW, rBGSW, bgswFinal.T, levels,
                                     cmap=plt.cm.get_cmap('coolwarm'), alpha=0.2,
                                     vmin=bgswMin, vmax=bgswMax)

            # difference array of current frame time frame_time_num+k to
            # position time frame_time_num
            dct = frame_time_num + k - pos_time_num
            # get index of closest to 0, use this for position
            timeind = np.argmin(abs(dct))
            stopAnimation_a = False
            stopAnimation_b = False

            if spacecraft == 'AB' or spacecraft == 'A':
                [slope_a, intercept_a, tangentPoints_a, angle_a, stopAnimation_a, shockAngles_a, diffMeasure_a] = plot_ellipse(ax, dayjump, pos, timeind,
                                                      cmeind_a, k, all_apex_f_a,
                                                      all_apex_w_a, all_apex_r_a,
                                                      all_apex_lon_a, all_apex_s_a,
                                                      all_apex_flag_a, frame_time_num,
                                                      et_time_num_interp_a,
                                                      et_elon_interp_a, et_time_num_a,
                                                      startcutFit_a, endcutFit_a, tangentPoints_a,
                                                      frontData_a, ELEvoHIResults, shockAngles_a, diffMeasure_a)

            if spacecraft == 'AB' or spacecraft == 'B':
                [slope_b, intercept_b, tangentPoints_b, angle_b, stopAnimation_b, shockAngles_b, diffMeasure_b] = plot_ellipse(ax, dayjump, pos, timeind,
                                                      cmeind_b, k, all_apex_f_b,
                                                      all_apex_w_b, all_apex_r_b,
                                                      all_apex_lon_b, all_apex_s_b,
                                                      all_apex_flag_b, frame_time_num,
                                                      et_time_num_interp_b,
                                                      et_elon_interp_b, et_time_num_b,
                                                      startcutFit_b, endcutFit_b, tangentPoints_b,
                                                      frontData_b, ELEvoHIResults, shockAngles_b, diffMeasure_b)
                
#            print('SAs: ', shockAngles_a)
            
            if plotSolarWind:
                if (spacecraft == 'A' or spacecraft == 'AB') and tangentPoints_a != [] and angle_a != np.NaN:
                    #print(tangentPoints_a)
                    tpR = tangentPoints_a[-1][0]
                    tpLon = tangentPoints_a[-1][1]
                    tpLon = tpLon - np.pi - rotAngle
                    lonPoint = -(tpLon * 180 / (2 * np.pi))
                    rPoint = tpR * rLen / 2 - startBGSW

                    tpWind = []
                    rPoint = np.int(rPoint)
                    for i in np.arange(rPoint-2, rPoint+2):
                        #for j in np.arange(90, 90):
                            tpWind.append(bgswData[i, np.int(lonPoint)])
                    #print('tpWind = ', tpWind)
                    
                    tpWind_a.append([frame_time_str, np.median(tpWind)])
                    tangAngle_a.append(angle_a)
      
                if (spacecraft == 'B' or spacecraft == 'AB') and tangentPoints_b != [] and angle_a != np.NaN:
                    tpR = tangentPoints_b[-1][0]
                    tpLon = tangentPoints_b[-1][1]
                    tpLon = tpLon - np.pi - rotAngle
                    lonPoint = -(tpLon * 180 / (2 * np.pi))
                    rPoint = tpR * rLen / 2 - startBGSW

                    tpWind = []
                    rPoint = np.int(rPoint)
                    for i in np.arange(rPoint-2, rPoint+2):
                        #for j in np.arange(90, 90):
                            tpWind.append(bgswData[i, np.int(lonPoint)])
                    #print('tpWind = ', tpWind)
                    
                    tpWind_b.append([frame_time_str, np.median(tpWind)])
                    tangAngle_b.append(angle_b)
            
            if spacecraft == 'AB':
                if False:
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
                s=50, c='brown', alpha=1, lw=0, zorder=3)
            ax.scatter(
                pos.messenger[1, timeind], pos.messenger[0, timeind],
                s=25, c='dimgray', marker='s', alpha=1, lw=0, zorder=3)
            ax.scatter(
                pos.sta[1, timeind], pos.sta[0, timeind],
                s=25, c='red', alpha=1, marker='s', lw=0, zorder=3)
            ax.scatter(
                pos.stb[1, timeind], pos.stb[0, timeind],
                s=25, c='royalblue', alpha=1, marker='s', lw=0, zorder=3)
            ax.scatter(
                pos.earth[1, timeind], pos.earth[0, timeind],
                s=50, c='mediumseagreen', alpha=1, lw=0, zorder=3)
            
            ax.scatter(np.deg2rad(30), 1,
                s=50, c='black', alpha=1, lw=0, zorder=3)
            #Earth = 0.47
            plt.figtext(0.698, 0.34, 'VSC1', fontsize=labelfontsize, ha='center')
            ax.scatter(np.deg2rad(-30), 1,
                s=50, c='black', alpha=1, lw=0, zorder=3)
            plt.figtext(0.698, 0.59, 'VSC2', fontsize=labelfontsize, ha='center')
            #ax.scatter(
            #    pos.mars[1, timeind], pos.mars[0, timeind],
            #    s=50, c='orangered', alpha=1, lw=0, zorder=3)
            #ax.scatter(
            #    pos.msl[1, timeind], pos.msl[0, timeind],
            #    s=25, c='magenta', marker='s', alpha=1, lw=0, zorder=3)
            #ax.scatter(
            #    pos.maven[1, timeind], pos.maven[0, timeind],
            #    s=25, c='steelblue', marker='s', alpha=1, lw=0, zorder=3)
            #ax.scatter(
            #    pos.rosetta[1, timeind], pos.rosetta[0, timeind],
            #    s=25, c='black', marker='s', alpha=1, lw=0, zorder=3)
            #ax.scatter(
            #    pos.ulysses[1, timeind], pos.ulysses[0, timeind],
            #    s=25, c='darkolivegreen', marker='s', alpha=1, lw=0, zorder=3)

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

            #for y in range(0, np.size(active_index_vex)):
            #    # access elements in tuple that is produced by where
            #    z = active_index_vex[0][y]
            #    # fadedays is maximum difference in time, and alpha from 0 to 1
            #    fadealpha = 1 - active_icme_vex[z] / (fadedays)
            #    ax.scatter(
            #        long[z], rdist[z], s=bmean[z], c='orange',
            #        alpha=fadealpha, zorder=4)

            #for y in range(0, np.size(active_index_sta)):
            #    z = active_index_sta[0][y]
            #    # 30 days is maximum difference in time, and alpha from 0 to 1
            #    fadealpha = 1 - active_icme_sta[z] / (fadedays)
            #    ax.scatter(long[z], rdist[z], s=bmean[z], c='red',
            #               alpha=fadealpha, zorder=4)

            #for y in range(0, np.size(active_index_stb)):
            #    z = active_index_stb[0][y]
            #    # 30 days is maximum difference in time, and alpha from 0 to 1
            #    fadealpha = 1 - active_icme_stb[z] / (fadedays)
            #    ax.scatter(
            #        long[z], rdist[z], s=bmean[z], c='royalblue',
            #        alpha=fadealpha, zorder=4)

            #for y in range(0, np.size(active_index_win)):
            #    z = active_index_win[0][y]
            #    # 30 days is maximum difference in time, and alpha from 0 to 1
            #    fadealpha = 1 - active_icme_win[z] / (fadedays)
            #    ax.scatter(
            #        long[z], rdist[z], s=bmean[z], c='mediumseagreen',
            #        alpha=fadealpha, zorder=4)

            #for y in range(0, np.size(active_index_mes)):
            #    z = active_index_mes[0][y]
            #    # 30 days is maximum difference in time, and alpha from 0 to 1
            #    fadealpha = 1 - active_icme_mes[z] / (fadedays)
            #    ax.scatter(
            #        long[z], rdist[z], s=bmean[z], c='dimgrey',
            #        alpha=fadealpha, zorder=4)

            #for y in range(0, np.size(active_index_uly)):
            #    z = active_index_uly[0][y]
            #    # 30 days is maximum difference in time, and alpha from 0 to 1
            #    fadealpha = 1 - active_icme_uly[z] / (fadedays)
            #    ax.scatter(
            #        long[z], rdist[z], s=bmean[z], c='darkolivegreen',
            #        alpha=fadealpha, zorder=4)

            # ##################### legend and additional text

            #plt.suptitle('ELEvoHI ensemble simulation ')

            # Sun
            ax.scatter(0, 0, s=100, c='yellow', edgecolors='yellow')
            #plt.figtext(0.51, 0.5, 'Sun', fontsize=10, ha='center')

            # Earth
        #    plt.figtext(0.51, 0.28, 'Earth', fontsize=10, ha='center')
            plt.figtext(0.698, 0.47, 'Earth', fontsize=labelfontsize, ha='center')

            plt.figtext(0.55, 0.1, coordSysString + ' longitude',
                        fontsize=labelfontsize, ha='left')

            plt.figtext(0.1 - 0.02, 0.02, 'Mercury', color='brown',
                        ha='center', fontsize=labelfontsize*1.2)
            plt.figtext(0.26 - 0.02, 0.02, 'MESSENGER', color='dimgray',
                        ha='center', fontsize=labelfontsize*1.2)
            plt.figtext(0.41 - 0.02, 0.02, 'Venus', color='orange',
                        ha='center', fontsize=labelfontsize*1.2)
            plt.figtext(0.55 - 0.02, 0.02, 'STEREO-A', color='red',
                        ha='center', fontsize=labelfontsize*1.2)
            plt.figtext(0.72 - 0.02, 0.02, 'STEREO-B', color='royalblue',
                        ha='center', fontsize=labelfontsize*1.2)
            plt.figtext(0.85 - 0.02, 0.02, 'Earth', color='mediumseagreen',
                        ha='center', fontsize=labelfontsize*1.2)
            #plt.figtext(0.68 - 0.02, 0.02, 'Mars', color='orangered',
            #            ha='center', fontsize=labelfontsize)
            #plt.figtext(0.78 - 0.02, 0.02, 'Maven', color='steelblue',
            #            ha='center', fontsize=labelfontsize)
            #plt.figtext(0.73 - 0.02, 0.02, 'MSL', color='magenta',
            #            ha='center', fontsize=labelfontsize)
            #plt.figtext(0.84 - 0.02, 0.02, 'Rosetta', color='black',
            #            ha='center', fontsize=labelfontsize)
            #plt.figtext(0.90 - 0.02, 0.02, 'Ulysses', color='darkolivegreen',
            #            ha='center', fontsize=labelfontsize)

            
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
                fmt='%d', fontsize=labelfontsize)  # , frac = 1.05)
            ax.set_theta_zero_location('E')
            ax.set_ylim(0, 1.0) #             ax.set_ylim(0, 2.0) # 

            # plot text for date extra so it does not move
            # year
            # plt.figtext(0.52,0.85,frame_time_str[0:13]+':00:00',
            # fontsize=13, ha='center')
            plt.figtext(0.5, 0.95, frame_time_str[0:4] + frame_time_str[5:7] + frame_time_str[8:10] + '   ' +
                        frame_time_str[11:13] + ':' + frame_time_str[14:16] + ':' + frame_time_str[17:19], fontsize=labelfontsize*1.8,
                        ha='center')

            ax.grid(True, linestyle='--', linewidth=0.5)
            ax.set_rticks([0.5, 1, 1.5])#, 2])
            ax.set_rlabel_position(135)
            ax.tick_params(labelsize=labelfontsize)

            # signature
            # plt.figtext(0.95,0.01/2,r'$C. M\ddot{o}stl, T. Amerstorfer$',
            # fontsize=4, ha='center')

            if plotSolarWind:
                cax = plt.axes([0.05, 0.2, 0.02, 0.6])
                cbar = plt.colorbar(cf, cax=cax, ticks=np.arange(200, 800, 50))
                ticklabs = cbar.ax.get_yticklabels()
                cbar.ax.set_yticklabels(ticklabs, fontsize=labelfontsize)
                cbar.set_label('Solar wind speed [km/s]', fontsize=labelfontsize)
                
            # ##################### save frame
            plt.savefig(
                current_event_dir + '/frames/elevohi_' + framestr + '.png',
                dpi=300)
            # clears plot window
            plt.clf()
            
            if spacecraft == 'A':
                if stopAnimation_a:
                    os.remove(current_event_dir + '/frames/elevohi_' + framestr + '.png')
                    break
            if spacecraft == 'B':
                if stopAnimation_b:
                    os.remove(current_event_dir + '/frames/elevohi_' + framestr + '.png')
                    break
            if spacecraft == 'AB':
                if stopAnimation_a and stopAnimation_b:
                    os.remove(current_event_dir + '/frames/elevohi_' + framestr + '.png')
                    break

                    
#            asdf
        # ########### end of loop
        
#        print('ShockAngles: ', np.size(shockAngles_a))
#        print('min: ', np.min(shockAngles_a))
#        print('max: ', np.max(shockAngles_a))
#        print('median: ', np.median(shockAngles_a))
#        print('mean: ', np.mean(shockAngles_a))

        
    
        if outPath is None:
            outpath = current_event_dir

        if not os.path.isdir(outpath):
            os.mkdir(outpath)

        outpath = outpath + '/'
        
        if spacecraft == 'AB':
            diffMeasure_a = np.array(diffMeasure_a)
            diffMeasure_b = np.array(diffMeasure_b)
            pickle.dump((tpWind_a, tpWind_b, et_time_num_interp_a, et_time_num_interp_b, tangAngle_a, 
                         tangAngle_b, tangentPoints_a, tangentPoints_b, shockAngles_a, shockAngles_b,
                         diffMeasure_a, diffMeasure_b),
                open(current_event_dir + "tpWind_AB.p", "wb"))
        if spacecraft == 'A':
            diffMeasure_a = np.array(diffMeasure_a)
            pickle.dump((tpWind_a, et_time_num_interp_a, tangAngle_a, tangentPoints_a, shockAngles_a, diffMeasure_a),
                open(current_event_dir + "tpWind_A.p", "wb"))
        if spacecraft == 'B':
            diffMeasure_b = np.array(diffMeasure_b)
            pickle.dump((tpWind_b, et_time_num_interp_b, tangAngle_a, tangentPoints_b, shockAngles_b, diffMeasure_b),
                open(current_event_dir + "tpWind_B.p", "wb"))


        # convert to jpg
        # os.system('ffmpeg -i "'+current_event_dir+'frames/elevo_%04d.png" '
        # +current_event_dir+'frames/elevo_%04d.jpg -y -loglevel quiet')
        # make mp4

        os.system(ffmpegPath + 'ffmpeg -r 20 -i "' + current_event_dir +
                  'frames/elevohi_%04d.png" -c:v libx264 -vf "fps=25,format=yuv420p" ' + outpath +
                  current_event + '_' + spacecraft +
                  '_ensemble_movie.mp4 -y -loglevel quiet')
        plt.close('all')
        
#        plot_BGSW_tangent(current_event_dir)

        print('Made movie.')
        print('End ElEvoHI animation program.')
        print('The run took: ', datetime.datetime.now() - startTime)




# ########################## MIT license




if __name__ == '__main__':

    eventslist = ['20100203']
    
    model = 'HUX'
    main(eventslist, spaceCraft='AB', scriptPath='HI_animate/events/DefFront/'+model+'/',
         catPath='HI_animate/cats/', readData=0, plotBGSW=True, bgswModel=model)

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
