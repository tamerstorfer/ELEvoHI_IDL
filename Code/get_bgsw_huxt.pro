; Name:		get_bgsw_huxt
;
; Purpose: 	Returns a range around the median speed of the ambient solar wind from a given region
;
; Calling sequence: winds = get_bgsw_huxt(bgswData, time[cut], time[ecut], scut, r_apex_sun[ecut], phi, lambda, sc)
;
; Parameters (input):
;			bgswData: ambient solar wind data
;			tstart: starttime of the DBM-fit
;           tend: endtime of the DBM-fit
;           startcut: startcut of the DBM-fit
;           endcut: endcut of the DBM-fit
;           phi: propagation direction
;           lambda: halfwidth of the CME
;           sc: spacecraft ('A' or 'B')
;
;           Pre-defined (as in 'get_bgsw'):
;               minimumSW = 233; km/s
;               MAError = 100; km/s
;               stSize = 25 ;km/s
;
; Parameters (output):
;			solarWind: returns 9 different ambient solar wind speed for the given area (dependent on lambda, phi, startcut, endcut)
;
;
; History:    2021/03: created (Juergen Hinterreiter)
;
; Authors:    Tanja Amerstorfer & Christian Moestl & Juergen Hinterreiter
;             Space Research Institute, Austrian Academy of Sciences
;			  Graz, Austria
; -

function get_bgsw_huxt, bgswData, tstart, tend, startcut, endcut, phi, lambda, sc

    tinit = tstart
    scut = round(startcut) - bgswData.r[0]
    ecut = round(endcut)- bgswData.r[0]

    pos_E=get_stereo_lonlat(tinit, 'Earth', system='HEE')
    pos_A=get_stereo_lonlat(tinit, 'Ahead', system='HEE')
    pos_B=get_stereo_lonlat(tinit, 'Behind', system='HEE')

    print, phi
    ;calculate direction from Earth
    if sc eq 'A' then begin
        sep=abs(pos_E[1]-pos_A[1])/!dtor  
        dir_E=(sep-phi)
    endif
    if sc eq 'B' then begin
        sep=abs(pos_E[1]-pos_B[1])/!dtor  
        dir_E=-(sep-phi)
    endif
  
    minLimit = (dir_E - lambda)*!dtor
    maxLimit = (dir_E + lambda)*!dtor
    
    lon = bgswData.lon
    lonSize = n_elements(lon)-1
    cr_lon_init = bgswData.cr_lon_init
    lons = lon;-cr_lon_init
    windStartTime = bgswData.tinit
    windStartTime = anytim(windStartTime)
    windTimes = bgswData.time + windStartTime
    EarthDir = lons[0]
    
    speeds = bgswData.v_grid
  
    timeMin = where(abs(windTimes-anytim(tstart)) eq min(abs(windTimes-anytim(tstart))))
    timeMax = where(abs(windTimes-anytim(tend)) eq min(abs(windTimes-anytim(tend))))
  
    lonMax = where(abs(lons-maxLimit) eq min(abs(lons-maxlimit)))
    lonsMinus = lons-lons[lonSize]
    lonMin = where(abs(lonsMinus-minLimit) eq min(abs(lonsMinus-minlimit)))
    
    cutoutSpeeds = speeds[0:lonMax, scut:ecut, timeMin:timeMax]
    cutoutSpeeds = [cutoutSpeeds, speeds[lonMin:lonSize, scut:ecut, timeMin:timeMax]]
    
    
    solarWind = median(cutoutspeeds)
    print, "Median background solar wind: ", solarWind
    
    minimumSW = 233; km/s
    MAError = 100; km/s
    stSize = 25 ;km/s
    nrElements = (2*MAError/stSize) + 1
    solarWind = findgen(nrElements)*stSize+solarWind-MAError

  
    ; remove all the SW values that are smaller than minimumSW
    ; and set the first value to minimumSW
    idMinSW = where(solarWind lt minimumSW)
    if idMinSW[0] ne -1 then begin
        remove, idMinSW, solarWind
        if min(solarWind)-stSize/2 gt minimumSW then solarWind = [minimumSW, solarWind]
    endif

    print, 'mean, stddev, min, max'
    print, [median(solarWind), stddev(solarWind), min(solarWind), max(solarWind)]
    
    return, solarwind
  
end