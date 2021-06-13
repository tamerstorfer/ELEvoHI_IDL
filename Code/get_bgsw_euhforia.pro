; Name:		get_bgsw_euhforia
;
; Purpose: 	Returns a range around the median speed of the ambient solar wind from a given region
;
; Calling sequence: winds = get_bgsw_euhforia(bgswData, time[cut], scut, r_apex_sun[ecut], phi, lambda, sc)
;
; Parameters (input):
;			bgswData: ambient solar wind data
;			tinit: time for the SC positions
;           startcut: startcut of the DBM-fit
;           endcut: endcut of the DBM-fit
;           maxPhi: propagation direction
;           halfwidth: halfwidth of the CME
;           sc: spacecraft ('A' or 'B')
;
; Parameters (output):
;			solarWind: returns 9 different ambient solar wind speed for the given area (dependent on lambda, phi, startcut, endcut)
;
; History:    2021/05: created (Juergen Hinterreiter)
;
; Authors:    Tanja Amerstorfer & Christian Moestl & Juergen Hinterreiter
;             Space Research Institute, Austrian Academy of Sciences
;			  Graz, Austria
; -

function get_bgsw_euhforia, bgswData, tinit, startcut, endcut, phi, halfWidth, sc

    tinitnum = anytim(tinit)
    
    lonDeg = bgswData.lon
    lonRes = lonDeg[1]-lonDeg[0] ;resolution of the longitude in °
    r = bgswData.r
    timeDiff = tinitnum - bgswData.tinit
    varr = bgswData.varr
    tDiff = timeDiff/(60.*60.*24.)
    rotSun = 27.27 ; days
    
    dSize = size(bgswData.varr)
    lonSize = dsize[1]
    
    shiftVal = lonSize/rotSun*tDiff ;only shift because of time difference to model output of HUX

    print, '!!!! Background Solar Wind from EUHFORIA!!!!'
    MAError = 100 ;/km/s
    stSize = 25 ;km/s
    minimumSW = 233 ;km/s

    pos_E=get_stereo_lonlat(tinit, 'Earth', system='HEE')
    pos_A=get_stereo_lonlat(tinit, 'Ahead', system='HEE')
    pos_B=get_stereo_lonlat(tinit, 'Behind', system='HEE')
    
    if sc eq 'A' then begin
        sep=abs(pos_E[1]-pos_A[1])/!dtor  
        dir_E=(sep-phi)
    endif
    if sc eq 'B' then begin
        sep=abs(pos_E[1]-pos_B[1])/!dtor  
        dir_E=-(sep-phi)
    endif 
    
    lonApex = lonDeg + shiftval - dir_e
    indApex = where(abs(lonApex) eq min(abs(lonApex)))
    rMin = where(abs(r-startcut) eq min(abs(r-startcut)))
    rMax = where(abs(r-endcut) eq min(abs(r-endcut)))
    print, indApex
    
    lonMin = indApex-halfwidth*lonRes
    lonMax = indApex+halfwidth*lonRes
  
    bgswRange = varr[lonMin:lonMax, rMin:rMax]
    windSpeed = median(bgswRange)
    
    print, "Median background solar wind: ", windSpeed
    nrElements = (2*MAError/stSize) + 1
    solarWind = findgen(nrElements)*stSize+windSpeed-MAError
    ; remove all the SW values that are smaller than minimumSW
    ; and set the first value to minimumSW
    idMinSW = where(solarWind lt minimumSW)
    if idMinSW[0] ne -1 then begin
        remove, idMinSW, solarWind
        if min(solarWind)-stSize/2 gt minimumSW then solarWind = [minimumSW, solarWind]
    endif
    
    return, solarWind
end
