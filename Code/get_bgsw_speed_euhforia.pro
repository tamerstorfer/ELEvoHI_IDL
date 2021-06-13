; Name:		get_bgsw_speed_euhforia
;
; Purpose: 	Returns the speed and density of the ambient solar wind of a specific location in the heliosphere for the given time
;
; Calling sequence:	sw_speed = get_bgsw_speed_euhforia(bgswdata=bgswdata, tdrag[j], lonNew[i], rbefore[i])
;
; Parameters (input):
;			bgswdata: ambient solar wind data for a full Carrington rotation from EUHFORIA
;			timeStep: time for which the speed is needed
;			lon: longitude for which the speed is needed
;           radius: radius for which the speed is needed
;
; Parameters (output):
;			[sw_speed, sw_density]: ambient solar wind speed and density for the given time, radial distance and longitude
;
; History:    2021/05: created (Juergen Hinterreiter)
;
; Authors:    Tanja Amerstorfer & Christian Moestl & Juergen Hinterreiter
;             Space Research Institute, Austrian Academy of Sciences
;			  Graz, Austria
; -
function get_bgsw_speed_euhforia, bgswdata=bgswdata, timeStep, lon, radius

    tDiff = (timeStep - bgswData.tinit)/(60.*60.*24.)
    rotSun = 27.27 ; days
    
    lonSize = (size(bgswData.varr))[1]
    lonNew = -1*lon
    
    lonEarth = bgswData.lon + lonSize/rotSun*tDiff
    
    indLon = where(abs(lonEarth-lonNew) eq min(abs(lonEarth-lonNew)))
    indR = where(abs(bgswData.r-radius) eq min(abs(bgswData.r-radius)))
    
    return, [bgswData.varr[indLon, indR], bgswData.narr[indLon, indR]*1.67d-24*1d15]
end