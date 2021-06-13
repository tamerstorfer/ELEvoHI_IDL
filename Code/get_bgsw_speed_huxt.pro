; Name:		get_bgsw_speed_huxt
;
; Purpose: 	Returns the speed of the ambient solar wind with respect to the radial distance and the given time
;
; Calling sequence:	sw_speed = get_bgsw_speed_huxt(bgswdata=bgswdata, tdrag[j], lonNew[i], rbefore[i])
;
; Parameters (input):
;			bgswdata: ambient solar wind data for a full Carrington rotation from HUXt
;			timeStep: time for which the speed is needed
;			lon: longitude for which the speed is needed
;           radius: radius for which the speed is needed
;
; Parameters (output):
;			speed: ambient solar wind speed for the given time, radial distance and longitude
;
; History:    2021/03: created (Juergen Hinterreiter)
;
; Authors:    Tanja Amerstorfer & Christian Moestl & Juergen Hinterreiter
;             Space Research Institute, Austrian Academy of Sciences
;			  Graz, Austria
; -
function get_bgsw_speed_huxt, bgswdata=bgswdata, timeStep, lon, radius

    windStartTime = bgswData.tinit
    windTimes = bgswData.time + windStartTime
    rad = round(radius)
    bgswLon = bgswData.lon
    
    lon = -1*lon*!dtor
    rad = rad-bgswData.r[0]
    
    if lon lt 0 then bgswLon = bgswLon - 2*!pi-bgswData.dlon
    

    timeInd = where(abs(windTimes-timeStep) eq min(abs(windTimes-timeStep)))
    lonInd = where(abs(bgswLon-lon) eq min(abs(bgswLon-lon)))

    return, bgswData.v_grid[lonInd, rad, timeInd]
end