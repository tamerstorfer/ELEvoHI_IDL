; tinit: time, when speed is needed
; lon: longitude with respect to Earth (0° is Earth direction, negative values are to the West, positive values to the East)
; rad: distance to the sun [r_sun]


; Name:		get_bgsw_speed
;
; Purpose: 	Returns the speed of the ambient solar wind with respect to the radial distance and the given time
;
; Calling sequence:	sw_speed = get_bgsw_speed(bgswdata=bgswdata, bgswTimeNum, tinitnum, lonNew[i], rbefore[i])
;
; Parameters (input):
;			bgswdata: ambient solar wind data for a full Carrington rotation
;			bgswTime: time for which the ambient solar wind data was created (needed for correct rotation)
;			tinit: time at which the speed is needed
;			lonInput: longitude with respect to Earth (0° is Earth direction, negative values are to the West, positive values to the East)
;			radInput: distance to the sun [r_sun]
;			plotWind: keyword to plot the ambient solar wind
;
; Parameters (output):
;			posSpeed: ambient solar wind speed for the given time, radial distance and longitude
;
; History:    2021/03: created (Juergen Hinterreiter)
;
; Authors:    Tanja Amerstorfer & Christian Moestl & Juergen Hinterreiter
;             Space Research Institute, Austrian Academy of Sciences
;			  Graz, Austria
; -


function get_bgsw_speed, bgswdata=bgswData, bgswTime, tinit, lonInput, radInput, plotWind=plotWind


    data = bgswData
    lon = lonInput/2.

    modelStartDist = 5. ; in R_sun
    earthPos = 30. ; in pixels
    rotSun = 27.27

    rad = radInput-modelStartDist

    daysLater = rotSun/4
    daysLater = 1

    if keyword_set(earthLon) then earthPos = earthLon

    tinitNum = tinit

    timeDiff = tinitNum - bgswTime
    tDiff = timeDiff/(60.*60.*24.)

    dSize = size(data)


    shiftVal = dSize[1]/rotSun*tDiff + lon; only shift because of time difference to model output of HUX
    ; add also shift because of longitude
;    shiftVal = shiftVal + lon

    posSpeed = data[earthPos + shiftval, rad]

    return, posSpeed
end
