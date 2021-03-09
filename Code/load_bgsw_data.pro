;+
;
; Name:       laod_bgsw_daat
;
; Purpose:    loads the ambient solar wind data from the WSA-HUX model for the given event
;
; Calling sequence: load_bgsw_data, bgsw_file, bgswData=bgswData, bgswTime=bgswTime
;
; Parameters (input):
;			  bgsw_file: Path to the ambient solar wind file
;			  bgswData: ambient solar wind speed for a full Carrington rotation
;			  bgswTime: Time for which the ambient solar wind was crated
;
;
; History:    2021/01: created (Juergen Hinterreiter)
;
; Authors:    Tanja Amerstorfer & Christian Moestl & Juergen Hinterreiter
;             Space Research Institute, Austrian Academy of Sciences
;			  Graz, Austria
; -

pro load_bgsw_data, file, bgswData=bgswData, bgswTime=bgswTime
    
    eventDate = strmid(file, strpos(file, 'vmap')-11, 8)
    array = ''
    line = ''
    i = 0

    OPENR, lun, file, /GET_LUN
    READF, lun, line 
    lineSplit = STRSPLIT(line, /EXTRACT)
    FREE_LUN, lun


    data = fltarr(n_elements(lineSplit), file_lines(file))
    OPENR, lun, file, /GET_LUN
    WHILE NOT EOF(lun) DO BEGIN
        READF, lun, line 
        lineSplit = STRSPLIT(line, /EXTRACT)
        data[*,i] = lineSplit
        i = i+1
    ENDWHILE
  ; Close the file and free the file unit
    FREE_LUN, lun

    ;n_elements(data[*,0]) (longitude)
    ;n_elements(data[0,*]) (Radius)

    strPos = strpos(file, '/vmap.txt')
    timeFile = strmid(file, 0, strPos-10) + 'dataset-overview.txt'
    timeFile = file_search(timeFile)

    yr = strmid(eventDate, 0, 4)
    mon = strmid(eventDate, 4, 2)
    dy = strmid(eventDate, 6, 2)
    tinitNum = anytim(yr+'-'+mon+'-'+dy)


    readcol, timefile, label, bgswDate, bgswTime, bgswCR, format='A,A,A,A', skipline=2

    bgswDateTime = bgswDate + ' ' + bgswTime
    timeDiff = abs(anytim(bgswDateTime) - tinitNum)
    idTimeDiff = where(timeDiff eq min(timeDiff))

    bgswTime = bgswDateTime[idTimeDiff]
    bgswData = data

end