
FUNCTION CIRCLE, xcenter, ycenter, radius
   points = (2 * !PI / 99.0) * FINDGEN(100)
   x = xcenter + radius * COS(points )
   y = ycenter + radius * SIN(points )
   RETURN, TRANSPOSE([[x],[y]])
   END


; Name:		get_bgsw
;
; Purpose: 	Returns a range around the median speed of the ambient solar wind from a given region
;
; Calling sequence: winds = get_bgsw_hux(bgsw_file, time[cut], scut, r_apex_sun[ecut], phi, phi, lambda, sc)
;
; Parameters (input):
;			file: ambient solar wind file
;			tinit: time for the SC positions
;           startcut: startcut of the DBM-fit
;           endcut: endcut of the DBM-fit
;           minPhi: minimum propagation direction
;           maxPhi: maximum propagation direction
;           halfwidth: halfwidth of the CME
;           sc: spacecraft ('A' or 'B')
;			lonInput: longitude with respect to Earth (0° is Earth direction, negative values are to the West, positive values to the East)
;			radInput: distance to the sun [r_sun]
;			plotWind: keyword to plot the ambient solar wind
;
; Parameters (output):
;			solarWind: returns 9 different ambient solar wind speed for the given area (dependent on lambda, phi, startcut, endcut)
;
; Keywords:
;		    stepSize: step size of the ambient solar wind output (default: 25 km/s)
;           MAE: mean absolute error (default 100 km/s)
;           savePlot: set to create plot
;           earthLon: Longitude of Earth in the Model (default 30 [corresponds to a longitude of 60°])
;           minSW: minimum solar wind speed which is the lower limit for the output (default 233 km/s)
;           plotPath: path where to plot when savePlot is set
;           saveData: save the data also to plotPath
;
; History:    2021/03: created (Juergen Hinterreiter)
;
; Authors:    Tanja Amerstorfer & Christian Moestl & Juergen Hinterreiter
;             Space Research Institute, Austrian Academy of Sciences
;			  Graz, Austria
; -

function get_bgsw_hux, file, tinit, startcut, endcut, minphi, maxphi, halfWidth, sc, stepSize=stepSize, MAE = MAE, savePlot=savePlot, earthLon = earthLon, minSW = minSW, plotPath = plotPath, saveData = saveData

  eventDate = strmid(file, strpos(file, 'vmap')-11, 8)

  print, '!!!! Background Solar Wind from WSA model !!!!'
  MAError = 100 ;/km/s
  modelStartDist = 5
  startcut = startcut - modelStartDist
  endcut = endcut - modelStartDist

  if startcut lt 0 then startcut = 0
  if endcut lt 0 then endcut = 0
  if keyword_set(MAE) then MAError = MAE

  stSize = 25 ;km/s
  if keyword_set(stepSize) then stSize = stepSize

  earthPos = 30 ; in pixels
  if keyword_set(earthLon) then earthPos = earthLon

  minimumSW = 233 ;km/s
  if keyword_set(minSW) then minimumSW = minSW



  pos_E=get_stereo_lonlat(tinit, 'Earth', system='HEE')
  pos_A=get_stereo_lonlat(tinit, 'Ahead', system='HEE')
  pos_B=get_stereo_lonlat(tinit, 'Behind', system='HEE')
    

  phi_mean = mean([minphi, maxphi])
  ;calculate direction from Earth
  ; divide by 2 because it is defined in degree
  if sc eq 'A' then begin
    sep=abs(pos_E[1]-pos_A[1])/!dtor  
    dir_E=fix((sep-phi_mean)/2)
  endif
  if sc eq 'B' then begin
    sep=abs(pos_E[1]-pos_B[1])/!dtor  
    dir_E=fix(-(sep-phi_mean)/2)
  endif  

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
    ;print, i
    data[*,i] = lineSplit
    i = i+1
  ENDWHILE
  ; Close the file and free the file unit
  FREE_LUN, lun

  dSize = size(data)

  ; shift matrix so 'Earth' is in center
  phiDiff = ((maxphi-minphi)/2)/2
  shiftVal = dSize[1]/2-earthPos+dir_E

  newdata = data
  newData = shift(newData, [shiftVal, 0])

  ; divide by 2 because it is defined in degree
  hw = halfWidth/2

  newbgsw = newdata[earthPos-hw+shiftVal-dir_E-phiDiff:earthPos+hw+shiftVal-dir_E+phiDiff, startcut:endcut]

  newswarr = reform(newbgsw, (size(newbgsw))[4])

  solarWind = median(newbgsw)

  maxVal = max(data)
  dataPolarPlot = newData
  if 1 eq 1 then begin
    newdata[earthPos-hw+shiftVal-dir_E-phiDiff,startcut:endcut] = maxVal
    newdata[earthPos-hw+shiftVal-dir_E-phiDiff:earthPos+hw+shiftVal-dir_E+phiDiff, startcut] = maxVal
    newdata[earthPos+hw+shiftVal-dir_E+phiDiff, startcut:endcut] = maxVal
    newdata[earthPos-hw+shiftVal-dir_E-phiDiff:earthPos+hw+shiftVal-dir_E+phiDiff, endcut] = maxVal
    ;newdata[earthPos+shiftval, 0:400] = maxVal
;    newData[earthPos+shiftval-dir_E, startcut:endcut] = min(data)
  endif

  ; shift data back
  newdata = shift(newdata, [-dir_E,0])
  
  if keyword_set(saveData) and keyword_set(plotPath) then begin
      strPos = strpos(file, '/vmap.txt')
      timeFile = strmid(file, 0, strPos-10) + 'dataset-overview.txt'
      timeFile = file_search(timeFile)

      yr = strmid(eventDate, 0, 4)
      mon = strmid(eventDate, 4, 2)
      dy = strmid(eventDate, 6, 2)
      evDateNum = anytim(yr + '-' + mon + '-' + dy)


      readcol, timefile, label, bgswDate, bgswTime, bgswCR, format='A,A,A,A', skipline=2

      bgswDateTime = bgswDate + ' ' + bgswTime
      timeDiff = abs(anytim(bgswDateTime) - evDateNum)
      idTimeDiff = where(timeDiff eq min(timeDiff))

      bgswStartTime = bgswDateTime[idTimeDiff]
      

      bgsw_data = newData
      save, bgsw_data, bgswStartTime, filename = plotPath + '/bgsw_data.sav'
      bgsw_data = dataPolarPlot
      save, bgsw_data, bgswStartTime, filename = plotPath + '/bgsw_data_noCutout.sav'
      bgsw_cutout = newbgsw
      save, bgsw_cutout, bgswStartTime, filename = plotPath + '/bgsw_cutout_only.sav'
  endif

  


  if keyword_set(savePlot) then begin
    plotP = ''
    if keyword_set(plotPath) then plotP = plotPath + '/'

    !p.multi = 0
    set_plot, 'ps'
    DEVICE, SET_FONT='Helvetica', /TT_FONT
    device, filename = plotP + 'bgsw_'+eventDate+'.ps', xsize = 15, ysize = 25, xoffset = 2, yoffset = 2, encaps = 0


    !p.multi = 0
    page_height = 27.94
    page_width = 21.59

    ; define position of plot
    ; bottom left corner (in cm)
    plot_left = 2
    
    ; define plot size (in cm)
    xsize = 18
    ysize = 10

    plot_bottom = page_height - ysize

    xticks = indgen(19)*20-180

    ; Plot BGSW 
    Loadct, 25
    cgimage, newdata, $
    position = [plot_left / page_width, plot_bottom / page_height, $
      (plot_left + xsize) / page_width, (plot_bottom + ysize) / page_height] 

    ; draw axes
    cgplot, [0], [0], $
    ;xcharsize = 1.0, $
    ;ycharsize = 1.0, $
    charsize = 1.0, $
    thick = 5, $
    ;xrange = [0, dSize[1]], $
    xrange = [dSize[1], -dSize[1]], $
    yrange = [modelStartDist, dSize[2]+modelStartDist], $
    xtitle = 'Longitude', $
    ;ytitle = 'Distance [R!DSun!N]', $
    ytitle = 'Distance [R'+sunsymbol()+']', $
    title = eventDate, $
    xstyle = 1, ystyle = 1, /nodata, /noerase, $
    position = [plot_left / page_width, plot_bottom / page_height, $
      (plot_left + xsize) / page_width, (plot_bottom + ysize) / page_height]

    cgColorbar, Range=[Min(data), Max(data)], /Vertical, /fit, /right, title = 'Speed [km s!E-1!N]', charsize = 1.0

    wspd = dataPolarPlot

    for i = 0, 179 do begin
      ;wspd[i,*] = dataPolarPlot[179-i, *]
      wspd[i,*] = newData[179-i, *]
    endfor
    
    lon = findgen(dSize[1])*2
    lat = findgen(dSize[2])

    lrad = lon/!radeg
    r = lat

    ; shift to look the same like Martins plots
    wspd = shift(wspd, [dSize[1]/2+1, 0])

    ysize = xsize*page_width/page_height
    plot_bottom = page_height - ysize
    

    ; Plot BGSW as polar plot
    diff = fix(max(data) - min(data))
    levs = indgen(diff)+min(data)
    loadct, 0
    POLAR_CONTOUR, wspd, lrad, r, /nodata, background = 255, color = cgcolor('black'), $
      position = [plot_left / page_width, plot_bottom / page_height, $
      (plot_left + xsize) / page_width, (plot_bottom + ysize) / page_height], xtitle = 'Distance [R'+sunsymbol()+']', ytitle = 'Distance [R'+sunsymbol()+']', title = eventDate
    loadct, 25
    POLAR_CONTOUR, wspd, lrad, r, /fill, levels = levs, /over, $
      position = [plot_left / page_width, plot_bottom / page_height, $
      (plot_left + xsize) / page_width, (plot_bottom + ysize) / page_height];, /isotropic, xgridstyle = 1, ygridstyle = 1
    ;cgColorbar, Range=[Min(data), Max(data)], /Vertical, position = [plot_left / page_width, plot_bottom / page_height, $
    ;  (plot_left + xsize) / page_width, (plot_bottom + ysize) / page_height], /normal, ncolors = 100
    cgColorbar, Range=[Min(data), Max(data)], /Vertical, position = [0.575, 0.815, 0.927, 0.835], /right, charsize = 1, title = 'Speed [km s!E-1!N]';, ncolors = 100
    ;cgColorbar, Range=[Min(data), Max(data)], /Vertical, position = [0.186, 0.9, 0.839, 0.95], /data;, ncolors = 100
    ;colorbar, range = [Min(data), Max(data)], /vertical, position = [-400, 400, 400, 500], /data

    diagVal = max(r)/sqrt(2)

    loadct, 0
    plots, [0, 0], [-max(r), max(r)], color = cgcolor('black'), /data, linestyle = 2
    plots, [-max(r), max(r)], [0, 0], color = cgcolor('black'), /data, linestyle = 2
    plots, [-diagVal, diagVal], [-diagVal, diagVal], color = cgcolor('black'), /data, linestyle = 2
    plots, [-diagVal, diagVal], [diagVal, -diagVal], color = cgcolor('black'), /data, linestyle = 2


    plots, circle(0, 0, max(r)), /data, color = cgcolor('black'), linestyle = 2
    plots, circle(0, 0, 215), /data, color = cgcolor('black'), linestyle = 2
    cgColorFill, CIRCLE(215, 0, 7), /data, color = cgcolor('black')
    xyouts, 215, 10, 'Earth', color = cgcolor('black'), /data, charsize = 1


    ; Plot histogram
    !p.multi = [0,1,3]
    cghistoplot, newswarr, title = 'Histogram of selected range of the BGSW', xtitle = 'Speed [km s!E-1!N]'
    !p.multi = 0

    device, /close
    set_plot, 'x'
  endif

  print, "Median background solar wind: ", solarWind
  nrElements = (2*MAError/stSize) + 1
  solarWind = findgen(nrElements)*stSize+solarWind-MAError



  
  ; remove all the SW values that are smaller than minimumSW
  ; and set the first value to minimumSW
  idMinSW = where(solarWind lt minimumSW)
  if idMinSW[0] ne -1 then begin
    remove, idMinSW, solarWind
    if min(solarWind)-stSize/2 gt minimumSW then solarWind = [minimumSW, solarWind]
  endif
  
  print, solarWind
  return, solarWind

end