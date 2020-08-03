;+
; Project     : STEREO - SSC
;
; Name        : SSC_PLOT_WHERE
;
; Purpose     : Plot the positions of the two STEREO spacecraft.
;
; Category    : STEREO, Orbit, Graphics
;
; Explanation : Plots the ecliptic positions of the two STEREO spacecraft on
;               top of an planetary orbital chart.  Early in the mission, it
;               will also zoom in on the Earth-Moon system.
;
; Syntax      : SSC_PLOT_WHERE, DATE
;
; Examples    : GET_UTC, UTC  &  SSC_PLOT_WHERE, UTC
;
; Inputs      : DATE    = The date/time to be plotted.  Can be in any format
;                         supported by ANYTIM2UTC.
;
; Opt. Inputs : None.
;
; Outputs     : A plot is created showing the spacecraft positions.
;
; Opt. Outputs: None.
;
; Keywords    : INNER     = If set, then expand out to the orbit of Mars.
;               OUTER     = If set, then expand out to the orbit of Saturn.
;
;               PARKER    = If set, then overplot the Parker spiral.
;
;               ULYSSES   = If set, then overplot orbit of Ulysses.
;               MESSENGER = If set, then overplot orbit of Messenger.  The
;                           orbit drawn is Keplerian, and does not include any
;                           deviations due to maneuvers or planet fly-bys.
;               ROSETTA   = If set, then overplot orbit of Rosetta, or Comet
;                           67P Churyumov-Gerasimenko, depending on date.
;               VESTA     = If set, then overplot orbit of Vesta (for DAWN)
;               CERES     = If set, then overplot orbit of Ceres (for DAWN)
;               MAVEN     = If set, then overplot orbit of MAVEN (projected).
;                           Only applicable during cruise phase.  Selects
;                           /INNER automatically.
;               BEPICOLOMBO = If set, then overplot orbit of BEPICOLOMBO
;                           (projected).  Only applicable during cruise phase.
;                           Selects /INNER automatically.
;               PSP       = If set, then overplot orbit of Parker Solar Probe.
;                           Selects /INNER automatically.
;               SOLO      = If set, then overplot orbit of Solar Orbiter.
;                           Selects /INNER automatically.
;
;               WHITE_BG  = If set, then plot against a white background.
;
;               XRANGE    = Explicit range along X (vertical) axis.  Ignored if
;                           /INNER, /OUTER, or certain missions selected.
;               YRANGE    = Explicit range along Y (horizontal) axis.  Ignored
;                           if /INNER, /OUTER, or certain missions selected.
;
;               Will also accept any keywords for LOAD_STEREO_SPICE or
;               SSC_OPLOT_PARKER.
;
; Calls       : LOAD_STEREO_SPICE, ANYTIM2UTC, GET_STEREO_COORD,
;               SSC_FORM_ORBIT, LINECOLOR, CIRCLE_SYM, SETVIEW, SETSCALE,
;               CONCAT_DIR
;
; Common      : None.
;
; Restrictions: None.
;
; Side effects: Loads the STEREO ephemerides if not already loaded.
;
; Prev. Hist. : None.
;
; History     : Version 1, 23-Jun-2006, William Thompson, GSFC
;               Version 2, 27-Jun-2006, William Thompson, GSFC
;                       Added keywords ULYSSES, MESSENGER, PARKER
;               Version 3, 05-Jul-2006, William Thompson, GSFC
;                       Plot positive values along Y axis (X_HEE)
;               Version 4, 19-Jul-2006, William Thompson, GSFC
;                       Call SSC_OPLOT_PARKER with /REVERSE
;                       Lightened color of Parker overplot
;                       Rearranged colors to be consistent with beacon plots
;               Version 5, 10-Aug-2006, William Thompson, GSFC
;                       Added keyword WHITE_BG.  Add axis labels.
;               Version 6, 03-Nov-2006, William Thompson, GSFC
;                       Add orbit trace on zoomed-in plot, sun vector.
;               Version 7, 09-Nov-2006, William Thompson, GSFC
;                       Manipulate colors based on background color.
;               Version 8, 17-Mar-2009, WTT, change logic for twin plots
;               Version 9, 24-Sep-2010, WTT, changed name of Ulysses file
;               Version 10, 4-Aug-2011, WTT, added VESTA and CERES keywords
;               Version 11, 8-Aug-2011, WTT, automatic scaling for VESTA/CERES
;               Version 12, 13-Sep-2011, WTT, better MU value
;               Version 13, 04-Feb-2013, WTT, added keyword ROSETTA
;               Version 14, 22-Aug-2018, WTT, added keyword PSP
;               Version 15, 28-Aug-2018, WTT, added keyword SOLO
;               Version 16, 05-Nov-2018, WTT, PSP, SO error checking
;               Version 17, 13-Nov-2019, WTT, added keyword BEPICOLOMBO
;               Version 18, 26-Feb-2019, WTT, change filename for BEPICOLOMBO
;
; Contact     : WTHOMPSON
;-
;
pro ssc_plot_where_elevo, date, inner=inner, outer=outer, ulysses=kulysses, $
                    messenger=kmessenger, vesta=kvesta, ceres=kceres, $
                    rosetta=krosetta, maven=kmaven, psp=kpsp, parker=parker, $
                    white_bg=white_bg, xrange=xrange, yrange=yrange, $
                    solo=ksolo, bepicolombo=kbepicolombo, _extra=_extra
on_error, 2
;
;  Make sure the ephemeris files are loaded.
;
load_stereo_spice, _extra=_extra
mu = 1.32712440018D11   ;G*Msun km^3/s^2
;
;  Get the positions of the relevant solar system bodies.
;
utc = anytim2utc(date, /ccsds)
mercury = get_stereo_coord(utc, /au, system='HEE', 'Mercury')
venus   = get_stereo_coord(utc, /au, system='HEE', 'Venus')
earth   = get_stereo_coord(utc, /au, system='HEE', 'Earth')
moon    = get_stereo_coord(utc, /au, system='HEE', 'Moon')
mars    = get_stereo_coord(utc, /au, system='HEE', 'Mars')
jupiter = get_stereo_coord(utc, /au, system='HEE', 'Jupiter Barycenter')
saturn  = get_stereo_coord(utc, /au, system='HEE', 'Saturn Barycenter')
uranus  = get_stereo_coord(utc, /au, system='HEE', 'Uranus Barycenter')
;
;  Get the positions of STEREO A & B
;
sta = get_stereo_coord(utc, /au, system='HEE', 'STEREO Ahead')
stb = get_stereo_coord(utc, /au, system='HEE', 'STEREO Behind')
;
;  Calculate the plot range.
;
if keyword_set(outer) then begin
    x = [-11,11]
    y = [-11,11]
end else if keyword_set(inner) then begin
    x = [-1.7,1.7]
    y = [-1.7,1.7]
end else if keyword_set(kceres) then begin
    x = [-3.1,3.1]
    y = [-3.1,3.1]
end else if keyword_set(kvesta) then begin
    x = [-2.6,2.6]
    y = [-2.6,2.6]
end else if keyword_set(kmaven) then begin
    x = [-1.7,1.7]
    y = [-1.7,1.7]
end else if keyword_set(kbepicolombo) then begin
    x = [-1.7,1.7]
    y = [-1.7,1.7]
end else if keyword_set(kpsp) then begin
    x = [-1.7,1.7]
    y = [-1.7,1.7]
end else if keyword_set(ksolo) then begin
    x = [-1.7,1.7]
    y = [-1.7,1.7]
end else begin
    if n_elements(yrange) eq 2 then x = yrange else $
      x = [sta[1],stb[1],0,0]
    if n_elements(xrange) eq 2 then y = xrange else $
      y = [sta[0],stb[0],0,earth[0]]
endelse
;
xmin=min(x,max=xmax)  &  xdelta = xmax-xmin
ymin=min(y,max=ymax)  &  ydelta = ymax-ymin
twin_plots = (ydelta gt (3*xdelta)) and (ydelta lt 1.5)
;
;  Form Earth's orbit.
;
orbit = ssc_form_orbit(utc, 'Earth', 'HEE', /au)
xearth = orbit[1,*]
yearth = orbit[0,*]
;
;  Load the color table.
;
loadct,0
grey  = !d.table_size * 3 / 4
grey2 = !d.table_size / 2
blue = 1        &  linecolor,blue,'blue'
yellow = 2      &  linecolor,yellow,'yellow'
red = 3         &  linecolor,red,'red'
orange = 4      &  linecolor,orange,'orange'
purple = 5      &  linecolor,purple,'purple'
green = 6       &  linecolor,green,'green'
;
;  If a white background is requested, then manipulate the !P.COLOR and
;  !P.BACKGROUND environment variables.
;
background_set = 0
if keyword_set(white_bg) and (!p.background eq 0) then begin
    background_set = 1
    background_save = !p.background
    color_save = !p.color
    !p.color = background_save
    !p.background = color_save
endif
;
;  If a black background, lighten up the blue color.
;
if !p.background le 0.2*!d.table_size then begin
    tvlct, rr, gg, bb, /get
    rr[blue] = 0.4 * !d.table_size
    gg[blue] = 0.6 * !d.table_size
    tvlct, rr, gg, bb
;
;  Otherwise, if a white background, then darken the yellow and grey colors.
;
end else if !p.background ge 0.8*!d.table_size then begin
    tvlct, rr, gg, bb, /get
    rr[yellow] = 0.9*rr[yellow]
    gg[yellow] = 0.9*gg[yellow]
    bb[yellow] = 0.9*bb[yellow]
    tvlct, rr, gg, bb
    grey = !d.table_size / 4
endif
;
;  Define a circle symbol.
;
circle_sym, /fill
;
;  Set up the plot area and scale.
;
if twin_plots then begin
    erase
    setview, 1, 3
    setscale, xmin*1.1, xmax*1.1, ymax, ymin, /noadjust
    plot, x, y, /nodata, xticks=3, xtickv=[-.2,0,.2], xtitle='Y (HEE)', $
      ytitle='X (HEE)'
end else begin
    setscale, xmin*1.1, xmax*1.1, ymax, ymin, /noadjust
    plot, x, y, /nodata, xtitle='Y (HEE)', ytitle='X (HEE)'
endelse
;
;  If requested, overplot the Parker spiral.
;
if keyword_set(parker) then ssc_oplot_parker, /au, line=1, color=grey2, $
  /reverse, _extra=_extra
;
;  Overplot the positions and orbits of the planets.
;
orbit = ssc_form_orbit(utc, 'Mercury', 'HEE', /au)
oplot, orbit[1,*], orbit[0,*], line=2, color=grey
xx = mercury[1]
yy = mercury[0]
if (xx ge !x.crange[0]) and (xx le !x.crange[1]) and (yy le !y.crange[0]) and $
  (yy ge !y.crange[1]) then begin
    plots, xx, yy, psym=8, color=grey
    if not keyword_set(outer) then xyouts, xx, yy, ' Mercury', color=grey
endif
;
orbit = ssc_form_orbit(utc, 'Venus', 'HEE', /au)
oplot, orbit[1,*], orbit[0,*], line=2, color=grey
xx = venus[1]
yy = venus[0]
if (xx ge !x.crange[0]) and (xx le !x.crange[1]) and (yy le !y.crange[0]) and $
  (yy ge !y.crange[1]) then begin
    plots, xx, yy, psym=8, color=grey
    if not keyword_set(outer) then xyouts, xx, yy, ' Venus', color=grey
endif
;
orbit = ssc_form_orbit(utc, 'Mars', 'HEE', /au)
oplot, orbit[1,*], orbit[0,*], line=2, color=grey
xx = mars[1]
yy = mars[0]
if (xx ge !x.crange[0]) and (xx le !x.crange[1]) and (yy le !y.crange[0]) and $
  (yy ge !y.crange[1]) then begin
    plots, xx, yy, psym=8, color=grey
    if not keyword_set(outer) then xyouts, xx, yy, ' Mars', color=grey
endif
;
if keyword_set(outer) then begin
    orbit = ssc_form_orbit(utc, 'Jupiter barycenter', 'HEE', /au)
    oplot, orbit[1,*], orbit[0,*], line=2, color=grey
    xx = jupiter[1]
    yy = jupiter[0]
    if (xx ge !x.crange[0]) and (xx le !x.crange[1]) and $
      (yy le !y.crange[0]) and (yy ge !y.crange[1]) then begin
        plots, xx, yy, psym=8, color=grey
        xyouts, xx, yy, ' Jupiter', color=grey
    endif
;
    orbit = ssc_form_orbit(utc, 'Saturn barycenter', 'HEE', /au)
    oplot, orbit[1,*], orbit[0,*], line=2, color=grey
    xx = saturn[1]
    yy = saturn[0]
    if (xx ge !x.crange[0]) and (xx le !x.crange[1]) and $
      (yy le !y.crange[0]) and (yy ge !y.crange[1]) then begin
        plots, xx, yy, psym=8, color=grey
        xyouts, xx, yy, ' Saturn', color=grey
    endif
;
    orbit = ssc_form_orbit(utc, 'Uranus barycenter', 'HEE', /au)
    oplot, orbit[1,*], orbit[0,*], line=2, color=grey
    xx = uranus[1]
    yy = uranus[0]
    if (xx ge !x.crange[0]) and (xx le !x.crange[1]) and $
      (yy le !y.crange[0]) and (yy ge !y.crange[1]) then begin
        plots, xx, yy, psym=8, color=grey
        xyouts, xx, yy, ' Uranus', color=grey
    endif
endif
;
;  If requested, overplot the positions of Ulysses and Messenger.
;
if keyword_set(kulysses) then begin
    file = concat_dir( getenv('STEREO_SPICE_OTHER'), $
                       'uly_05_10_01_07_02_01.bsp')
    get_stereo_spice_range, file, tai0, tai1, scid, /tai
    cspice_furnsh, file
    tai = utc2tai(utc)
    edate = tai2utc(tai0 > tai < tai1, /ccsds)
    orbit = ssc_form_orbit(utc, 'Ulysses', 'HEE', /au, edate=edate)
    oplot, orbit[1,*], orbit[0,*], line=2, color=purple
    if (tai ge tai0) and (tai le tai1) then begin
        ulysses = get_stereo_coord(utc, /au, system='HEE', 'Ulysses')
    end else begin
        state = get_stereo_coord(edate, 'Ulysses', system='HAE')
        cspice_utc2et, edate, et
        cspice_oscelt, state, et, mu, elts
        cspice_utc2et, utc, et
        cspice_conics, elts, et, ulysses
        convert_stereo_coord, utc, ulysses, 'HAE', 'HEE'
        ulysses = ulysses / 1.4959787D8
    endelse
    cspice_unload, file
    xx = ulysses[1]
    yy = ulysses[0]
    if (xx ge !x.crange[0]) and (xx le !x.crange[1]) and $
      (yy le !y.crange[0]) and (yy ge !y.crange[1]) then begin
        plots, xx, yy, psym=8, color=purple
        xyouts, xx, yy, ' Ulysses', color=purple
    endif
endif
;
if keyword_set(kmessenger) then begin
    file = concat_dir( getenv('STEREO_SPICE_OTHER'), $
                       'msgr_20040803_20120401_od051.bsp')
    get_stereo_spice_range, file, tai0, tai1, scid, /tai
    cspice_furnsh, file
    tai = utc2tai(utc)
    edate = tai2utc(tai0 > tai < tai1, /ccsds)
    orbit = ssc_form_orbit(utc, 'Messenger', 'HEE', /au, edate=edate)
    oplot, orbit[1,*], orbit[0,*], line=2, color=orange
    if (tai ge tai0) and (tai le tai1) then begin
        messenger = get_stereo_coord(utc, /au, system='HEE', 'Messenger')
    end else begin
        state = get_stereo_coord(edate, 'Messenger', system='HAE')
        cspice_utc2et, edate, et
        cspice_oscelt, state, et, mu, elts
        cspice_utc2et, utc, et
        cspice_conics, elts, et, messenger
        convert_stereo_coord, utc, messenger, 'HAE', 'HEE'
        messenger = messenger / 1.4959787D8
    endelse
    cspice_unload, file
    xx = messenger[1]
    yy = messenger[0]
    if (xx ge !x.crange[0]) and (xx le !x.crange[1]) and $
      (yy le !y.crange[0]) and (yy ge !y.crange[1]) then begin
        plots, xx, yy, psym=8, color=orange
        xyouts, xx, yy, ' Messenger', color=orange
    endif
endif
;
;  If requested, overplot the orbit of Rosetta.  After the rendezvous with
;  Comet 67P Churyumov-Gerasimenko, use that orbit instead.
;
if keyword_set(krosetta) then begin
    file = concat_dir( getenv('STEREO_SPICE_OTHER'), 'rosetta.bsp')
    body = 'ROSETTA'
    get_stereo_spice_range, file, tai0, tai1, scid, /tai
    tai = utc2tai(utc)
    if tai gt tai1 then begin
        file = concat_dir( getenv('STEREO_SPICE_OTHER'), $
                           '67P_CHURY_GERAS_2004_2016.BSP')
        body = 'CHURYUMOV-GERASIMENKO'
    endif
    cspice_furnsh, file
    edate = tai2utc(tai0 > tai < tai1, /ccsds)
    orbit = ssc_form_orbit(utc, body, 'HEE', /au, edate=edate)
    oplot, orbit[1,*], orbit[0,*], line=2, color=purple
    if (tai ge tai0) and (tai le tai1) then begin
        rosetta = get_stereo_coord(utc, /au, system='HEE', body)
    end else begin
        state = get_stereo_coord(edate, body, system='HAE')
        cspice_utc2et, edate, et
        cspice_oscelt, state, et, mu, elts
        cspice_utc2et, utc, et
        cspice_conics, elts, et, rosetta
        convert_stereo_coord, utc, rosetta, 'HAE', 'HEE'
        rosetta = rosetta / 1.4959787D8
    endelse
    cspice_unload, file
    xx = rosetta[1]
    yy = rosetta[0]
    if (xx ge !x.crange[0]) and (xx le !x.crange[1]) and $
      (yy le !y.crange[0]) and (yy ge !y.crange[1]) then begin
        plots, xx, yy, psym=8, color=purple
        xyouts, xx, yy, ' Rosetta', color=purple
    endif
endif
;
;  If requested, overplot the position of MAVEN.  This is only done during the
;  cruise phase.
;
if keyword_set(kmaven) then begin
    file = concat_dir( getenv('STEREO_SPICE_OTHER'), $
                       'trj_c_131118-141004_p00_cpwsr2_130328.bsp')
    body = '-202'
    get_stereo_spice_range, file, tai0, tai1, scid, /tai
    tai = utc2tai(utc)
    if (tai ge tai0) and (tai le tai1) then begin
        cspice_furnsh, file
        edate = tai2utc(tai, /ccsds)
        orbit = ssc_form_orbit(utc, body, 'HEE', /au, edate=edate)
        oplot, orbit[1,*], orbit[0,*], line=2, color=orange
        maven = get_stereo_coord(utc, /au, system='HEE', body)
        cspice_unload, file
        xx = maven[1]
        yy = maven[0]
        if (xx ge !x.crange[0]) and (xx le !x.crange[1]) and $
          (yy le !y.crange[0]) and (yy ge !y.crange[1]) then begin
            plots, xx, yy, psym=8, color=orange
            xyouts, xx, yy, ' MAVEN', color=orange
        endif
    endif
endif
;
;  If requested, overplot the position of BEPICOLOMBO.  This is only done
;  during the cruise phase.
;
if keyword_set(kbepicolombo) then begin
    file = concat_dir( getenv('STEREO_SPICE_OTHER'), 'bc_mpo_fcp_latest.bsp')
    body = '-121'
    get_stereo_spice_range, file, tai0, tai1, scid, /tai
    tai = utc2tai(utc)
    tai1 = utc2tai('2024-09-05T08') ;Enter orbit with Mercury
    if (tai ge tai0) and (tai le tai1) then begin
        cspice_furnsh, file
        edate = tai2utc(tai, /ccsds)
        orbit = ssc_form_orbit(utc, body, 'HEE', /au, edate=edate)
        oplot, orbit[1,*], orbit[0,*], line=2, color=orange
        bepicolombo = get_stereo_coord(utc, /au, system='HEE', body)
        cspice_unload, file
        xx = bepicolombo[1]
        yy = bepicolombo[0]
        if (xx ge !x.crange[0]) and (xx le !x.crange[1]) and $
          (yy le !y.crange[0]) and (yy ge !y.crange[1]) then begin
            plots, xx, yy, psym=8, color=orange
            xyouts, xx, yy, ' BEPICOLOMBO', color=orange
        endif
    endif
endif
;
;  If requested, overplot the position of PSP.
;
plot_psp = keyword_set(kpsp)
catch, error_status
if error_status ne 0 then begin
    plot_psp = 0
    catch, /cancel
endif
;
if plot_psp then begin
    if utc2tai(utc) ge utc2tai('13-Aug-2018') then begin
        orbit = ssc_form_orbit(utc, 'PSP', 'HEE', /au)
        oplot, orbit[1,*], orbit[0,*], line=2, color=orange
        psp = get_sunspice_coord(utc, /au, system='HEE', 'PSP')
        xx = psp[1]
        yy = psp[0]
        if (xx ge !x.crange[0]) and (xx le !x.crange[1]) and $
          (yy le !y.crange[0]) and (yy ge !y.crange[1]) then begin
            plots, xx, yy, psym=8, color=orange
            xyouts, xx, yy, ' PSP', color=orange
        endif
    endif
endif
;
;  If requested, overplot the position of Solar Orbiter.
;
plot_solo = keyword_set(ksolo)
catch, error_status
if error_status ne 0 then begin
    plot_solo = 0
    catch, /cancel
endif
;
if plot_solo then begin
    if utc2tai(utc) ge utc2tai('13-Aug-2018') then begin
        orbit = ssc_form_orbit(utc, 'SOLO', 'HEE', /au)
        oplot, orbit[1,*], orbit[0,*], line=2, color=purple
        solo = get_sunspice_coord(utc, /au, system='HEE', 'SOLO')
        xx = solo[1]
        yy = solo[0]
        if (xx ge !x.crange[0]) and (xx le !x.crange[1]) and $
          (yy le !y.crange[0]) and (yy ge !y.crange[1]) then begin
            plots, xx, yy, psym=8, color=purple
            xyouts, xx, yy, ' SO', color=purple
        endif
    endif
endif
;
;  If requested, overplot the positions of Vesta and Ceres.
;
if keyword_set(kvesta) then begin
    file = concat_dir( getenv('STEREO_SPICE_OTHER'), 'vesta_1900_2100.bsp')
    get_stereo_spice_range, file, tai0, tai1, scid, /tai
    cspice_furnsh, file
    tai = utc2tai(utc)
    edate = tai2utc(tai0 > tai < tai1, /ccsds)
    orbit = ssc_form_orbit(utc, 'Vesta', 'HEE', /au, edate=edate)
    oplot, orbit[1,*], orbit[0,*], line=2, color=purple
    if (tai ge tai0) and (tai le tai1) then begin
        vesta = get_stereo_coord(utc, /au, system='HEE', 'Vesta')
    end else begin
        state = get_stereo_coord(edate, 'Vesta', system='HAE')
        cspice_utc2et, edate, et
        cspice_oscelt, state, et, mu, elts
        cspice_utc2et, utc, et
        cspice_conics, elts, et, vesta
        convert_stereo_coord, utc, vesta, 'HAE', 'HEE'
        vesta = vesta / 1.4959787D8
    endelse
    cspice_unload, file
    xx = vesta[1]
    yy = vesta[0]
    if (xx ge !x.crange[0]) and (xx le !x.crange[1]) and $
      (yy le !y.crange[0]) and (yy ge !y.crange[1]) then begin
        plots, xx, yy, psym=8, color=purple
        xyouts, xx, yy, ' Vesta', color=purple
    endif
endif
;
if keyword_set(kceres) then begin
    file = concat_dir( getenv('STEREO_SPICE_OTHER'), 'ceres_1900_2100.bsp')
    get_stereo_spice_range, file, tai0, tai1, scid, /tai
    cspice_furnsh, file
    tai = utc2tai(utc)
    edate = tai2utc(tai0 > tai < tai1, /ccsds)
    orbit = ssc_form_orbit(utc, 'Ceres', 'HEE', /au, edate=edate)
    oplot, orbit[1,*], orbit[0,*], line=2, color=orange
    if (tai ge tai0) and (tai le tai1) then begin
        ceres = get_stereo_coord(utc, /au, system='HEE', 'Ceres')
    end else begin
        state = get_stereo_coord(edate, 'Ceres', system='HAE')
        cspice_utc2et, edate, et
        cspice_oscelt, state, et, mu, elts
        cspice_utc2et, utc, et
        cspice_conics, elts, et, ceres
        convert_stereo_coord, utc, ceres, 'HAE', 'HEE'
        ceres = ceres / 1.4959787D8
    endelse
    cspice_unload, file
    xx = ceres[1]
    yy = ceres[0]
    if (xx ge !x.crange[0]) and (xx le !x.crange[1]) and $
      (yy le !y.crange[0]) and (yy ge !y.crange[1]) then begin
        plots, xx, yy, psym=8, color=orange
        xyouts, xx, yy, ' Ceres', color=orange
    endif
endif
;
;  Overplot the position and orbit of Earth.
;
oplot, xearth, yearth, line=2
oplot, [0, 0], [!y.crange[0],0], line=1
plots, symsize=2, 0, earth[0], psym=8, color=green
if not keyword_set(outer) then xyouts, 0, earth[0], '  Earth', color=green
;
;  Overplot the Sun.
;
plots, symsize=3, 0, 0, psym=8, color=yellow
if not keyword_set(outer) then xyouts, 0, 0, '  Sun', color=yellow
;
;  Overplot the positions of the two STEREO spacecraft.  Draw lines from the
;  Sun through each spacecraft.
;
oplot, [0,20*sta[1]], [0,20*sta[0]], line=1
oplot, [0,20*stb[1]], [0,20*stb[0]], line=1
plots, symsize=1, sta[1], sta[0], psym=8, color=red
if not keyword_set(outer) then xyouts, sta[1], sta[0], ' A', charsize=2, $
  color=red
plots, symsize=1, stb[1], stb[0], psym=8, color=blue
if not keyword_set(outer) then xyouts, stb[1], stb[0], ' B', charsize=2, $
  color=blue
;
;  Reset the scale.
;
setscale
;
;  Early in the mission, also zoom in on the Earth-Moon system.  Set up the
;  plot area and scale.
;
if twin_plots then begin
    setview, 1.5, 1.5
    mdist = 0.0027
    xmax = 1.2*max(abs([sta[1],stb[1],mdist]))
    ymax = 1.2*max(abs([earth[0]-sta[0],earth[0]-stb[0],mdist]))
    setscale, xmax, -xmax, -ymax, ymax, /noadjust
    plot, -x, -y, /nodata, xtickname=replicate(' ',30), xtick_get=xt, $
      xtitle='Y (GSE)', ytitle='X (GSE)', title='To Sun'
;
;  Add an arrow for Sun location.
;
    yy = !y.crange
    yy[0] = 0.1*yy[0] + 0.9*yy[1]
    arrow, 0, yy[0], 0, yy[1], color=yellow, /data
;
;  Label the plot ticks along the X axis.
;
    while n_elements(xt) gt 5 do begin
        nt = n_elements(xt)/4
        delta = 2*(xt[1]-xt[0])
        xt = (indgen(2*nt+1)-nt)*delta
    endwhile
    if ((xt[1]-xt[0]) lt 0.01) and (n_elements(xt) gt 3) then begin
        nt = n_elements(xt)/4
        delta = 2*(xt[1]-xt[0])
        xt = (indgen(2*nt+1)-nt)*delta
    endif
    axis, xaxis=0, xticks=n_elements(xt)-1, xtickv=xt
;
;  Overplot the orbit of Earth.  Draw lines from the Sun through Earth, and
;  through both spacecraft.
;
    oplot, -xearth, -yearth+earth[0], line=2
    oplot, [0,0], !y.crange, line=1
    oplot, -[0,2*sta[1]], -[0,2*sta[0]]+earth[0], line=1
    oplot, -[0,2*stb[1]], -[0,2*stb[0]]+earth[0], line=1
;
;  Overplot the orbit and position of the Moon.
;
    orbit = ssc_form_orbit(utc, 'Moon', 'HEE', /au)
    oplot, -orbit[1,*], -orbit[0,*]+earth[0], line=1, color=grey
    plots, symsize=1, -moon[1], -moon[0]+earth[0], psym=8, color=grey
    if !x.crange[1] lt 0.03 then $
      xyouts, -moon[1], -moon[0]+earth[0], ' Moon', color=grey
;
;  Overplot, the position of Earth.
;
    plots, symsize=3, 0, 0, psym=8, color=green
    xyouts, 0, 0, '   Earth', color=green
;
;  Get A&B positions for 7 days before and after.
;
    tai = utc2tai(utc)
    delta = 7 * 86400.d0
    tai0 = (tai - delta) > utc2tai('2006-10-26T05:01')
    tai1 =  tai + delta
    tai = tai0 + (tai1-tai0)*dindgen(1001)/1000
    orbit = get_stereo_coord(tai, 'ahead', /au, system='GEI', /novelocity)
    convert_stereo_coord, utc, orbit, 'GEI', 'GSE'
    oplot, orbit[1,*], orbit[0,*], color=red
;
    orbit = get_stereo_coord(tai, 'behind', /au, system='GEI', /novelocity)
    convert_stereo_coord, utc, orbit, 'GEI', 'GSE'
    oplot, orbit[1,*], orbit[0,*], color=blue
;
;  Overplot the positions of STEREO A&B.
;
    plots, symsize=1, -sta[1], -sta[0]+earth[0], psym=8, color=red
    xyouts, -sta[1], -sta[0]+earth[0], ' A', charsize=2, color=red
    plots, symsize=1, -stb[1], -stb[0]+earth[0], psym=8, color=blue
    xyouts, -stb[1], -stb[0]+earth[0], ' B', charsize=2, color=blue
;
;  Reset the scale and view.
;
    setscale
    setview
endif
;
if background_set then begin
    !p.color = color_save
    !p.background = background_save
endif
;
end

