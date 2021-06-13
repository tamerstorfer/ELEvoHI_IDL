;+
;
; Name:       dbmfit_V1 - testversion to fix problem with negative gamma

;
; Purpose:    to fit a kinematical profile (time-distance) of a CME observed by heliospheric imagers;
;             More information can be found in Rollett et al. (2016) and Amerstorfer et al. (2018)
;
;
; Calling sequence: part of ELEvoHI package
;
;
; Authors:    Tanja Amerstorfer & Christian Moestl & Juergen Hinterreiter
;             Space Research Institute, Austrian Academy of Sciences
;			  Graz, Austria
;
; History:
;			20190201: added keyword nightly to avoid plotting of DBMfit (Juergen Hinterreiter)
;
;
;(c) 2018 T. Amerstorfer, The software is provided "as is", without warranty of any kind.
;         When using ELEvoHI for a publication, please cite Rollett et al. (2016, ApJ) and Amerstorfer et al. (2018, Space Weather).
;		  Please add in the acknowledgements section of your article, where the ELEvoHI package can be obtained (figshare doi).
;         We are happy if you could send a copy of the article to tanja.amerstorfer@oeaw.ac.at.
; -

PRO dbmfit, time, r_apex, r_error, sw, dir, runnumber, tinit, rinit, vinit, swspeed, drag_parameter, fitend, lambda, phi, startcut=startcut, endcut=endcut, silent=silent, nightly=nightly, bgsw, bgswData=bgswData, spendcut=spendcut

au=149597870.
r_sun=695700.

NaN=!Values.F_NAN

if n_elements(startcut) ne 0 then scut=startcut
if n_elements(endcut) ne 0 then ecut=endcut

;define common block for amoeba fit
common distfit, X, Y, r_init, v_init, sw_speed

undefine, rin
undefine, vin

ti = (anytim(time)-anytim(time[0]))
s = DERIV(ti, r_apex*au)
r_apex_sun=r_apex*au/r_sun

speed_errhi = DERIVSIG(ti, r_apex*au, 0.0, r_error[0,*]*au)
speed_errlo = DERIVSIG(ti, r_apex*au, 0.0, r_error[1,*]*au)

yrmax=max(s+speed_errhi)

if keyword_set(nightly) ne 1 then begin
    loadct, 0

	white_bg=1
	if (white_bg eq 1) and (!p.background eq 0) then begin
		background_save = !p.background
		color_save = !p.color
		!p.color = background_save
		!p.background = color_save
	endif else begin
		!p.background = 0
		!p.color = long(16777215)
	endelse
endif

!P.MULTI=0

;manually define starting point of fit if needed (-> cut off area where Lorentz force is still dominant)
if n_elements(startcut) eq 0 then begin
    window, 1
    plot, r_apex_sun[1:n_elements(r_apex_sun)-2], s[1:n_elements(r_apex_sun)-2], psym=1, symsize=1.2, yr=[0,yrmax], title='ElCon distance-speed profile', xtit='Heliocentric distance [Rsun]', ytit='Speed [km s!E-1!N]', charsize=1.5
    errplot, r_apex_sun[1:n_elements(r_apex_sun)-2], s[1:n_elements(r_apex_sun)-2]-speed_errlo[1:n_elements(r_apex_sun)-2], s[1:n_elements(r_apex_sun)-2]+speed_errhi[1:n_elements(r_apex_sun)-2]

    print, 'Klick on starting point of DBM fit!'
    cursor, rin, vin, /data, /down

    ;Find closest data point:
    c=where(r_apex_sun lt rin, count)

    if count ne 0 then begin
    cnear=rin-r_apex_sun[c[n_elements(c)-1]]
    endif else begin
    j=0
    scut=r_apex_sun[j]
    endelse

    c1=where(r_apex_sun gt rin, count1)

    if count1 ne 0 then begin
    c1near=r_apex_sun[c1[0]]-rin
    endif else begin
    j=n_elements(r_apex_sun)-1
    scut=r_apex_sun[j]
    endelse

    if count ne 0 and count1 ne 0 then begin
    if cnear lt c1near then begin
     scut=r_apex[c[n_elements(c)-1]]*au/r_sun
     j=c[n_elements(c)-1]
    endif
    if c1near lt cnear then begin
      scut=r_apex[c1[0]]*au/r_sun
      j=c1[0]
    endif
    endif

    a = findgen(17) * (!pi*2/16.)
    usersym, cos(A), sin(A)

    oplot, [r_apex_sun[j],r_apex_sun[j]], [s[j],s[j]], psym=8, symsize=3

    cut=j

    r_init=scut*r_sun
    v_init=s[cut]

endif else begin
    cut=startcut
    scut=r_apex[cut]*au/r_sun
    r_init=scut*r_sun
    v_init=s[cut]
endelse

print, 'Initial heliocentric distance for DBM fit in Rsun: ', scut
print, 'Initial speed for DBM fit in km/s: ', v_init

undefine, rout
undefine, vout

!P.MULTI=0

;manually define end point of fit if needed (-> cut off where measurements get unreliable, i.e. high acceleration or deceleration at the end of the track!)
if n_elements(endcut) eq 0 then begin
    window, 1
    ti = (anytim(time)-anytim(time[0]))
    s = DERIV(ti, r_apex*au)
    speed_errhi = DERIVSIG(ti, r_apex*au, 0.0, r_error[0,*]*au)
    speed_errlo = DERIVSIG(ti, r_apex*au, 0.0, r_error[1,*]*au)

    r_apex_sun=r_apex*au/r_sun
    plot, r_apex_sun[1:n_elements(r_apex_sun)-2], s[1:n_elements(r_apex_sun)-2], psym=1, symsize=1.2, yr=[0,yrmax], title='ElCon distance-speed profile', xtit='Heliocentric distance [Rsun]', ytit='Speed [km s!E-1!N]', charsize=1.5
    errplot, r_apex_sun[1:n_elements(r_apex_sun)-2], s[1:n_elements(r_apex_sun)-2]-speed_errlo[1:n_elements(r_apex_sun)-2], s[1:n_elements(r_apex_sun)-2]+speed_errhi[1:n_elements(r_apex_sun)-2]

    oplot, [r_apex_sun[cut],r_apex_sun[cut]], [s[cut],s[cut]], psym=8, symsize=3

    print, 'Klick on end point of DBM fit!'
    cursor, rout, vout, /data, /down

    ;find closest data point:
    c=where(r_apex_sun lt rout, count)
    if count ne 0 then begin
    cnear=rout-r_apex_sun[c[n_elements(c)-1]]
    endif else begin
    j=0
    ecut=r_apex_sun[j]
    endelse

    c1=where(r_apex_sun gt rout, count1)
    if count1 ne 0 then begin
    c1near=r_apex_sun[c1[0]]-rout
    endif else begin
    j=n_elements(r_apex_sun)-1
    ecut=r_apex_sun[j]
    endelse

    if count ne 0 and count1 ne 0 then begin
    if cnear lt c1near then begin
     ecut=r_apex[c[n_elements(c)-1]]*au/r_sun
     j=c[n_elements(c)-1]
    endif
    if c1near lt cnear then begin
      ecut=r_apex[c1[0]]*au/r_sun
      j=c1[0]
    endif
    endif

    a = findgen(17) * (!pi*2/16.)
    usersym, cos(A), sin(A)

    oplot, [r_apex_sun[j],r_apex_sun[j]], [s[j],s[j]], psym=8, symsize=3

    ecut=j

    ecuts=r_apex[ecut]*au/r_sun
    print, 'Cut-off of DBM fit in Rsun: ', ecuts
endif else begin
    ecuts=r_apex[ecut]*au/r_sun
    print, 'Cut-off of DBM fit in Rsun: ', ecuts
endelse

if runnumber eq 1 and keyword_set(nightly) eq 0 then begin
    a = findgen(17) * (!pi*2/16.)
    usersym, cos(A), sin(A)

    window, 2
    ti = (anytim(time)-anytim(time[0]))
    s = DERIV(ti, r_apex*au)
    speed_errhi = DERIVSIG(ti, r_apex*au, 0.0, r_error[0,*]*au)
    speed_errlo = DERIVSIG(ti, r_apex*au, 0.0, r_error[1,*]*au)

    r_apex_sun=r_apex*au/r_sun
    plot, r_apex_sun[1:n_elements(r_apex_sun)-2], s[1:n_elements(r_apex_sun)-2], psym=1, symsize=1.2, yr=[0,yrmax], title='ElCon distance-speed profile', xtit='Heliocentric distance [Rsun]', ytit='Speed [km s!E-1!N]', charsize=1.5
    errplot, r_apex_sun[1:n_elements(r_apex_sun)-2], s[1:n_elements(r_apex_sun)-2]-speed_errlo[1:n_elements(r_apex_sun)-2], s[1:n_elements(r_apex_sun)-2]+speed_errhi[1:n_elements(r_apex_sun)-2]

    oplot, [r_apex_sun[startcut],r_apex_sun[startcut]], [s[startcut],s[startcut]], psym=8, symsize=3
    oplot, [r_apex_sun[endcut],r_apex_sun[endcut]], [s[endcut],s[endcut]], psym=8, symsize=3

    x2jpeg, dir+'/Dbmfit_cuts.jpg'
endif

;**********************

t = (anytim(time[cut:ecut])-anytim(time[cut])) ;resolution in seconds
X = t
Y = r_apex[cut:ecut]*au

fitend=r_apex[ecut]

;calculate fit for a range of background solar wind speeds

winds = findgen(19)*25+250

if bgsw eq 'insitu' then begin
	;chose background solar wind speed from L1 (or STEREO) data

	startt=where(anytim(sw.time) le anytim(time[0]))
	endt=where(anytim(sw.time) gt anytim(time[n_elements(time)-1]))

	bg_speed=mean(sw[startt[n_elements(startt)-1]:endt[0]].vtot, /NaN)
	bg_speed_std=stddev(sw[startt[n_elements(startt)-1]:endt[0]].vtot, /NaN)

	min_bg_speed = min(sw[startt[n_elements(startt)-1]:endt[0]].vtot)
	max_bg_speed = max(sw[startt[n_elements(startt)-1]:endt[0]].vtot)

	print, 'mean:', bg_speed
	print, 'stddev:', bg_speed_std
	print, 'min:', min_bg_speed
	print, 'max:', max_bg_speed

	winds=fltarr(5)
	winds[0] = min_bg_speed
	winds[1] = min_bg_speed+(bg_speed-min_bg_speed)/2
	winds[2] = bg_speed
	winds[3] = bg_speed+(max_bg_speed-bg_speed)/2
	winds[4] = max_bg_speed
endif


if strupcase(bgsw) eq 'HUX' then begin
	; run ELEvoHI with the data from the modeled background solar wind
	datadir=getenv('DATA_DIR')
	event = strmid(dir, strpos(dir, '/', /reverse_search)-10, 11)
    bgsw_file = datadir + 'bgsw_WSA/' + event + 'vmap.txt'
    sc = strmid(event, 9, 1)
	winds = get_bgsw_hux(bgsw_file, time[cut], scut, r_apex_sun[ecut], phi, phi, lambda, sc);, /saveData, plotPath = dir)
endif

if strupcase(bgsw) eq 'HUXT' then begin
	; run ELEvoHI with the data from the modeled background solar wind
	sc = strmid(dir, strpos(dir, '/', /reverse_search)-1, 1)
    winds = get_bgsw_huxt(bgswData, time[cut], time[ecut], scut, r_apex_sun[ecut], phi, lambda, sc)
endif

if strupcase(bgsw) eq 'EUHFORIA' then begin
    sc = strmid(dir, strpos(dir, '/', /reverse_search)-1, 1)
    winds = get_bgsw_euhforia(bgswData, time[cut], scut, r_apex_sun[ecut], phi, lambda, sc)
endif

fitauall=fltarr(n_elements(winds),n_elements(y))
fitspeedall=fltarr(n_elements(winds),n_elements(y))

fitpara=fltarr(n_elements(winds),3)

for i=0, n_elements(winds)-1 do begin

	resi=NaN

	sw_speed=float(fix(winds[i]))

	print, '------------------------------------'
	print, 'Background solar wind speed [km/s]: ', winds[i], format='(A,1x, I3)'

	;chose right sign for gamma parameter
	if v_init ge sw_speed then begin
	   gasi=1
	   print, 'Initial speed larger then background solar wind speed!'
    endif
    if v_init lt sw_speed then begin
	   gasi=-1
	   print, 'Initial speed smaller then background solar wind speed!'
    endif

	;print, 'STOP 2'
	;stop

	A=[gasi*1.0e-07]

	;do the fitting

	;if gasi gt 0 then begin
	RES = 0
	RES = AMOEBA(1.0e-05,FUNCTION_NAME='fitdbm',FUNCTION_VALUE=values,P0=A,scale=[-1e-8,1e-8])

	fitres=fltarr(2)
	fitres[0]=res[0]
	fitres[1]=sw_speed

	;if n_elements(res) eq 1 then begin
	;  print, 'DBM fit failed to converge'
	;  continue
	;endif

	;undefine, fitpara

;	if signum(gasi) ne signum(res[0]) then begin
;	  print, 'Fit not valid: drag-parameter has wrong sign!'
;	  print, 'i: ', i
;	  drag_parameter=NaN
;	  ;stop
;	  ;if i lt n_elements(winds)-1 then continue else return
;	endif ;else begin
;	  ;fitpara=fltarr(n_elements(winds),3)
;	;endelse

	;calculate function values
	;if gasi gt 0 then begin
	fit = (1/fitres[0])*alog(1 + fitres[0]*(v_init - fitres[1])*X) + fitres[1]*X + r_init
	;endif else begin
	 ; fit = (-1/fitres[0])*alog(1 - fitres[0]*(v_init - fitres[1])*X) + fitres[1]*X + r_init
	;endelse

	fitauall[i,*]=fit/au

	;if gasi lt 0 and signum(fitres[0]) lt 0 then begin
	;  print, 'check DBM-fit procedure!'
	;  return
	;endif

	;mean residual per curve
	;if fit is not converging (e.g. because of improper bg solar wind speed) resi should be NaN.

	resi = mean(abs(y-fit))

    ;mean residual of last three points fitted
    ;resi = mean(abs([y[n_elements(y)-3], y[n_elements(y)-2], y[n_elements(y)-1]]-[fit[n_elements(y)-3], fit[n_elements(y)-2], fit[n_elements(y)-1]]))

	if resi gt 2*r_sun then begin
	   print, 'Mean residual greater 2R_sun'
	   resi = NaN
	endif

	;if i eq 0 then fitpara=fltarr(n_elements(winds),3)
	if signum(gasi) ne signum(res[0]) then begin
	  print, 'Fit not valid: drag-parameter has wrong sign!'
	  print, 'i: ', i
	  drag_parameter = NaN
	  fitpara[i,0] = NaN
	  fitpara[i,1] = sw_speed
	  fitpara[i,2] = resi/r_sun
	endif else begin
	  fitpara[i,0] = res
	  fitpara[i,1] = sw_speed
	  fitpara[i,2] = resi/r_sun
	endelse

	;calculate fitspeed
	fitspeed = DERIV(X, fit)
	fitspeedall[i,*]=fitspeed

	fit_au=fit/au

	if keyword_set(silent) ne 1 then begin

		loadct, 0

		if !p.background eq 0 then begin
			background_save = !p.background
			color_save = !p.color
			!p.color = background_save
			!p.background = color_save
		endif else begin
			!p.background = 0
			!p.color = long(16777215)
	    endelse

		window, 2

		!P.MULTI = [0,1,2]

		utplot, time, r_apex_sun, psym=1, timerange=[time[0],time[n_elements(time)-1]]
		uterrplot, time, r_apex_sun+r_error*au/r_sun, r_apex_sun-r_error*au/r_sun
		outplot, time[cut:*], fit_au*au/r_sun

		utplot, time[1:n_elements(time)-2], s[1:n_elements(time)-2], psym=1, yr=[0,yrmax], timerange=[time[0],time[n_elements(time)-1]]
		uterrplot, time[1:n_elements(time)-2], s[1:n_elements(time)-2]+speed_errhi[1:n_elements(time)-2], s[1:n_elements(time)-2]-speed_errlo[1:n_elements(time)-2]
		outplot, time[cut:*], fitspeed

		cleanplot, /silent

	endif

	!P.MULTI=0

	print, 'Drag parameter [1/km]: ', fitpara[i,0], format='(A,11x,E12.3)'
	print, 'Mean residual [solar radii]: ', fitpara[i,2], format='(A,8x,F4.2)'

	print, '------------------------------------'

	;plot outcome
    if keyword_set(silent) ne 1 then begin

    	if res ne -1 then begin

    		!P.MULTI=[0,1,2]

    		num=strmid(string(i),5,3)

    		drag=0
    		if gasi eq -1 then begin
    			drag=strmid(string(fitpara[i,0]),1,4)
    		endif else begin
    			drag=strmid(string(fitpara[i,0]),2,3)
    		endelse

    		dragexp=strmid(string(fitpara[i,0]),10,3)
    		bgspeed=strmid(string(fitpara[i,1]),6,3)
    		meanresi=strmid(string(fitpara[i,2]),5,4)

    		;sw=strtrim(string(fix(sw_speed)),2)

    		c=strtrim(string(cut),2)

    		set_plot, 'ps'
    		device, filename=dir+'dbmfit_'+strtrim(string(i), 2)+'.eps', /encapsulated, /color, XSIZE=30, YSIZE=20, BITS_PER_PIXEL=16

    		utplot, time, r_apex_sun, psym=1, charsize=1.4, thick=2, symsize=2, ytit='Heliocentric Distance [R!D!9n!N!X]'
    		uterrplot, time, r_apex_sun+r_error[0,*]*au/r_sun, r_apex_sun-r_error[1,*]*au/r_sun
    		outplot, time[cut:*], fit_au*au/r_sun, thick=2
    		xyouts, 0.68, 0.745, '!4c!X: '+drag+'x10!U'+dragexp+'!N'+' km!U-1!N', /norm, charsize=1.6
    		xyouts, 0.68, 0.7, 'sw speed: '+bgspeed+' km s!E-1!N', /norm, charsize=1.6
    		xyouts, 0.68, 0.655, 'mean residual: '+meanresi+' R!D!9n!N!X', /norm, charsize=1.6

    		utplot, time[1:n_elements(time)-2], s[1:n_elements(time)-2], timerange=[time[0],time[n_elements(time)-1]], psym=1, yr=[0,yrmax], charsize=1.4, thick=2, symsize=2, ytit='CME Apex Speed [km s!E-1!N]'
    		uterrplot, time[1:n_elements(time)-2], s[1:n_elements(time)-2]+speed_errhi[1:n_elements(time)-2], s[1:n_elements(time)-2]-speed_errlo[1:n_elements(time)-2], skip=2
    		outplot, time[cut:*], fitspeed, thick=2

    		DEVICE, /CLOSE
    		SET_PLOT, 'X'

    		!P.MULTI=0

    	endif
    endif

endfor

if isa(fitpara) eq 0 then begin
 print,  '****************************************************'
 print, '*For these parameters no valid DBM fit is possible!*'
 print, '*             Change initial distance!             *'
 print,  '****************************************************'
 tinit=0
 rinit=0
 vinit=0
 swspeed=0
 drag_parameter=0

 return
endif

counti=0

cu=where(fitpara eq 0., count)
if count gt 0 then fitpara[cu]=NaN

index = WHERE(fitpara[*,2] eq min(fitpara[*,2], /NaN), counti)

dragrangepos=3d-7
dragrangeneg=-3d-7 ;allows the drag parameter to be valid within a range of -3d-7 and 3d-7 1/km
dragrangemin = 0.2d-8
dragrangemax = 2d-7


bestIndex = -1
smallestResidual = 100.
for i = 0, n_elements(winds)-1 do begin
;	if (fitpara[i,2] lt smallestResidual) and (fitpara[i,0] ge dragrangeneg and fitpara[i,0] le dragrangepos) then begin
	if (fitpara[i,2] lt smallestResidual) and (abs(fitpara[i,0]) ge dragrangemin and abs(fitpara[i,0]) le dragrangemax) then begin
		smallestResidual = fitpara[i,2]
		bestIndex = i
	endif
endfor

index = bestIndex

if index eq -1 then begin
	print, 'No DBM fit possible!'
	tinit=NaN
	rinit=NaN
	vinit=NaN
	swspeed=NaN
	drag_parameter=NaN
	return
endif

counti=fix(counti)

if counti ne 0 then begin
	print, '------------best result:------------'
	print, '------------------------------------'
	print, 'Background solar wind speed [km/s]: ', fix(fitpara[index,1]), format='(A,1x, I3)'
	print, 'Drag parameter [1/km]: ', fitpara[index,0], format='(A,11x,E12.3)'
	print, 'Mean residual [solar radii]: ', fitpara[index,2], format='(A,8x,F4.2)'
	print, '------------------------------------'
	print, '------------------------------------'

	tinit=time[cut]

	gamma_=fitpara[index,0]
	drag_parameter=gamma_[0]

	sw_sp=fix(fitpara[index,1])
	solarwind_speed=sw_sp[0]
	swspeed=solarwind_speed

	mean_res_=fitpara[index,2]
	mean_residual=mean_res_[0]

	rinit=r_init/r_sun
	vinit=v_init

	dbmfile=dir+'dbmfit_results.sav'
	save, tinit, rinit, vinit, drag_parameter, solarwind_speed, mean_residual, cut, ecut, filename=dbmfile
	;print, 'Best result saved under ', dir, 'dbmfit_results.sav'

	;plot result
	if finite(fitpara[index, 0]) then begin

    	if keyword_set(nightly) ne 1 then begin
    		loadct, 0

    		drag=0
    		if gasi eq -1 then begin
    			drag=strmid(string(fitpara[index,0]),1,4)
    		endif else begin
    			drag=strmid(string(fitpara[index,0]),2,3)
    		endelse

    		dragexp=strmid(string(fitpara[index,0]),10,3)
    		bgspeed=strmid(string(fitpara[index,1]),6,3)
    		meanresi=strmid(string(fitpara[index,2]),5,4)

    		;white_bg=1
    		;if (white_bg eq 1) and (!p.background eq 0) then begin
    		if !p.background eq 0 then begin
    			background_save = !p.background
    			color_save = !p.color
    			!p.color = background_save
    			!p.background = color_save
    		endif else begin
    			 !p.background = 0
    			 !p.color = long(16777215)
    		endelse

    		!P.MULTI=[0,1,2]

    		utplot, time, r_apex_sun, psym=1, timerange=[time[0],time[n_elements(time)-1]], ytit='Heliocentric Distance [R!D!9n!N!X]', charsize=1.4
    		uterrplot, time, r_apex_sun+r_error[0,*]*au/r_sun, r_apex_sun-r_error[1,*]*au/r_sun
    		outplot, time[cut:*], fitauall[index,*]*au/r_sun
    		xyouts, 0.68, 0.745, '!4c!X: '+drag+'x10!U'+dragexp+'!N'+' km!U-1!N', /norm, charsize=1.6
    		xyouts, 0.68, 0.7, 'sw speed: '+bgspeed+' km s!E-1!N', /norm, charsize=1.6
    		xyouts, 0.68, 0.655, 'mean residual: '+meanresi+' R!D!9n!N!X', /norm, charsize=1.6

    		utplot, time[1:n_elements(time)-2], s[1:n_elements(time)-2], psym=1, yr=[0,yrmax], timerange=[time[0],time[n_elements(time)-1]], ytit='CME Apex Speed [km s!E-1!N]', charsize=1.4
    		uterrplot, time[1:n_elements(time)-2], s[1:n_elements(time)-2]+speed_errhi[1:n_elements(time)-2], s[1:n_elements(time)-2]-speed_errlo[1:n_elements(time)-2]
    		outplot, time[cut:*], fitspeedall[index,*]

    		!P.MULTI=0
    	endif
	endif
endif

spEndCut = fitspeedall[index,n_elements(y)-1]

startcut=cut
endcut=ecut

end
