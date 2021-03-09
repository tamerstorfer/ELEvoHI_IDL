;+
;
; Name:       deformable_front
;
; Purpose:    makes the front of the CME deformable by adjusting the kinematics to the ambient solar wind conditions
;			  for each ensemble member a .sav-file is stored with all the needed information
;			  at the end of elevohi, all these .sav-files are combined to one final .sav file ('frontDataAll.sav')
;
; Calling sequence: deformable_front, lambda, f, phi, kappa, tinit, fitend, swspeed, drag_parameter, defFrontStartTime, speedEndCut, sc, bgswdata, bgswTimeNum, runnumber, resdir, realtime=realtime
;
; Parameters (input):
;			  lambda: angular half width within ecliptic in degree
;			  f: inverse ellipse aspect ratio (b/a)
;			  phi: apex direction from the HI observer in degree
;			  kappa: angular half width with respect to the latitude
;			  tinit: initial time of the CME
;			  fitend: distance of the last fit from DBM-fitting
;			  swspeed: solar wind speed from DBM-fitting
;			  drag_parameter: drag parameter from DBM-fitting
;			  defFrontStartTime: start time for the deformable front
;			  speedEndCut: speed of the solar wind at fitend
;			  sc: HI observer ['A' or 'B']
;			  bgswdata: ambient solar wind speed from WSA/HUX model combination
;			  bgswTimeNum: time for which the ambient solar wind speed was created
;			  runnumber: number of run in the ensemble mode
;			  resdir: directory for the .sav-files
;			  realtime: keyword to use for realtime (meaning STA is already east of Earth when looking from Earth to Sun)
;
; History:    2021/03: created (Juergen Hinterreiter)
;
; Authors:    Tanja Amerstorfer & Christian Moestl & Juergen Hinterreiter
;             Space Research Institute, Austrian Academy of Sciences
;			  Graz, Austria
; -

pro deformable_front, lambda, f, phi, kappa, tinit, fitend, swspeed, drag_parameter, defFrontStartTime, speedEndCut, sc, bgswdata, bgswTimeNum, runnumber, resdir, realtime=realtime
	
	test = 0
	if test eq 1 then begin
		print, '!!!!!!!!!!!!!!!!!!!!!'
		print, '!!!!!!!!  TEST  !!!!!'
		print, '!!!!!!!!!!!!!!!!!!!!!'
		lambda = float(30.0)
		f = float(0.7)
		phi = float(92)
		kappa = float(32)
		tinit = '20-Mar-2010 14:30:39.103'
		fitend = float(0.61958951)
		swspeed = 481
		drag_parameter = float(-8.34943e-09)
		defFrontStartTime = anytim('21-Mar-2010 06:24:14.207')
		speedEndCut = float(361.357)
		sc = 'B'
		runnumber = 1
		resdir = '/home/jhinterreiter/ELEvoHI/PredictedEvents/20100319_B/results/'
		bgswData = bgsw_data
		
		restore, resdir + 'bgsw_Data_new.sav', /verbose
		bgswtimenum = anytim(bgswstarttime)
		
		
			
		help, lambda
		help, f
		help, phi
		help, kappa
		help, tinit
		help, fitend
		help, swspeed
		help, drag_parameter
		help, speedEndCut
		help, runnumber	
	endif
	runtimeStart = systime(/seconds)
	au=double(149597870.)
	r_sun=double(695700.)

	
	numberOfEllPoints = 201
	timeResolution = 10 ; in minutes
	runDays = 5
	numberTimeArray = runDays*24*60/timeResolution
	
	
	if finite(tinit) ne 0 and tinit ne 0 then begin
		pos_E=get_stereo_lonlat(tinit, 'Earth', system='HEE')
		pos_A=get_stereo_lonlat(tinit, 'Ahead', system='HEE')
		pos_B=get_stereo_lonlat(tinit, 'Behind', system='HEE')


		;calculate direction from Earth

		if sc eq 'A' then begin
			sep=abs(pos_E[1]-pos_A[1])/!dtor
			dir_E=sep-phi
		endif
		if sc eq 'B' then begin
			sep=abs(pos_E[1]-pos_B[1])/!dtor
			dir_E=-(sep-phi)
		endif
		
		if keyword_set(realtime) then begin
			print, 'KW realtime set?'
			if sc eq 'B' then begin
				sep=abs(pos_E[1]-pos_B[1])/!dtor
				dir_E=sep-phi
			endif
			if sc eq 'A' then begin
				sep=abs(pos_E[1]-pos_A[1])/!dtor
				dir_E=-(sep-phi)
			endif
		endif
		print, 'runnumber: ', runnumber

		distAU = double(fitend) ; [AU]
		distRSun = double(fitend*au/r_sun) ; [R_sun]
		distKM = double(fitend*au) ; [km]
		areaFact = double(1);(0.25*au)*(0.25*au)
		earthangle = -1 * dir_e

		theta=atan(f^2*tan(lambda*!dtor))
		omega=sqrt(cos(theta)^2*(f^2-1)+1)
		;if this factor is set to other than 1 one can make very wide ellipses around the Sun
		;if necessary
		factor=1  ; another possible free parameter
		b=distRSun*omega*sin(lambda*!dtor)/(cos(lambda*!dtor-theta)+omega*sin(lambda*!dtor))*factor
		a=b/f
		c=distRSun-b

		xangle=sin((earthangle)*!dtor)
		yangle=cos((earthangle)*!dtor)
		ellipse_center=[xangle,yangle]*c

		pos_ang = 90 - earthangle; angle with repect to Earth
		xc = ellipse_center[0]; x center of ellipse
		yc = ellipse_center[1]; y center of ellipse
		rmin = a;
		rmax = b;


		kap = kappa*!dtor
		cLat = distRSun*tan(kap)

		csArea = a*cLat*!pi*r_sun*r_sun ; [km^2]

  
		rhoSW = double(get_bgsw_density(fitend, swspeed, mass_dens=1)) ; [g/km^3]
		dp = abs(double(drag_parameter)) ; [km^-1]
		mass = csArea*rhoSW/dp ; [g]

;		print, 'EC: ', eC
		print, 'SWspeed [km/s]: ', swspeed
		print, 'rhoSW [g/km^3]: ', rhoSW
		print, 'dp [km^-1]:', dp
		print, 'csArea [km^2]: ', csArea
		print, 'mass [g]: ', mass

		dragNew = csArea*rhoSW/mass
;		print, 'drag new: '
;		print, dragNew
  

		;to draw parts of ellipse phi should go from 0 to 135 and then from 225 to 360
		phiEll=findgen(numberOfEllPoints)+260
   
		phiEll=phiEll*!dpi/180.

		;original version
		;phi = 2*!pi*(findgen(npoints)/(npoints-1))       ;Divide circle into Npoints
		ang = pos_ang/!RADEG                            ;Position angle in radians
		cosang = cos(ang)
		sinang = sin(ang)

		x =  rmax*cos(phiEll)              ;Parameterized equation of ellipse
		y =  rmin*sin(phiEll)

		xprime = xc + x*cosang - y*sinang      ;Rotate to desired position angle
		yprime = yc + x*sinang + y*cosang


		R_ellipse = sqrt(xprime*xprime + yprime*yprime)
		lon_ell = atan(yprime, xprime)

   
		;print, r_ell
		lonNew = (lon_ell - lon_ell[n_elements(lon_ell)/2])*180/!dpi + earthangle
		;print, lonNew

		vbefore = fltarr(n_elements(r_ellipse))
		dragRuns = dblarr(n_elements(r_ellipse))
		densRuns = dblarr(n_elements(r_ellipse))

		;vbefore[*] = vinit
		vbefore[*] = speedEndcut

		rbefore = r_ellipse
		tini = defFrontStartTime
		tinitnum = tini

		test = 0
		tdrag = dblarr(numberTimeArray)
		if test eq 1 then tdrag = dindgen(15)
		tdrag[0] = tinitnum

		timeLen = n_elements(tdrag)
		ellLen = n_elements(r_ellipse)

		frontArr = dblarr(timelen, ellLen)
		vArr = dblarr(timelen, ellLen)
		dragArr = dblarr(timelen, ellLen)
		densArr = dblarr(timelen, ellLen)


		;print, 'CME mass: ', mass
		;stop
		for j = 0, n_elements(tdrag)-1 do begin
		;for j = 1, 3-1 do begin
			;tdrag[j] = j*10.*60. + anytim(tinitnum)

			; 10 for 10 minute time resolution

			tdrag[j] = (j+1)*timeResolution*60. + anytim(tini)
			if test eq 1 then tdrag[j] = (j+1)*240.*60. + anytim(tini)

			; cLat is the latitduinal extent of the CME front
			; taken is the radius of the apex to calculate the latitudinal extent
			cLat = rbefore[n_elements(rbefore)/2]*tan(kap)
			areaRun = cLat*a*!pi*r_sun*r_sun
			bgswdataTest = bgswData
			if min(rbefore) lt 250 then begin
				for i=0, n_elements(r_ellipse)-1 do begin
					sw_speed = get_bgsw_speed(bgswdata=bgswdata, bgswTimeNum, tinitnum, lonNew[i], rbefore[i])
					sw_density = double(get_bgsw_density(rbefore[i]/au*r_sun, sw_speed, mass_dens=1))

					accsign = 1

					gammaparam = areaRun*sw_density/mass
					
					dragRuns[i] = gammaparam
					densRuns[i] = sw_density

					background_wind = sw_speed

					if vbefore[i] lt background_wind then accsign = -1
		 
					rnew=(accsign/(gammaparam))*alog(1+(accsign*(gammaparam)*((vbefore[i]-background_wind)*(tdrag[j]-tinitnum))))+background_wind*(tdrag[j]-tinitnum)+rbefore[i]*r_sun

			 ;if vbefore[i] gt background_wind and accsign*gammaparam lt 0 then begin
			 ;  accsign = -1*accsign
				;print, '!!!!!!!!!!!!!!!!!!!!!!!'
				;print, '!!!accsign changed!!!!!'
			 ;endif
					vnew=(vbefore[i]-background_wind)/(1+(accsign*(gammaparam)*((vbefore[i]-background_wind)*(tdrag[j]-tinitnum))))+background_wind
					if finite(rnew) eq 0 then begin
						print, '!!!!!!!!!!!!!!!'
						print, '!!!!!!!!!!!!!!!'
						print, '!! nan value !!'
						print, '!!!!!!!!!!!!!!!'
						print, '!!!!!!!!!!!!!!!'

						print, 'i: ', i
						print, 'gamma param: ', gammaparam
						print, 'sw dens: ', sw_density
						print, 'vbefore: ', vbefore[i]
						print, 'bgsw: ', background_wind
						print, 'rbefore, ', rbefore[i]

						print, '!!!!!!!!'
						print, 'gamma original: ', drag_parameter
						print, 'dens original: ', rhoSW
						stop

						;rnew = 1
			 		endif
			 		rnew = rnew/r_sun

			 		rbefore[i] = rnew
			 		vbefore[i] = vnew			 
				endfor

				frontArr[j, *] = rbefore
				vArr[j, *] = vbefore
				dragArr[j, *] = dragRuns
				densArr[j, *] = densRuns
			endif else begin
				j = n_elements(tdrag)-1
			endelse
			tinitnum = tdrag[j]
		endfor


		timearr = tdrag
	
		noEarthHit = 0
		indMinLon = where(abs(lonnew) eq min(abs(lonnew)))
		if min(abs(lonnew)) gt 1 then noEarthHit = 1
		radsEarth = frontarr[*, indMinLon]
		;distEarth = au/r_sun
		distEarth = pos_E[0]/r_sun
		indMinRad = where(abs(radsEarth - distEarth) eq min(abs(radsEarth - distEarth)))

		distEarth = radsEarth[indMinRad]
		vEarth = vArr[indMinRad, indMinLon]
		arrtimeEarth = anytim(timearr[indMinRad], /ccsds)
		if 0 eq 1 then begin
			print, 'distEarth: ', radsEarth[indMinRad]
			print, 'Earth Arr: ', anytim(timearr[indMinRad], /ccsds)
			print, 'min Lon: ', lonnew[indMinLon]
			print, 'DistRSun: ', distRSun
			print, 'area: ', csArea
			print, 'mass: ', mass
		endif

		indEarthDirection = indMinLon
		longitude = lonnew
		timearr = anytim(timearr, /ccsds)

		if noEarthHit eq 1 then begin
			indEarthDirection = -1
			arrtimeEarth = -1
			vEarth = -1
		endif

		runArea = csArea
		runMass = mass
		runRho = rhoSW
		distMass = distRSun
		runDP = dp
		filenameFront = 'frontData_'+string(runnumber, format='(I003)')+'.sav'
		save, timearr, frontarr, vArr, dragArr, densArr, longitude, indEarthDirection, distEarth, arrtimeEarth, vEarth, phi, lambda, f, runarea, runmass, distmass, runRho, runDP, filename = resdir + filenameFront


		runtimeEnd = systime(/seconds)
	
		print, 'single run neededs ', (runtimeEnd-runtimeStart), ' seconds', format='(A20, F5.1, A9)'
 	endif
	
end