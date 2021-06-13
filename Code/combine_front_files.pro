;+
;
; Name:       combine_front_files
;
; Purpose:    Creates one big .sav-file from all the individual .sav files, which can then be used to create a movie showing the CME propagation.
;
; Calling sequence: combine_front_files, dir
;
; Parameters (input):
;			  dir: Path to the individual .sav files and also directory of the combined .sav file
;
;
; History:    2021/02: created (Juergen Hinterreiter)
;
; Authors:    Tanja Amerstorfer & Christian Moestl & Juergen Hinterreiter
;             Space Research Institute, Austrian Academy of Sciences
;			  Graz, Austria
; -
pro combine_front_files, dir

;	dir = '/home/jhinterreiter/ELEvoHI/PredictedEvents/20100319_B/results/'

	files=file_search(dir, 'frontData_???.sav')
	nfiles = n_elements(files)

	if nfiles gt 0 then begin
		restore, files[0]
	
		run=fix(strmid(files[0], strpos(files[0], '.')-3, 3))
		
		frontKin = {runnumber:run,    $
				timearr:timearr,		$
				frontarr:frontarr,		$
				varr:varr,		$
				dragarr:dragarr,		$
				densarr:densarr,		$
				swarr:swarr,		$
				longitude:longitude,	$
				indEarthDirection:indEarthdirection[0], $
				distEarth:distEarth[0], $
				arrtimeearth:arrtimeearth[0],$
				vEarth:vearth[0], $
				phi:phi, $
				lambda:lambda, $
				f:f, $
				area:runarea[0], $
				mass:runmass[0], $
				distmass:distmass[0], $
				rho:runRho[0], $
				dragParameter:runDP[0]}
	
		frontKins = [frontKin]
		
		for i = 1,  nfiles-1 do begin
			restore, files[i]
			run=fix(strmid(files[i], strpos(files[i], '.')-3, 3))
			frontKin = {runnumber:run,    $
				timearr:timearr,		$
				frontarr:frontarr,		$
				varr:varr,		$
				dragarr:dragarr,		$
				densarr:densarr,		$
				swarr:swarr,		$
				longitude:longitude,	$
				indEarthDirection:long(indEarthdirection[0]), $
				distEarth:distEarth[0], $
				arrtimeearth:string(arrtimeearth[0]),$
				vEarth:double(vearth[0]), $
				phi:phi, $
				lambda:lambda, $
				f:f, $
				area:runarea[0], $
				mass:runmass[0], $
				distmass:distmass[0], $
				rho:runRho[0], $
				dragParameter:runDP[0]}
				
			frontkins = [frontkins, frontkin]
		endfor
	endif

	for i = 0, nfiles-1 do begin
		file_delete, files[i]
	endfor
	
	save, frontKins, filename=dir+'frontDataAll.sav'
	print, 'front data combined'
end