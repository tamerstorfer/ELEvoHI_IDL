;+
;
; Name:       combine_movie_files
;
; Purpose:    Creates one big .sav-file from all the individual .sav files, which can then be used to create a movie showing the CME propagation.
;
; Calling sequence: combine_movie_files, dir
;
; Parameters (input):
;			  dir: Path to the individual .sav files and also directory of the combined .sav file
;
;
; History:    2019/08: created (Juergen Hinterreiter)
;
; Authors:    Tanja Amerstorfer & Christian Moestl & Juergen Hinterreiter
;             Space Research Institute, Austrian Academy of Sciences
;			  Graz, Austria
; -

pro combine_movie_files, dir
;dir = '/home/jhinterreiter/ELEvoHI/PredictedEvents/20090623_A/ForMovie/'

	files=file_search(dir, 'formovie????.sav')

	restore, files[0]
	timesteps=n_elements(elevo_kin.all_apex_r)


	data_size=n_elements(files)*timesteps

	apex_r=dblarr(data_size)
	apex_t=strarr(data_size)
	apex_lat=fltarr(data_size)
	apex_lon=fltarr(data_size)
	apex_s=strarr(data_size)
	apex_f=fltarr(data_size)
	apex_w=fltarr(data_size)
	runnumber=intarr(data_size)
	colorflag=intarr(data_size)

	;j=1 n_elements(files)-1

	for i=0,  n_elements(files)-1 do begin
		restore, files[i]

		run=fix(strmid(files[i], strpos(files[i], '.')-4, 4))

		apex_r[i*timesteps:(round(i*timesteps+timesteps-1))]=elevo_kin.all_apex_r
		apex_t[i*timesteps:(round(i*timesteps+timesteps-1))]=elevo_kin.all_apex_t
		apex_lat[i*timesteps:(round(i*timesteps+timesteps-1))]=elevo_kin.all_apex_lat
		apex_lon[i*timesteps:(round(i*timesteps+timesteps-1))]=elevo_kin.all_apex_lon
		apex_s[i*timesteps:(round(i*timesteps+timesteps-1))]=elevo_kin.all_apex_s
		apex_f[i*timesteps:(round(i*timesteps+timesteps-1))]=elevo_kin.all_apex_f
		apex_w[i*timesteps:(round(i*timesteps+timesteps-1))]=elevo_kin.all_apex_w
		runnumber[i*timesteps:(round(i*timesteps+timesteps-1))]=intarr(timesteps)+run
		colorflag[i*timesteps:(round(i*timesteps+timesteps-1))]=intarr(timesteps)+1
	endfor


	elevo_kin = {runnumber:intarr(data_size),    $
			all_apex_r:dblarr(data_size),		$
		    all_apex_t:strarr(data_size),		$
		    all_apex_lat:fltarr(data_size),	$
	        all_apex_lon:fltarr(data_size), $
	        all_apex_s:strarr(data_size), $
	        all_apex_f:fltarr(data_size),$
	        all_apex_w:fltarr(data_size), $
	        colorflag:intarr(data_size)}


	elevo_kin.runnumber = runnumber
	elevo_kin.all_apex_r = apex_r
	elevo_kin.all_apex_t = apex_t
	elevo_kin.all_apex_lat = apex_lat
	elevo_kin.all_apex_lon = apex_lon
	elevo_kin.all_apex_s = apex_s
	elevo_kin.all_apex_f = apex_f
	elevo_kin.all_apex_w  = apex_w
	elevo_kin.colorflag = colorflag

	save, elevo_kin, filename=dir+'formovie_all_flag.sav'
	print, 'formovie_all saved'


end
