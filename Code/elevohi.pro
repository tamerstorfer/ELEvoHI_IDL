;+
;
; Name:       elevohi
;
; Purpose:    predicting arrival times and speeds of coronal mass ejections based on heliospheric imager observations;
;             Detailed information can be found in Rollett et al. (2016) and Amerstorfer et al. (2018)
;
; Calling sequence: elevohi, save_results=save_results, statistics=statistics, silent=silent
;
; Parameters (input):
;             ELEvoHI reads a config-file (ELEvoHI_config.dat), in which several directories are defined.
;             This file should be checked at first to make sure everything is in its place.
;             Next, the input file (input.txt) is read in.
;			  The input.txt file needs to be modified for the event of interest.
;			  Parameters required by input.txt file (please see readme-file for clear instructions):
;							eventdate......'YYYYMMDD' first measurement in HI
;							sc.............HI observer ['A' or 'B']
;							source.........'helcats' ['other']
;							phi............apex direction from the HI observer in degree
;							f..............inverse ellipse aspect ratio (b/a)
;							lambda.........angular half width within ecliptic in degree
;							insitu.........['Earth', 'A', 'B'] most likely position at 1 AU for the CME to be detected
;										   important for background solar wind conditions
;
;
;		     (output):
;			  prediction....IDL structure containing model prediction results (information on the structure: IDL> help, pred, /str)
;						    The prediction results are additionally displayed in the terminal; if set, also the validation is displayed.
;
; Keywords:
;							save_results......set TRUE to produce an IDL save file containing the model parameters and the prediction
;							statistics........set TRUE to prepare results for usage of Python visualization (only possible for ensemble run)
;							silent............set TRUE to avoid plotting of every iteration of DBMfit
;							nightly...........set TRUE to avoid plotting of the DBMfit
;							forMovie..........set TRUE to save the CME parameters (.sav-files can then be used to create a movie)
;							deformableFront...set TRUE to additionally create .sav-files with the information of the deformed CME front
;             realtime.......set TRUE to get the correct angle for STEREO-A
;             bgsw...........set to 'stat' for range of background solar wind (250 - 700 km/s (25 km/s steps))
;                            set to 'HUX', 'HUXt' for to use background solar wind from model from HUX or HUXt
;                            set to 'insitu' for in-situ sw
;
; Required functions and procedures:  	see readme
;
; Data needed:		HI time-elongation tracks derived in the EU/FP7 HELCATS project are part of this package as a single IDL .sav-file.
;					In situ data are needed for background solar wind conditions and also have been prepared by the HELCATS team.
;					These in situ data are part of this package as well in form of IDL .sav-files for Wind, STEREO-A and STEREO-B.
;					The HELCATS packages can be accessed at https://www.helcats-fp7.eu.
;
;
; History:    2016: ELEvoHI 1.0 - single-run (Rollett et al., 2016)
;			  2018: ELEvoHI 1.1 - ensemble-mode (Amerstorfer et al., 2018)
;             2018/10: ELEvoHI 1.2 - GCS ecliptic cut implemented (by J. Hinterreiter)
;             2019/01: uploaded to github
;			  2019/02: keyword 'nightly' added (Juergen Hinterreiter)
;			  2019/08: keyword 'forMovie' added
;
; Authors:    Tanja Amerstorfer & Christian Moestl & Juergen Hinterreiter
;             Space Research Institute, Austrian Academy of Sciences
;			  Graz, Austria
;
;(c) 2018 T. Amerstorfer, The software is provided "as is", without warranty of any kind.
;         When using ELEvoHI for a publication, please cite Rollett et al. (2016, ApJ) and Amerstorfer et al. (2018, Space Weather).
;		  Please add in the acknowledgements section of your article, where the ELEvoHI package can be obtained (figshare doi, github-link).
;         We are happy if you could send a copy of the article to tanja.amerstorfer@oeaw.ac.at.
; -
PRO elevohi, save_results=save_results, statistics=statistics, silent=silent, nightly=nightly, forMovie=forMovie, realtime=realtime, bgsw=bgsw, deformableFront=deformableFront


if ~keyword_set(bgsw) then bgsw = 'stat'
read_config_file

path=getenv('ELEvoHI_DIR')
data=getenv('DATA_DIR')
gcs_path=getenv('EAGEL_DIR')

au=double(149597870.)
r_sun=double(695700.)

r_start_min = 200
r_end_max = 0
phi_min = 200
phi_max = 0
lam_max = 0
eventTime = '0'

set_plot, 'x'

duration_start=systime(/seconds)

fnam=path+'Code/elevohi_input.txt'

str = STRARR(200)
OPENR, 10, fnam
    dummy = ''
    i = 0L
    WHILE NOT(EOF(10)) DO BEGIN
      READF, 10, dummy
      str[i] = dummy
      i = i + 1
    ENDWHILE
CLOSE, 10

sc=str[4]

evstr=strsplit(str[7], '/', /extract)

;check if GCS fitting is wanted

if n_elements(evstr) eq 2 then begin
    GCS=1
endif else begin
    GCS=0
endelse

eventdate=evstr[0]
eventdateSC = eventdate+'_'+sc

;produce name for event directory
dir=path+'PredictedEvents/'+eventdateSC+'/'

;file_delete, dir, /allow_n, /recursive

resdir=dir+'results/'
if FILE_TEST(resdir, /dir) ne 1 then file_mkdir, resdir

if keyword_set(forMovie) then begin
    forMovieDir = dir+'ForMovie/'
    if FILE_TEST(forMovieDir, /dir) ne 1 then file_mkdir, forMovieDir
endif

gcsdone=0

datetime=eventdate+'T00:00:00'

nfiles=0

if GCS eq 1 then begin
    ;check if directory already exists
    gcsresdir=gcs_path+'results/EAGEL4ELEvoHI/'+eventdateSC+'/'
    filetest = FILE_TEST(gcsresdir, /dir)
    GCSFiles=['']
    if filetest eq 1 then begin
        GCSFiles=file_search(gcsresdir, '*.sav', count=nfiles)
    endif else begin
        spawn, 'mkdir '+gcsresdir
        gsdone=1
        print, 'Initializing GCS fitting tool...'
        EAGEL, eventdateSC, datetime=datetime
        parafile=file_search(gcsresdir+'EAGEL_results_*.sav')
    	  if n_elements(parafile) gt 1 then begin
    		    for p=0, n_elements(parafile)-1 do begin
    				    print, p+1, ': ', parafile[p]
    		    endfor

            fnum=''
            read, fnum
            parafile = parafile[fnum-1]
        endif
    endelse

    if filetest eq 1 and nfiles eq 0 then begin
        gsdone=1
        print, 'Initializing GCS fitting tool...'
        EAGEL, eventdatesc, datetime=datetime
        parafile=file_search(gcsresdir+'EAGEL_results_*.sav')
    	  if n_elements(parafile) gt 1 then begin
      			for p=0, n_elements(parafile)-1 do begin
      				  print, p+1, ': ', parafile[p]
      			endfor

      			fnum=''
      			read, fnum
      			parafile = parafile[fnum-1]
        endif
    endif

    if GCSFiles[0] ne '' then begin
        print, 'Do you want to restore existing GCS results? Please type number of file or "n".'
        for p=0, nfiles-1 do begin
            print, p+1, ': ', GCSFiles[p]
        endfor

        fnum=''
        read, fnum

        if fnum eq 'n' then begin
            gcsdone=1
            print, 'Initializing GCS fitting tool...'
            EAGEL, eventdatesc, datetime=datetime
            parafile=file_search(gcsresdir+'EAGEL_results_*.sav')
            stop
        		if n_elements(parafile) gt 1 then begin
          			for p=0, n_elements(parafile)-1 do begin
          				  print, p+1, ': ', parafile[p]
          			endfor

          			fnum=''
          			read, fnum
          			parafile = parafile[fnum-1]
        		endif
        endif else begin
            parafile=GCSFiles[fnum-1]
            gcsdone=1
        endelse
    endif
endif

if GCS eq 0 then begin
    if str[14] eq '' or str[18] eq '' or str[21] eq '' then begin
        print, 'CME front shape parameters not properly set.'

        ;check if GCSCut-directory already exists
        gcsresdir=gcs_path+'results/EAGEL4ELEvoHI/'+eventdatesc+'/'
        filetest = FILE_TEST(gcsresdir, /dir)
        GCSFiles=['']
        if filetest eq 1 then begin
            GCSFiles=file_search(gcsresdir, '*.sav', count=nfiles)
        endif else begin
            spawn, 'mkdir '+gcsresdir
        endelse

        if GCSFiles[0] ne '' then begin
            print, 'Do you want to restore existing GCS results? Please type number of file or "n".'
            for p=0, nfiles-1 do begin
                print, p+1, ': ', GCSFiles[p]
            endfor

            fnum=''
            read, fnum

            if fnum eq 'n' then begin
                gcsdone=1
              	print, 'Initializing GCS fitting tool...'
              	EAGEL, eventdatesc, datetime=datetime
              	parafile=file_search(gcsresdir+'EAGEL_results_*.sav')
            		if n_elements(parafile) gt 1 then begin
                		for p=0, n_elements(parafile)-1 do begin
                		    print, p+1, ': ', parafile[p]
                		endfor

                		fnum=''
                		read, fnum
                		parafile = parafile[fnum-1]
                endif
            endif else begin
              	parafile=GCSFiles[fnum-1]
              	gcsdone=1
            endelse
        endif else begin
            print, 'Do you want to do GCS fitting? (y/n)'
            b=''
            read, b

            if b eq 'y' then begin
                gcsdone=1
                print, 'Initializing GCS fitting tool...'
                EAGEL, eventdatesc, datetime=datetime
                parafile=file_search(gcsresdir+'EAGEL_results_*.sav')
                if n_elements(parafile) gt 1 then begin
                		for p=0, n_elements(parafile)-1 do begin
                        print, p+1, ': ', parafile[p]
                		endfor
                		fnum=''
                		read, fnum
                		parafile = parafile[fnum-1]
                endif
            endif else print, 'Aborting ELEvoHI...'
        endelse
    endif
endif

sourcestring=strsplit(str[10], '/', /extract)
source=sourcestring[0]

;start and end of dbmfit
if n_elements(sourcestring) eq 3 then begin
    startcut=sourcestring[1]
    endcut_original=sourcestring[2]
endif

insitu=str[24]
if str[27] ne '' and str[29] eq '' then arr=[[str[27], strmid(str[28], 0, 17)]]
if str[27] ne '' and str[29] ne '' and str[31] eq '' then arr=[[str[27], strmid(str[28], 0, 17)],[str[29], strmid(str[30], 0, 17)]]
if str[27] ne '' and str[29] ne '' and str[31] ne '' then arr=[[str[27], strmid(str[28], 0, 17)],[str[29], strmid(str[30], 0, 17)], [str[31], strmid(str[32], 0, 17)]]

;initialize counting variable for non-converging DBMfits
nofit=0
nofit_para=[0,0,0]
;initialize counting variable for converging DBMfits
fitworks=0

journal, path+'logfile.log'

case source of
    'helcats': begin
        print, 'Source file from HELCATS'
        read_hi, eventdate, sc, time, elon, elon_err, filen, /save_file, /silent
        restore, filen, /verb

        track = {track_date: time, elon: elon, elon_std: elon_err, sc: sc}
        save, track, filename = dir+eventdateSC+'_ccsds.sav'
        end
    'user-defined': begin
        print, 'User-defined HI input file'
        restore, str[11], /verb
        time=track.track_date
        elon=track.elon
        elon_err=track.elon_stdd
        ;sc=track.sc

        track.track_date = anytim(track.track_date, /ccsds)
        save, track, filename = dir+eventdateSC+'_ccsds.sav'
        end
    else: print, 'Define HI input file!'
endcase

res=stereo_rsun(time[0],sc,distance=distance)
d=distance[0]/au ; Sun-s/c distance in AU

if bgsw eq 'insitu' then begin
    case insitu of
        'Earth': begin
        		 insitu_file=data+'DATACAT/WIND_2007to2018_HEEQ.sav'
        		 restore, insitu_file, /verb
        		 sw=wind
        end
      	'A':     begin
      			 insitu_file=data+'DATACAT/STA_2007to2015_SCEQ.sav'
      			 restore, insitu_file, /verb
      			 sw=sta
      	end
      	'B':     begin
      			 insitu_file=data+'DATACAT/STB_2007to2014_SCEQ.sav'
      			 restore, insitu_file, /verb
      			 sw=stb
      	end
  	else: print, 'In situ spacecraft not defined!'
  endcase
endif

ensemble=0
;check if ELEvoHI is in ensemble mode:

if gcsdone eq 0 then begin
    fstr=strsplit(str[14], '/', /extract)
    phistr=strsplit(str[18], '/', /extract)
    lambdastr=strsplit(str[21], '/', /extract)
endif else begin
    restore, parafile, /verb

    ;define parameter ranges

    ;halfangle - lambda
    lgcs=round(halfangle/10.)*10
    lstart=string(lgcs-10., format='(I2)')
    lend=string(lgcs+10., format='(I2)')
    lambdastr=[lstart,lend,'5']

    insert_line, fnam, 21, lstart+'/'+lend+'/5'

    ;phi
    case sc of
      	; only get even number for the apex direction
        'A': angle=(floor(apexsta)+1)/2*2
        'B': angle=(floor(apexstb)+1)/2*2
    endcase

    if angle le 10 then anglestart='1' else anglestart=string(angle-10, format='(I3)')

    angleend=string(angle+10, format='(I3)')
    phistr=[anglestart,angleend,'2']

    insert_line, fnam, 18, anglestart+'/'+angleend+'/2'

    ;f - is read in from input file and not used from GCS ecliptic cut
    fstr=strsplit(str[14], '/', /extract)

    insert_line, fnam, 14, fstr[0]+'/'+fstr[1]+'/'+fstr[2]
endelse

if n_elements(fstr) eq 3 or n_elements(phistr) eq 3 or n_elements(lambdastr) ge 3 then ensemble=1

if keyword_set(save_results) then begin
    ;check if directory already exists
    filetest = FILE_TEST(dir, /dir)

    ;make directory for analyzed event
    if filetest eq 0 then spawn, 'mkdir '+dir

    ;copy elevohi input file in event directory
    spawn, 'cp '+fnam+' '+dir

    if ensemble eq 1 then save_elevohi_e, fnam, dir, '0', '0', /new
endif

if ensemble eq 1 then begin
    print, '========================'
    print, 'ELEvoHI in ensemble mode'
    print, '========================'
endif

if n_elements(phistr) eq 3 then begin
    phistart=float(phistr[0])
    phiend=float(phistr[1])
    deltaphi=float(phistr[2])

    n_phi=fix((phiend-phistart)/deltaphi+1)
    phi_arr=findgen(n_phi, start=phistart, increment=deltaphi)
endif else begin
    phistart=float(phistr[0])
    phiend=phistart
    deltaphi=0
    n_phi=1
    phi_arr=1
endelse

if n_elements(fstr) eq 3 then begin
    fstart=float(fstr[0])
    fend=float(fstr[1])
    deltaf=float(fstr[2])

    n_f=fix((fend-fstart)/deltaf+1)
    f_arr=findgen(n_f, start=fstart, increment=deltaf)
endif else begin
    fstart=float(fstr[0])
    fend=fstart
    deltaf=0
    n_f=1
    f_arr=1
endelse

if n_elements(lambdastr) eq 3 then begin
    lambdastart=float(lambdastr[0])
    lambdaend=float(lambdastr[1])
    deltalambda=float(lambdastr[2])
    n_lambda=fix((lambdaend-lambdastart)/deltalambda+1)
    lambda_arr=findgen(n_lambda, start=lambdastart, increment=deltalambda)
endif

if n_elements(lambdastr) eq 1 then begin
    lambdastart=float(lambdastr[0])
    lambdaend=lambdastart
    deltalambda=0
    n_lambda=1
    lambda_arr=1
endif

kappa = -1
if n_elements(lambdastr) eq 4 then begin
    lambdastart=float(lambdastr[0])
    lambdaend=float(lambdastr[1])
    deltalambda=float(lambdastr[2])
    kappa=float(lambdastr[3])

    n_lambda=fix((lambdaend-lambdastart)/deltalambda+1)
    lambda_arr=findgen(n_lambda, start=lambdastart, increment=deltalambda)
endif

if n_elements(lambdastr) eq 2 then begin
    lambdastart=float(lambdastr[0])
    lambdaend=lambdastart
    deltalambda=0
    n_lambda=1
    lambda_arr=1
    kappa=float(lambdastr[1])
endif

;iterating runs starts here

runnumber = 0

if ensemble eq 1 then begin
    phiCenter = (phistart+phiend)/2
    lambdaCenter = (lambdastart+lambdaend)/2
    ; fCenter fixed set
    fCenter = 0.7
endif


if strupcase(bgsw) eq 'HUX' then begin
    event = strmid(dir, strpos(dir, '/', /reverse_search)-10, 11)
    bgsw_file = data + 'bgsw_WSA/' + event + 'vmap.txt'
    sc = strmid(event, 9, 1)
    
    load_bgsw_hux, bgsw_file, bgswData=bgswData, bgswTime=bgswTime
    bgsw_data = bgswData
    bgswStartTime = bgswTime
    bgswTimeNum = anytim(bgswTime)
    
    ;    bgsw_data = congrid(bgswdata, (size(bgswData))[1]*2, (size(bgswData))[2], /interp)
        ; interpolate the ambient solar wind data: longitude resolution 0.5° radial resoulution 0.5 R_sun
    ;    bgsw_data = congrid(bgswdata, (size(bgswData))[1]*4, (size(bgswData))[2]*2, /interp)

    print, 'Size bgsw data: ', size(bgsw_data)
    
    save, bgsw_data, bgswStartTime, filename = resdir + 'bgsw_Data.sav'
endif

if strupcase(bgsw) eq 'HUXT' then begin

    event = strmid(dir, strpos(dir, '/', /reverse_search)-10, 8)
    file = data + 'bgsw_HUXt/' + event + '.hdf5'
    bgswData = load_bgsw_huxt(file)
endif

if strupcase(bgsw) eq 'EUHFORIA' then begin
    event = strmid(dir, strpos(dir, '/', /reverse_search)-10, 8)
    bgswfile = data + 'bgsw_EUHFORIA/' + event + '.h5'
    
    bgswData = load_bgsw_euhforia(bgswfile)
    save, bgswdata, filename = resdir + 'bgswData_EUHFORIA.sav'
endif


for k=0, n_phi-1 do begin
    for l=0, n_f-1 do begin
        for m=0, n_lambda-1 do begin

            if n_elements(phi_arr) ne 1 then begin
                phi = phi_arr[k]
                ;print, phi
            endif else phi = phistart

            if n_elements(f_arr) ne 1 then begin
                f = f_arr[l]
                ;print, f
            endif else f = fstart

            if n_elements(lambda_arr) ne 1 then begin
                lambda=lambda_arr[m]
                ;print, lambda
            endif else lambda=lambdastart

            ;use for the first iteration the values in the middle of the parameter range
            ;change the values from the first run with those from the "middle" run
            if ensemble eq 1 then begin
                if f eq fCenter and phi eq phiCenter and lambda eq lambdaCenter then begin
                    f = fstart
                    phi = phistart
                    lambda = lambdaStart
                endif

                if k eq 0 and l eq 0 and m eq 0 then begin
                    f = fCenter
                    phi = phiCenter
                    lambda = lambdaCenter
                endif
            endif

            runnumber = runnumber + 1

            print, '*****'
            print, 'f=', f
            print, 'phi=', phi
            print, 'lambda=', lambda
            print, '*****'
            print, 'runnumber=', runnumber

            elon_err=fltarr(n_elements(elon))

            ;produce elongation measurement errors with values err_HI1=+/-0.1 deg, err_HI2=+/-0.4 deg (see Rollett et al. 2013)

            for i=0, n_elements(elon)-1 do begin
                if elon[i] lt 24. then elon_err[i]=0.1 else elon_err[i]=0.4
            endfor

            ;Here, the time-elongation track is converted in a time-distance track.
            ;An elliptical shape for the CME front is assumed, as well as constant
            ;direction of motion and a fixed (pre-defined) half width and ellipse aspect ratio.
            ;for more information please see Rollett et al. (2016, ApJ)

            ;apex heliocentric distance in AU
            elcon, elon, d, phi, lambda, f, r_ell

            ;error in AU
            elcon, elon+elon_err, d, phi, lambda, f, r_errhi
            elcon, elon-elon_err, d, phi, lambda, f, r_errlo

            r_err=fltarr(2,n_elements(r_ell))

            r_err[0,*]=r_errhi-r_ell
            r_err[1,*]=r_ell-r_errlo

            print, r_ell[2]

            save, time, r_ell, r_err, phi, lambda, f, filename=dir+'elcon_results.sav'

            ;next step is fitting the time-distance profile using the DBM
            endcut = endcut_original
            
            ec = endCut
            ;ec = (fix(startcut) + fix(endcut))/2

            print, 'SC: ', startcut
            print, 'EC: ', endCut
            print, 'EC: ', ec
            endcut = ec

            dbmfit, time, r_ell, r_err, sw, dir, runnumber, tinit, rinit, vinit, swspeed, drag_parameter, fitend, lambda, phi, startcut=startcut, endcut=ec, silent=silent, nightly=nightly, bgsw, bgswData = bgswData, spEndCut=spEndCut

            if keyword_set(deformableFront) then begin
              	if strupcase(bgsw) ne 'HUX' and strupcase(bgsw) ne 'HUXT' and strupcase(bgsw) ne 'EUHFORIA' then begin
                    print, 'Ambient solar wind must be used from model'
                    print, 'set bgsw="HUX" or "HUXt" or "EUHFORIA"'
                    stop
              	endif
              	
              	if finite(tinit) ne 0 and tinit ne 0 then begin
                    if kappa eq -1 then begin
                        print, 'Latitudinal extent of the CME not defined!!!'
                        print, 'Check elevohi_input.txt!'
                        stop
                    endif
                    print, 'Calculate deformable front: '

                    print, 'kappa: ', kappa

                    deformable_front, bgsw, lambda, f, phi, kappa, tinit, fitend, swspeed, drag_parameter, anytim(time[eC]), spEndcut, sc, bgswdata, bgswTimeNum, runnumber, resdir, realtime=realtime

                endif
            endif

            if isa(startcut) eq 1 and isa(endcut) eq 1 then begin
                startmin = r_ell[startcut]*au/r_sun
                endmax = r_ell[endcut]*au/r_sun

                if r_start_min gt startmin then r_start_min = startmin
                if r_end_max lt endmax then r_end_max = endmax
                if phi_min gt phi then phi_min = phi
                if phi_max lt phi then phi_max = phi
                if lam_max lt lambda then lam_max = lambda
                if finite(tinit) ne 0 then eventTime = tinit
            endif

            ;print, 'Gamma after fitting:'
            ;print, drag_parameter

            ;count and save number of converging and non-converging fits

            a=[phi, f, lambda]

            if drag_parameter eq 0 or finite(drag_parameter) eq 0 then begin
                if nofit eq 0 then nofit_para=a else nofit_para=[[nofit_para],[a]]
                nofit=nofit+1
                continue
            endif else begin
                fitworks=fitworks+1
            endelse

            elevo_input, sc, lambda, 1./f, phi, tinit, rinit, vinit, swspeed, drag_parameter, dir, realtime=realtime
            elevo, dir, pred, elevo_kin, runnumber

            if keyword_set(forMovie) then begin
                elevo_kin.all_apex_s = sc
                save, elevo_kin, startcut, endcut, filename=forMovieDir+'formovie'+string(runnumber, format='(I0004)')+'.sav'
            endif

            if n_elements(arr) ne 0 then begin
                print, '------------------------------------'
                print, '*****************************************************'
                print, '*Differences of predicted and detected arrival times*'
                print, '*"-" means predicted to arrive earlier than detected*'
                print, '*****************************************************'

                j = n_elements(arr[0,*])

                for i=0, j-1 do begin
                    case arr[0,i] of
                    'MES': begin
                        da_mes=!VALUES.F_NAN
                        if finite(pred.mes_time) then begin
                            da_mes = (anytim(pred.mes_time) - anytim(arr[1,i]))/3600.
                            print, '**********MESSENGER***********'
                            print, '*', round(da_mes*100)/100., 'hours', '     *', format='(A,5x,F6.2,2x,A,5x,A)'
                            print, '******************************'
                        endif else print, 'No arrival predicted at MESSENGER!'
                        end
                    'VEX': begin
                        da_mvex=!VALUES.F_NAN
                        if finite(pred.vex_time) then begin
                            da_vex = (anytim(pred.vex_time) - anytim(arr[1,i]))/3600.
                            print, '*********Venus Express********'
                            print, '*', round(da_vex*100)/100., 'hours', '     *', format='(A,5x,F6.2,2x,A,5x,A)'
                            print, '******************************'
                        endif else print, 'No arrival predicted at Venus Express!'
                        end
                    'Earth': begin
                        da_earth=!VALUES.F_NAN
                        if finite(pred.wind_time) then begin
                            da_earth = (anytim(pred.wind_time) - anytim(arr[1,i]))/3600.
                            print, '************Earth*************'
                            print, '*', round(da_earth*100)/100., 'hours', '     *', format='(A,5x,F6.2,2x,A,5x,A)'
                            print, '******************************'
                        endif else print, 'No arrival predicted at Wind!'
                        end
                    'A': begin
                        da_sta=!VALUES.F_NAN
                        if finite(pred.sta_time) then begin
                            da_sta = (anytim(pred.sta_time) - anytim(arr[1,i]))/3600.
                            print, '***********STEREO-A***********'
                            print, '*', round(da_sta*100)/100., 'hours', '     *', format='(A,5x,F6.2,2x,A,5x,A)'
                            print, '******************************'
                        endif else print, 'No arrival predicted at STEREO-A!'
                        end
                    'B': begin
                        da_stb=!VALUES.F_NAN
                        if finite(pred.stb_time) then begin
                            da_stb = (anytim(pred.stb_time) - anytim(arr[1,i]))/3600.
                            print, '***********STEREO-B***********'
                            print, '*', round(da_stb*100)/100., 'hours', '     *', format='(A,5x,F6.2,2x,A,5x,A)'
                            print, '******************************'
                        endif else print, 'No arrival predicted at STEREO-B!'
                        end
                    'SOLO': begin
                        da_solo=!VALUES.F_NAN
                        if finite(pred.solo_time) then begin
                            da_solo = (anytim(pred.solo_time) - anytim(arr[1,i]))/3600.
                            print, '***********Solar Orbiter***********'
                            print, '*', round(da_solo*100)/100., 'hours', '     *', format='(A,5x,F6.2,2x,A,5x,A)'
                            print, '******************************'
                        endif else print, 'No arrival predicted at Solar Orbiter!'
                        end
                    'PSP': begin
                        da_psp=!VALUES.F_NAN
                        if finite(pred.psp_time) then begin
                            da_psp = (anytim(pred.psp_time) - anytim(arr[1,i]))/3600.
                            print, '***********Parker Solar Probe***********'
                            print, '*', round(da_psp*100)/100., 'hours', '     *', format='(A,5x,F6.2,2x,A,5x,A)'
                            print, '******************************'
                        endif else print, 'No arrival predicted at Parker Solar Probe!'
                        end
                    'BEPI': begin
                        da_bepi=!VALUES.F_NAN
                        if finite(pred.bepi_time) then begin
                            da_bepi = (anytim(pred.bepi_time) - anytim(arr[1,i]))/3600.
                            print, '***********BEPI***********'
                            print, '*', round(da_bepi*100)/100., 'hours', '     *', format='(A,5x,F6.2,2x,A,5x,A)'
                            print, '******************************'
                        endif else print, 'No arrival predicted at BEPI!'
                        end
                        else: begin
                            print, 'Check in situ s/c in input file!'
                        end
                    endcase
                endfor
            endif

            if not isa(da_mes) then da_mes=!VALUES.F_NAN
            if not isa(da_vex) then da_vex=!VALUES.F_NAN
            if not isa(da_earth) then da_earth=!VALUES.F_NAN
            if not isa(da_sta) then da_sta=!VALUES.F_NAN
            if not isa(da_stb) then da_stb=!VALUES.F_NAN
            if not isa(da_solo) then da_solo=!VALUES.F_NAN
            if not isa(da_psp) then da_psp=!VALUES.F_NAN
            if not isa(da_bepi) then da_bepi=!VALUES.F_NAN

            dt_all=[da_mes, da_vex, da_earth, da_sta, da_stb, da_solo, da_psp, da_bepi]

            if ensemble eq 1 and keyword_set(save_results) then begin
              save_elevohi_e, fnam, dir, pred, dt_all
            endif

            if lambda eq lambdaend then break

        endfor

        if f eq fend then break

    endfor

    if phi eq phiend then break

endfor

;iteration of runs end here

if fitworks eq 0 then begin
    print, 'For this CME with the chosen parameters no prediction is possible.'
    journal
    if keyword_set(nightly) ne 1 then stop
endif

if ensemble eq 1 and keyword_set(save_results) then begin
    print, 'Ensemble results saved at '+dir+'eELEvoHI_results.txt'
endif

if keyword_set(statistics) and ensemble eq 1 then begin
    elevohi2sav, dir, path=path
    if ensemble eq 1 then begin
        print, 'Ensemble results prepared for Python statistics.'
    endif else print, 'No statistics for single run.'
endif

print, 'number of runs:'
print, nofit+fitworks
print, ' '
print, 'number of runs without result:'
print, '(no DBM fit possible)'
print, nofit

;insert number of runs into results file
;insert_line, dir, fitworks, nofit
if keyword_set(statistics) and ensemble eq 1 then begin
    file=dir+'eELEvoHI_results.txt'
    linenumber1=40
    data_insert1='# '+trim(fitworks+nofit)
    linenumber2=43
    data_insert2='# '+trim(nofit)
    insert_line, file, linenumber1, data_insert1
    insert_line, file, linenumber2, data_insert2
    save, nofit_para, filename=dir+'invalidFits.sav'
endif


if ensemble ne 1 then print, 'No estimation of uncertainty in single run mode!'

journal

if keyword_set(save_results) then begin
    ;copy log-file in event directory
    spawn, 'cp '+path+'logfile.log ' + dir
endif

if bgsw eq 'HUX' then begin
    event = strmid(dir, strpos(dir, '/', /reverse_search)-10, 11)
    bgsw_file = data + 'bgsw_WSA/' + event + 'vmap.txt'
    sc = strmid(event, 9, 1)
    wind = get_bgsw_hux(bgsw_file, eventTime, r_start_min, r_end_max, phi_min, phi_max, lam_max, sc, /savePlot, plotPath = dir, /saveData)
endif

if keyword_set(forMovie) then begin
    combine_movie_files, forMovieDir
endif

if keyword_set(deformableFront) then begin
    combine_front_files, resdir
endif

duration_end=systime(/seconds)
print, nofit+fitworks, ' runs needed ', (duration_end-duration_start)/60., ' minutes.', format='(I4, A13, F5.1, A9)'

if keyword_set(nightly) ne 1 then stop

end
