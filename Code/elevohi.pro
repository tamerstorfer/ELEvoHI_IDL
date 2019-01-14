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
;							save_results...set TRUE to produce an IDL save file containing the model parameters and the prediction
;							statistics.....set TRUE to prepare results for usage of Python visualization (only possible for ensemble run)
;							silent.........set TRUE to avoid plotting of every iteration of DBMfit
;
; Required functions and procedures:  	see readme
;
; Data needed:		HI time-elongation tracks derived in the EU/FP7 HELCATS project are part of this package as a single IDL .sav-file.
;					In situ data are needed for background solar wind conditions and also have been prepared by the HELCATS team.
;					These in situ data are part of this package as well in form of IDL .sav-files for Wind, STEREO-A and STEREO-B. 
;					The HELCATS packages can be accessed at https://www.helcats-fp7.eu.
;
;
; History:    2018: ELEvoHI
;             2019/01: uploaded to github
; 
; Authors:    Tanja Amerstorfer & Christian Mšstl & JŸrgen Hinterreiter
;             Space Research Institute, Austrian Academy of Sciences
;			  Graz, Austria
;
;(c) 2018 T. Amerstorfer, The software is provided "as is", without warranty of any kind.
;         When using ELEvoHI for a publication, please cite Rollett et al. (2016, ApJ) and Amerstorfer et al. (2018, Space Weather).
;		  Please add in the acknowledgements section of your article, where the ELEvoHI package can be obtained (figshare doi, github-link).
;         We are happy if you could send a copy of the article to tanja.amerstorfer@oeaw.ac.at.
; -
PRO elevohi, save_results=save_results, statistics=statistics, silent=silent

read_config_file

path=getenv('ELEvoHI_DIR')
data=getenv('DATA_DIR')
gcs_path=getenv('EAGEL_DIR')



au=149597870.
r_sun=695700.

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
    eventdate=evstr[0]
endif else begin
    GCS=0
    eventdate=evstr[0]
endelse

;produce name for event directory
dir=path+'PredictedEvents/'+eventdate+'/'

resdir=dir+'results/'
filetest = FILE_TEST(resdir, /dir)
if filetest ne 1 then file_mkdir, resdir



gcsdone=0

datetime=eventdate+'T00:00:00'

nfiles=0

if GCS eq 1 then begin
  ;check if directory already exists
  gcsresdir=gcs_path+'results/EAGEL4ELEvoHI/'+eventdate+'/'
  filetest = FILE_TEST(gcsresdir, /dir)
  GCSFiles=['']   
  if filetest eq 1 then begin
      GCSFiles=file_search(gcsresdir, '*.sav', count=nfiles)      
  endif else begin
      spawn, 'mkdir '+gcsresdir
      gsdone=1
      print, 'Initializing GCS fitting tool...'
      EAGEL, eventdate, datetime=datetime
      parafile=gcsresdir+'EAGEL_results_'+eventdate+'*.sav'
  endelse
  
  if filetest eq 1 and nfiles eq 0 then begin
      gsdone=1
      print, 'Initializing GCS fitting tool...'
      EAGEL, eventdate, datetime=datetime
      parafile=gcsresdir+'EAGEL_results_'+eventdate+'*.sav'
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
        EAGEL, eventdate, datetime=datetime
        parafile=gcsresdir+'EAGEL_results_'+eventdate+'*.sav'
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
      gcsresdir=gcs_path+'results/EAGEL4ELEvoHI/'+eventdate+'/'
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
			EAGEL, eventdate, datetime=datetime
			parafile=gcsresdir+'EAGEL_results_'+eventdate+'*.sav'
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
			 EAGEL, eventdate, datetime=datetime
			 parafile=gcsresdir+'EAGEL_results_'+eventdate+'*.sav'      
		  endif else print, 'Aborting ELEvoHI...'
	   endelse
	endif
endif


sourcestring=strsplit(str[10], '/', /extract)
source=sourcestring[0]

;start and end of dbmfit
if n_elements(sourcestring) eq 3 then begin
  startcut=sourcestring[1]
  endcut=sourcestring[2]
endif

insitu=str[24]
if str[27] ne '' and str[29] eq '' then arr=[[str[27], str[28]]]
if str[27] ne '' and str[29] ne '' and str[31] eq '' then arr=[[str[27], str[28]],[str[29], str[30]]]
if str[27] ne '' and str[29] ne '' and str[31] ne '' then arr=[[str[27], str[28]],[str[29], str[30]], [str[31], str[32]]]

;initialize counting variable for non-converging DBMfits
nofit=0
;initialize counting variable for converging DBMfits
fitworks=0

journal, path+'logfile.log'


if source eq 'helcats' then begin
	print, 'Source file from HELCATS'
	read_hi, eventdate, sc, time, ymean, ystdd, filen, /silent
	restore, filen, /verb
	elon=ymean
	elon_err=ystdd
endif else begin
    restore, data+'STEREO/HItracks/'+eventdate+'.sav', /verb
endelse

res=stereo_rsun(time[0],sc,distance=distance)
d=distance[0]/au ; Sun-s/c distance in AU



case insitu of
	'Earth': begin
			 insitu_file=data+'HELCATS/WIND_2007to2015_HEEQ.sav'
			 restore, insitu_file, /verb
			 sw=wind
	end
	'A':     begin
			 insitu_file=data+'HELCATS/STA_2007to2015_SCEQ.sav'
			 restore, insitu_file, /verb
			 sw=sta
	end
	'B':     begin
			 insitu_file=data+'HELCATS/STB_2007to2014_SCEQ.sav'	
			 restore, insitu_file, /verb
			 sw=stb
	end
	else: print, 'In situ spacecraft not defined!'
endcase

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
     'A': angle=round(apexsta)
     'B': angle=round(apexstb)
   endcase
   
   if angle le 10 then anglestart='1' else anglestart=string(angle-10, format='(I3)')  
   angleend=string(angle+10, format='(I3)')
   phistr=[anglestart,angleend,'2']   
   
   insert_line, fnam, 18, anglestart+'/'+angleend+'/2'
   
   ;f - is fixed and not read from GCS ecliptic cut
   fstr=['0.7','1','0.1']
   
   insert_line, fnam, 14, fstr[0]+'/'+fstr[1]+'/'+fstr[2]
   
   
endelse



if n_elements(fstr) eq 3 or n_elements(phistr) eq 3 or n_elements(lambdastr) eq 3 then ensemble=1



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

endif else begin
  lambdastart=float(lambdastr[0])
  lambdaend=lambdastart
  deltalambda=0
  n_lambda=1
  lambda_arr=1

endelse

  
for k=0, n_phi-1 do begin 
  
  if n_elements(phi_arr) ne 1 then begin
    phi=phi_arr[k]
    print, phi
  endif else phi=phistart

    for l=0, n_f-1 do begin
  
      if n_elements(f_arr) ne 1 then begin
        f=f_arr[l]
        print, f
      endif else f=fstart
      
            for m=0, n_lambda-1 do begin
  
              if n_elements(lambda_arr) ne 1 then begin
                lambda=lambda_arr[m]
                print, lambda
              endif else lambda=lambdastart
  
print, '*****'
print, 'f=', f
print, 'phi=', phi
print, 'lambda=', lambda
print, '*****'


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

dbmfit, time, r_ell, r_err, sw, dir, tinit, rinit, vinit, swspeed, drag_parameter, fitend, startcut=startcut, endcut=endcut, silent=silent

print, 'Gamma after fitting:'
print, drag_parameter

;count and save number of converging and non-converging fits

a=[phi, f, lambda]
nofit_para=[0,0,0]

if tinit eq 0 then begin
    if nofit eq 0 then nofit_para=a else nofit_para=[[nofit_para],[a]]
    nofit=nofit+1
    ;if phi eq phiend then break
    ;if fstart eq fend then break
  continue
endif else begin
  fitworks=fitworks+1
endelse



elevo_input, sc, lambda, 1./f, phi, tinit, rinit, vinit, swspeed, drag_parameter, dir

elevo, dir, pred



if n_elements(arr) ne 0 then begin
  print, '------------------------------------'
  print, '*****************************************************'
  print, '*Differences of predicted and detected arrival times*'
  print, '*"-" means predicted to arrive earlier than detected*'
  print, '*****************************************************'  

  j=n_elements(arr[0,*])
  check=0

    for i=0, j-1 do begin

      case arr[0,i] of
        'MES': begin 
                 da_mes=!VALUES.F_NAN
                 if finite(pred.mes_time) then begin
        		  da_mes = (anytim(pred.mes_time) - anytim(arr[1,i]))/3600.
        		    print, '**********MESSENGER***********'
        			print, '*', round(da_mes*100)/100., 'hours', '     *', format='(A,5x,F6.2,2x,A,5x,A)'
        			print, '******************************'
        			check=1
        		 endif else print, 'No arrival predicted at MESSENGER!'
        end	
        'VEX': begin
                 da_mvex=!VALUES.F_NAN        
                 if finite(pred.vex_time) then begin
        		  da_vex = (anytim(pred.vex_time) - anytim(arr[1,i]))/3600.
        		    print, '*********Venus Express********'
        			print, '*', round(da_vex*100)/100., 'hours', '     *', format='(A,5x,F6.2,2x,A,5x,A)'
        			print, '******************************'
        			check=1
        		 endif else print, 'No arrival predicted at Venus Express!'
        end	
        'Earth': begin
                   da_earth=!VALUES.F_NAN
                   if finite(pred.wind_time) then begin
        			da_earth = (anytim(pred.wind_time) - anytim(arr[1,i]))/3600.
        			print, '************Earth*************'
        			print, '*', round(da_earth*100)/100., 'hours', '     *', format='(A,5x,F6.2,2x,A,5x,A)'
        			print, '******************************'  
        			check=1
        		   endif else print, 'No arrival predicted at Wind!' 
        end      			       
        'A': begin 
                 da_sta=!VALUES.F_NAN
                 if finite(pred.sta_time) then begin
                    da_sta = (anytim(pred.sta_time) - anytim(arr[1,i]))/3600.
        		    print, '***********STEREO-A***********'
        			print, '*', round(da_sta*100)/100., 'hours', '     *', format='(A,5x,F6.2,2x,A,5x,A)'
        			print, '******************************'
        			check=1
        	     endif else print, 'No arrival predicted at STEREO-A!'
        end	        
        'B': begin 
               da_stb=!VALUES.F_NAN
               if finite(pred.stb_time) then begin
        		  da_stb = (anytim(pred.stb_time) - anytim(arr[1,i]))/3600.
        		    print, '***********STEREO-B***********'
        			print, '*', round(da_stb*100)/100., 'hours', '     *', format='(A,5x,F6.2,2x,A,5x,A)'
        			print, '******************************'   
        			check=1
        	   endif else print, 'No arrival predicted at STEREO-B!'
        end     
        else: begin
        		  if check ne 1 then print, 'Check in situ s/c in input file!'
        end
      endcase
      
    
    endfor
  
endif  

if not isa(da_mes) then da_mes=!VALUES.F_NAN
if not isa(da_vex) then da_vex=!VALUES.F_NAN
if not isa(da_earth) then da_earth=!VALUES.F_NAN
if not isa(da_sta) then da_sta=!VALUES.F_NAN
if not isa(da_stb) then da_stb=!VALUES.F_NAN

dt_all=[da_mes,da_vex,da_earth,da_sta, da_stb]

if ensemble eq 1 and keyword_set(save_results) then begin
  save_elevohi_e, fnam, dir, pred, dt_all
endif

     
        if lambda eq lambdaend then break
        
      endfor
    
    if f eq fend then break

  endfor

  if phi eq phiend then break

endfor

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
	linenumber1=34
	data_insert1='# '+trim(fitworks+nofit)
	linenumber2=37
	data_insert2='# '+trim(nofit)

	insert_line, file, linenumber1, data_insert1
	insert_line, file, linenumber2, data_insert2

	save, nofit_para, filename=dir+'invalidFits.sav'

endif

duration_end=systime(/seconds)

print, nofit+fitworks, ' runs needed ', (duration_end-duration_start)/60., ' minutes.', format='(I4, A13, F5.1, A9)'

if ensemble ne 1 then print, 'No estimation of uncertainty in single run mode!'


journal


if keyword_set(save_results) then begin
  
  ;copy log-file in event directory
  spawn, 'cp '+path+'logfile.log ' +dir

endif




stop

end
