; prepare ELEvoHI ensemble results for visualization with Python


PRO elevohi2sav, dir, path=path

restore, path+'Code/ASCII_template.sav'
eELEvoHI=read_ascii(dir+'eELEvoHI_results.txt', template=temp)

;nanmes=where(eELEvoHI.arrtime_mes eq 'NaN', countmes)
;if countmes gt 0 then eELEvoHI.arrtime_mes[nanmes]=!VALUES.F_NAN

;empty output folder for current run (serves as download directory for Python visualization)
spawn, 'rm -Rf '+path+'PredictedEvents/current/*'

;read in ELEvoHI input file
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
eventdate=strmid(str[7], 0, 8)
source=str[10]
insitu=str[24]
f_in=str[14]
phi_in=str[18]
lambda_in=str[21]
arr=['','']
if str[27] ne '' and str[29] eq '' then arr=[[str[27], str[28]]]
if str[27] ne '' and str[29] ne '' and str[31] eq '' then arr=[[str[27], str[28]],[str[29], str[30]]]
if str[27] ne '' and str[29] ne '' and str[31] ne '' then arr=[[str[27], str[28]],[str[29], str[30]], [str[31], str[32]]]


;calculate transit time for each run
;transit time is time from tinit to predicted arrival time

print, '==========================='
print, 'ELEvoHI Ensemble Prediction'
print, '---------------------------'

if arr[0] ne '' then begin
	j=n_elements(arr[0,*])
endif else begin
	print, "Target must be defined in input file!"
	return
	END
	


for w=0, j-1 do begin

      ;check if directory already exists
      filetest = FILE_TEST(dir+'results/'+arr[0,w]+'/', /dir) 
      ;make directory for analyzed event
      if filetest eq 0 then file_mkdir, dir+'results/'+arr[0,w]+'/'
      

      case arr[0,w] of
        'MES': begin
                  mnan = where(finite(eELEvoHI.arrtime_mes), mcount)
                  tt = (anytim(eELEvoHI.arrtime_mes) - anytim(eELEvoHI.tinit))/3600.                  
               end
        'VEX': tt = (anytim(eELEvoHI.arrtime_vex) - anytim(eELEvoHI.tinit))/3600.
        'Earth': tt = (anytim(eELEvoHI.arrtime_earth) - anytim(eELEvoHI.tinit))/3600.
        'A': tt = (anytim(eELEvoHI.arrtime_sta) - anytim(eELEvoHI.tinit))/3600.	        
        'B': tt = (anytim(eELEvoHI.arrtime_stb) - anytim(eELEvoHI.tinit))/3600.    
        else: print, 'In situ spacecraft name not valid.'
      endcase
      
      case arr[0,w] of
        'MES': begin
                 mes_nan = where(finite(eelevohi.arrtime_mes), countmes, /null)
                 if countmes ne 0 then begin
                     tt = (anytim(eELEvoHI.arrtime_mes[mes_nan]) - anytim(eELEvoHI.tinit[mes_nan]))/3600.
				 endif
               end
        'VEX': begin
                 vex_nan = where(finite(eelevohi.arrtime_vex), countvex, /null)
                 if countvex ne 0 then begin
                     tt = (anytim(eELEvoHI.arrtime_vex[vex_nan]) - anytim(eELEvoHI.tinit[vex_nan]))/3600.
				 endif
               end
        'Earth': begin
                 earth_nan = where(finite(eelevohi.arrtime_earth), countearth, /null)
                 if countearth ne 0 then begin
                     tt = (anytim(eELEvoHI.arrtime_earth[earth_nan]) - anytim(eELEvoHI.tinit[earth_nan]))/3600.
				 endif           
                 end
        'A': begin
                 sta_nan = where(finite(eelevohi.arrtime_sta), countsta, /null)
                 if countsta ne 0 then begin
                     tt = (anytim(eELEvoHI.arrtime_sta[sta_nan]) - anytim(eELEvoHI.tinit[sta_nan]))/3600.
				 endif
             end        
        'B': begin
                 stb_nan = where(finite(eelevohi.arrtime_stb), countstb, /null)
                 if countstb ne 0 then begin
                     tt = (anytim(eELEvoHI.arrtime_stb[stb_nan]) - anytim(eELEvoHI.tinit[stb_nan]))/3600.
				 endif
             end
        else: print, 'In situ spacecraft name not valid.'
      endcase


;split tt into groups for color coding histograms in Python

	tt_d = (max(tt)-min(tt))/8.
	tt_start = min(tt)
	labels=strarr(8)

	for i=0, 7 do begin

	  bin=tt_start+(i+1)*tt_d
	  result=where(tt lt bin)
  
	  if i eq 7 then result=where(tt le bin)
  
	  case i of
		0: tt1=result
		1: tt2=result
		2: tt3=result
		3: tt4=result
		4: tt5=result
		5: tt6=result
		6: tt7=result
		7: tt8=result
	  endcase
	
	if i eq 0 then begin
	   labels[i]='$'+strtrim(string(fix(min(tt[result]))),2)+' - '+strtrim(string(fix(max(tt[result]))),2)+'$ h'
	   labelstart=strtrim(string(fix(max(tt[result]))),2)
	endif else begin
	   labels[i]='$'+labelstart+' - '+strtrim(string(fix(max(tt[result]))),2)+'$ h'
	   labelstart=strtrim(string(fix(max(tt[result]))),2)
	endelse

	endfor

    save, labels, filename=dir+'results/'+arr[0,w]+'/labels.sav'


;split other parameters into groups according to transit time

	rinit1=eELEvoHI.rinit[tt1]
	vinit1=eELEvoHI.vinit[tt1]
	gamma1=eELEvoHI.gamma[tt1]
	sw1=eELEvoHI.bg_sw_speed[tt1]

	rinit2=eELEvoHI.rinit[tt2]
	vinit2=eELEvoHI.vinit[tt2]
	gamma2=eELEvoHI.gamma[tt2]
	sw2=eELEvoHI.bg_sw_speed[tt2]

	rinit3=eELEvoHI.rinit[tt3]
	vinit3=eELEvoHI.vinit[tt3]
	gamma3=eELEvoHI.gamma[tt3]
	sw3=eELEvoHI.bg_sw_speed[tt3]

	rinit4=eELEvoHI.rinit[tt4]
	vinit4=eELEvoHI.vinit[tt4]
	gamma4=eELEvoHI.gamma[tt4]
	sw4=eELEvoHI.bg_sw_speed[tt4]

	rinit5=eELEvoHI.rinit[tt5]
	vinit5=eELEvoHI.vinit[tt5]
	gamma5=eELEvoHI.gamma[tt5]
	sw5=eELEvoHI.bg_sw_speed[tt5]

	rinit6=eELEvoHI.rinit[tt6]
	vinit6=eELEvoHI.vinit[tt6]
	gamma6=eELEvoHI.gamma[tt6]
	sw6=eELEvoHI.bg_sw_speed[tt6]

	rinit7=eELEvoHI.rinit[tt7]
	vinit7=eELEvoHI.vinit[tt7]
	gamma7=eELEvoHI.gamma[tt7]
	sw7=eELEvoHI.bg_sw_speed[tt7]

	rinit8=eELEvoHI.rinit[tt8]
	vinit8=eELEvoHI.vinit[tt8]
	gamma8=eELEvoHI.gamma[tt8]
	sw8=eELEvoHI.bg_sw_speed[tt8]

	rinit=eELEvoHI.rinit
	vinit=eELEvoHI.vinit
	gamma=eELEvoHI.gamma
	sw=eELEvoHI.bg_sw_speed

    ;check if directory already exists
    ;filetest = FILE_TEST(dir+'results/'+arr[0,w]+'/', /dir) 

    ;make directory for analyzed event
    ;if filetest eq 0 then file_mkdir, dir+'results/'+arr[0,w]+'/'

	save, rinit, rinit1, rinit2, rinit3, rinit4, rinit5, rinit6, rinit7, rinit8, filename=dir+'results/'+arr[0,w]+'/rinit.sav'
	save, vinit, vinit1, vinit2, vinit3, vinit4, vinit5, vinit6, vinit7, vinit8,  filename=dir+'results/'+arr[0,w]+'/vinit.sav'
	save, gamma, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8, filename=dir+'results/'+arr[0,w]+'/gamma.sav'
	save, sw, sw1, sw2, sw3, sw4, sw5, sw6, sw7, sw8, filename=dir+'results/'+arr[0,w]+'/sw.sav'
	save, tt, tt1, tt2, tt3, tt4, tt5, tt6, tt7, tt8, filename=dir+'results/'+arr[0,w]+'/transittimes.sav'

	residuals=eELEvoHI.mean_residual

	save, residuals, tt, filename=dir+'results/'+arr[0,w]+'/residuals.sav'
	
	save, eELEvoHI, filename=dir+'eELEvoHI_results.sav'
	
	;produce data for prediction result
	
	print, '=========================='
	
      case arr[0,w] of
        'MES': begin
                 mes_nan = where(finite(eelevohi.arrtime_mes), countmes, /null)
                 if countmes ne 0 then begin
					 arrtime_mean = anytim(mean(anytim(eELEvoHI.arrtime_mes[mes_nan])), /ccsds)
					 arrstdd = stddev(anytim(eELEvoHI.arrtime_mes[mes_nan]))/3600.
					 arrtime_median = anytim(median(anytim(eELEvoHI.arrtime_mes[mes_nan])), /ccsds)
					 arrivaltimes = eELEvoHI.arrtime_mes
					 insitu_sc='MES'
					 save, arrivaltimes, insitu_sc, filename=dir+'results/'+arr[0,w]+'/arrivaltimes.sav'
					 print, 'Arrival time at MESSENGER'
					 print, '-------------------------'
				 endif else begin
			         print, 'No arrival predicted at MESSENGER!'
				     print, '-------------------------'
				 endelse
               end
        'VEX': begin
                 vex_nan = where(finite(eelevohi.arrtime_vex), countvex, /null)
                 if countvex ne 0 then begin
					 arrtime_mean = anytim(mean(anytim(eELEvoHI.arrtime_vex[vex_nan])), /ccsds)
					 arrstdd = stddev(anytim(eELEvoHI.arrtime_vex[vex_nan]))/3600.
					 arrtime_median = anytim(median(anytim(eELEvoHI.arrtime_vex[vex_nan])), /ccsds)
					 arrivaltimes = eELEvoHI.arrtime_vex
					 insitu_sc='VEX'
					 save, arrivaltimes, insitu_sc, filename=dir+'results/'+arr[0,w]+'/arrivaltimes.sav'
					 print, 'Arrival time at Venus Express'
					 print, '-----------------------------'
				 endif else begin
				 	 print, 'No arrival predicted at Venus Express!'
				     print, '-------------------------'
				 endelse
               end
        'Earth': begin
                 earth_nan = where(finite(eelevohi.arrtime_earth), countearth, /null)
                 if countearth ne 0 then begin
					 arrtime_mean = anytim(mean(anytim(eELEvoHI.arrtime_earth[earth_nan])), /ccsds)
					 arrstdd = stddev(anytim(eELEvoHI.arrtime_earth[earth_nan]))/3600.
					 arrtime_median = anytim(median(anytim(eELEvoHI.arrtime_earth[earth_nan])), /ccsds)
					 arrivaltimes = eELEvoHI.arrtime_earth
					 insitu_sc='Earth'
					 save, arrivaltimes, insitu_sc, filename=dir+'results/'+arr[0,w]+'/arrivaltimes.sav'
					 print, 'Arrival time at Earth'
					 print, '-------------------------'     
				 endif else begin
				 	 print, 'No arrival predicted at Earth!'
				     print, '-------------------------'
				 endelse            
                 end
        'A': begin
                 sta_nan = where(finite(eelevohi.arrtime_sta), countsta, /null)
                 if countsta ne 0 then begin
					 arrtime_mean = anytim(mean(anytim(eELEvoHI.arrtime_sta[sta_nan])), /ccsds)
					 arrstdd = stddev(anytim(eELEvoHI.arrtime_sta[sta_nan]))/3600.
					 arrtime_median = anytim(median(anytim(eELEvoHI.arrtime_sta[sta_nan])), /ccsds)
					 arrivaltimes = eELEvoHI.arrtime_sta
					 insitu_sc='A'
					 save, arrivaltimes, insitu_sc, filename=dir+'results/'+arr[0,w]+'/arrivaltimes.sav'
					 print, 'Arrival time at STEREO-Ahead'
					 print, '-------------------------' 
				 endif else begin
				 	 print, 'No arrival predicted at STEREO-A!'
				     print, '-------------------------'
				 endelse
             end        
        'B': begin
                 stb_nan = where(finite(eelevohi.arrtime_stb), countstb, /null)
                 if countstb ne 0 then begin
					 arrtime_mean = anytim(mean(anytim(eELEvoHI.arrtime_stb[stb_nan]), /nan), /ccsds)
					 arrstdd = stddev(anytim(eELEvoHI.arrtime_stb[stb_nan]))/3600.
					 arrtime_median = anytim(median(anytim(eELEvoHI.arrtime_stb[stb_nan])), /ccsds)
					 arrivaltimes = eELEvoHI.arrtime_stb
					 insitu_sc='B'
					 save, arrivaltimes, insitu_sc, filename=dir+'results/'+arr[0,w]+'/arrivaltimes.sav'
					 print, 'Arrival time at STEREO-Behind'
					 print, '-------------------------' 
				 endif else begin
				 	 print, 'No arrival predicted at STEREO-B!'
				     print, '-------------------------'
				 endelse
             end
        else: print, 'In situ spacecraft name not valid.'
      endcase

 if isa(arrtime_mean) and isa(arrtime_median) then begin
   if isa(arrtime_mean) and finite(arrtime_mean) and isa(arrtime_median) and finite(arrtime_median) then begin
 
	 print, 'Ensemble mean:'
	 print, arrtime_mean, 'UT +/-', arrstdd, 'hours', format='(A16,A6,1X,F5.1,1X,A5)'
	 print, 'Ensemble median:'
	 print, arrtime_median, 'UT +/-', arrstdd, 'hours', format='(A16,A6,1X,F5.1,1X,A5)'  
	 print, '=========================='
     
     if str[28] ne '' then begin
       posSpeed = strpos(str[28], '/')
	   arrTime = str[28]
	   if posSpeed ne -1 then begin
	   	arrTime = strmid(arrTime, 0, posSpeed)
	   endif
       print, 'In situ arrival time (observed): ', arrTime+' UT', format='(A32,1X,A20)'
       print, '--------------------------'
       print, 'Difference between observation and prediction:'
       print, 'dt (mean): ', round((anytim(arrtime_mean)-anytim(arrTime))/360.)/10., 'h', format='(A11,F5.1,1X,A1)'
       print, 'dt (median): ', round((anytim(arrtime_median)-anytim(arrTime))/360.)/10., 'h', format='(A13,F5.1,1X,A1)'
       print, '=========================='
     endif
     
 
	 shadelimlow=anytim(anytim(arrtime_median)-arrstdd*3600., /ccsds)
	 shadelimhigh=anytim(anytim(arrtime_median)+arrstdd*3600., /ccsds)
 
	 save, arrtime_mean, arrtime_median, arrstdd, shadelimlow, shadelimhigh, filename=+dir+'results/'+arr[0,w]+'/prediction.sav'
 
	 cutstring1=strmid(arrtime_median, 0, 10)
	 cutstring2=strmid(arrtime_median, 11, 5)
 
	 arrplotmedian=cutstring1+' '+cutstring2
 
	 cutstring11=strmid(arrtime_mean, 0, 10)
	 cutstring21=strmid(arrtime_mean, 11, 5)
 
	 arrplotmean=cutstring11+' '+cutstring21
 
	 arrerr=string((round(arrstdd*100.)/100.), format='(F5.1)')
 
 	 insituArrSpeed = ''
	 if anytim(arr[1,w]) ne 0 then begin
	   posSpeed = strpos(arr[1,w], '/')
	   arrTime = arr[1,w]
	   if posSpeed ne -1 then begin
	   	arrTime = strmid(arrTime, 0, posSpeed)
	   	insituArrSpeed = string(strmid(arr[1,w], posSpeed+1, strlen(arr[1,w])-posSpeed+1))
	   endif
	   cutstring12=strmid(anytim(arrTime, /ccsds), 0, 10)
	   cutstring22=strmid(anytim(arrTime, /ccsds), 11, 5)
	   insituarr=cutstring12+' '+cutstring22   
	   save, arrplotmedian, arrplotmean, arrerr, insituarr, insituArrSpeed, filename=dir+'results/'+arr[0,w]+'/plottimes.sav'
	 endif else begin
	   insituarr=['']
	   save, arrplotmedian, arrplotmean, arrerr, insituarr, insituArrSpeed, filename=dir+'results/'+arr[0,w]+'/plottimes.sav'
	 endelse
  endif

  
 endif else print, '=========================='


 file_mkdir, path+'/PredictedEvents/current/'+eventdate
 
 spawn, 'cp -R '+dir+'results/'+arr[0,w]+'/ '+path+'/PredictedEvents/current/'+eventdate
 
endfor



end