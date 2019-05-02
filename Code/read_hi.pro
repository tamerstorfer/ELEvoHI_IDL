
;readin HELCATS hi tracks

PRO read_hi, eventdate, sc, time, ymean, ystdd, filen, _EXTRA = ex, silent=silent

path=getenv('ELEvoHI_DIR')
data=getenv('DATA_DIR')

fpath=data+'HELCATS/HCME_WP3_V03_TE_PROFILES'

result=file_search(fpath+'/HCME_'+sc+'__'+eventdate+'_*.dat')
restore, fpath+'/track_template.sav', /verb



if result[0] ne '' then begin

if n_elements(result) eq 1 then begin
  prof=read_ascii(result[0], template=temp)
  print, 'Single event on that day!'
endif else begin
  nu=strtrim(string(n_elements(result)),2)
  print, nu+' events on this day!'
  print, 'Please select an event!'
  print, 'Type 1-'+nu
  
  en=0
  read, en
  n=strtrim(string(en),2)
  print, 'Selected event: '+n
  
  prof=read_ascii(result[en-1], template=temp)
endelse

endif else begin
  print, 'No event on this day!'
  stop
endelse
  

;calculate mean of the different tracks

num_tracks=max(prof.track_id)
num_index=intarr(num_tracks+1)

for i=0, num_tracks do begin

  index=where(prof.track_id eq i)
  print, 'i: ', i
  print, 'number of elements:', n_elements(index)
  
  case i of
    0: begin
         tr_time0=prof.time[index]
         tr_elon0=prof.elongation[index]
         num_index[i]=n_elements(index)
       end
    1: begin
         tr_time1=prof.time[index]
         tr_elon1=prof.elongation[index]
         num_index[i]=n_elements(index)
       end
    2: begin
         tr_time2=prof.time[index]
         tr_elon2=prof.elongation[index]
         num_index[i]=n_elements(index)
       end
    3: begin
         tr_time3=prof.time[index]
         tr_elon3=prof.elongation[index]
         num_index[i]=n_elements(index)
       end
    4: begin
         tr_time4=prof.time[index]
         tr_elon4=prof.elongation[index]
         num_index[i]=n_elements(index)
       end
    5: begin
         tr_time5=prof.time[index]
         tr_elon5=prof.elongation[index]
         num_index[i]=n_elements(index)
       end
    6: begin
         tr_time6=prof.time[index]
         tr_elon6=prof.elongation[index]
         num_index[i]=n_elements(index)
       end
    else: print, 'More than 6 tracks available - only 6 are used!'
  endcase

endfor

;new time axis for all tracks:
gu=where(num_index eq max(num_index))

index=where(prof.track_id eq gu[0])

new_time=anytim(prof.time[index])
new_elon=fltarr(n_elements(new_time), num_tracks+1)

for i=0, num_tracks do begin

  case i of
    0: begin
         res = interpol(tr_elon0, anytim(tr_time0), new_time)
         new_elon[*,0]=res
       end
    1: begin
         res = interpol(tr_elon1, anytim(tr_time1), new_time)
         new_elon[*,1]=res
       end     
    2: begin
         res = interpol(tr_elon2, anytim(tr_time2), new_time)
         new_elon[*,2]=res
       end   
    3: begin
         res = interpol(tr_elon3, anytim(tr_time3), new_time)
         new_elon[*,3]=res
       end   
    4: begin
         res = interpol(tr_elon4, anytim(tr_time4), new_time)
         new_elon[*,4]=res
       end  
    5: begin
         res = interpol(tr_elon5, anytim(tr_time5), new_time)
         new_elon[*,5]=res
       end   
    6: begin
         res = interpol(tr_elon6, anytim(tr_time6), new_time)
         new_elon[*,6]=res
       end       
  endcase
endfor


if keyword_set(silent) eq 0 then begin

	utplot, anytim(new_time, /ccsds), new_elon[*,0], _EXTRA = ex

	for i=0, num_tracks do begin

	 outplot, anytim(new_time, /ccsds), new_elon[*,i], linestyle=i

	endfor

endif

	ymean=fltarr(n_elements(new_time))
	ystdd=fltarr(n_elements(new_time))

	for i=0, n_elements(new_time)-1 do begin

	  ymean[i]=mean(new_elon[i,*])
	  ystdd[i]=stddev(new_elon[i,*])

	endfor

if keyword_set(silent) eq 0 then begin

	outplot, anytim(new_time, /ccsds), ymean

endif

;filename to save mean track:
time=anytim(new_time, /ccsds)

if n_elements(result) eq 1 then u=strtrim(string(1),2) else u=n

fil=strmid(time[0],0,10)+'_0'+u
filen=path+'PredictedEvents/'+eventdate+'_'+sc+'/helcatsHI_track.sav'

save, time, ymean, ystdd, filename=filen
print, 'track saved under...'
print, filen





end