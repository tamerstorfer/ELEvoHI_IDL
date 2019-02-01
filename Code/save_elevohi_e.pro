;+
; 
; Name:       save_ELEvoHI_e.pro
; 
; Purpose:    This procedure saves the ELEvoHI ensemble results, i.e. all input and output parameters
;			  to a .txt file. 
;		      
; 
; Called by: ELEvoHI_V1.pro
;
;-

PRO save_elevohi_e, fnam, dir, pred, dt_all, new=new 

path=getenv('ELEvoHI_DIR')
data=getenv('DATA_DIR')

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

evdate=str[7]

evstr=strsplit(str[7], '/', /extract)

;check if GCS fitting is wanted 

if n_elements(evstr) eq 2 then begin
    GCS=1
    eventdate=evstr[0]
endif else begin
    GCS=0
    eventdate=evstr[0]
endelse


sourcestring=strsplit(str[10], '/', /extract)
source=sourcestring[0]
insitu=str[24]
f_in=str[14]
phi_in=str[18]
lambda_in=str[21]
if str[27] ne '' and str[29] eq '' then arr=[[str[27], str[28]]]
if str[27] ne '' and str[29] ne '' and str[31] eq '' then arr=[[str[27], str[28]],[str[29], str[30]]]
if str[27] ne '' and str[29] ne '' and str[31] ne '' then arr=[[str[27], str[28]],[str[29], str[30]], [str[31], str[32]]]

if keyword_set(new) then begin

d=strarr(41)
  
d[0]='####################################################'
d[1]='#      Results of ELEvoHI ensemble prediction      #'
d[2]='#'
d[3]='# Event date (first appearance in HI)'
d[4]='# '+evdate
d[5]='#'
d[6]='# observing spacecraft (HI)'
d[7]='# STEREO-'+sc
d[8]='#'
d[9]='# in situ spacecraft used for solar wind speed'
d[10]='# (usually most likely in situ direction)'
d[11]='# '+insitu
d[12]='#'
d[13]='# in situ arrival times'
d[14]='# MES'
d[15]='#'
d[16]='# VEX'
d[17]='#'
d[18]='# Earth'
d[19]='#'
d[20]='# STEREO-A'
d[21]='#'
d[22]='# STEREO-B' 
d[23]='#'
d[24]='#'
d[25]='# range of parameters of CME front shape'
d[26]='# phi'
d[27]='# '+strtrim(string(phi_in),2)
d[28]='# f'
d[29]='# '+strtrim(string(f_in),2)
d[30]='# lamda'
d[31]='# '+strtrim(string(lambda_in),2)
d[32]='#'
d[33]='# number of runs (total)' 
d[34]='#' 
d[35]='# number of runs without result' 
d[36]='# (no DBM fit possible)' 
d[37]='#'
d[38]='####################################################'
d[39]='phi[deg] f lambda[deg] elongation_min[deg] elongation_max[deg] startcut endcut bg_sw_speed[km/s] drag-parameter[e-7/km] t_init[UT] r_init[solar radii] v_init[km/s] mean_residual[r_sun] arr_time(MES) arr_speed(MES) arr_time(VEX) arr_speed(VEX)'
d[40]='arr_time(Earth) arr_speed(Earth) arr_time(STEREO-A) arr_speed(STEREO-A) arr_time(STEREO-B) arr_speed(STEREO-B) dt(MES) dt(VEX) dt(Earth) dt(STEREO-A) dt(STEREO-B)'



if n_elements(arr) ne 0 then begin

  j=n_elements(arr[0,*])

    for i=0, j-1 do begin

      case arr[0,i] of
        'MES': d[15]='# '+arr[1,i]
        'VEX': d[17]='# '+arr[1,i]
        'Earth': d[19]='# '+arr[1,i]
        'A': d[21]='# '+arr[1,i]	        
        'B': d[23]='# '+arr[1,i]    
        else: print, 'In situ spacecraft name not valid.'
      endcase

    endfor
  
endif  

for o=0, n_elements(d)-1 do begin
  if o eq 0 then app=0 else app=1
  openw, lun, dir+'/eELEvoHI_results.txt', /get_lun, append=app
  printf, lun, d[o]
  close, lun
  free_lun, lun
endfor

 return
 end


;add results



if source eq 'helcats' then begin
    restore, path+'PredictedEvents/'+eventdate+'/helcatsHI_track.sav'
	elon=ymean
	res=stereo_rsun(time[0],sc,distance=distance)
	elon_err=ystdd
endif else begin
    restore, data+'STEREO/HItracks/'+eventdate+'.sav'
    ;elon=track.track_y
endelse


restore, dir+'elcon_results.sav'

restore, dir+'dbmfit_results.sav'

data=strarr(28)
data[0]=string(phi, format='(I3)')
data[1]=string(f, format='(F3.1)')
data[2]=string(lambda, format='(I3)')
data[3]=string(round(elon[cut]*100)/100., format='(F5.2)')
data[4]=string(round(elon[ecut]*100)/100., format='(F6.2)')
data[5]=string(cut, format='(I2)')
data[6]=string(ecut, format='(I3)')
data[7]=strtrim(string(solarwind_speed),2)
data[8]=string(round(drag_parameter*1e9)/100., format='(F5.2)')
data[9]=strtrim(strmid(anytim(tinit, /ccsds),0,16), 2)
data[10]=trim(round(rinit*100)/100.)
data[11]=trim(round(vinit))
data[12]=mean_residual
data[13]=strtrim(strmid(pred.mes_time,0,16), 2)
data[15]=strtrim(strmid(pred.vex_time,0,16), 2)
data[17]=strtrim(strmid(pred.wind_time,0,16), 2)
data[19]=strtrim(strmid(pred.sta_time,0,16), 2)
data[21]=strtrim(strmid(pred.stb_time,0,16), 2)

if finite(pred.mes_speed) then data[14]=trim(round(pred.mes_speed)) else data[14]=strtrim(string(pred.mes_speed),2)
if finite(pred.vex_speed) then data[16]=trim(round(pred.vex_speed)) else data[16]=strtrim(string(pred.vex_speed),2)
if finite(pred.wind_speed) then data[18]=trim(round(pred.wind_speed)) else data[18]=strtrim(string(pred.wind_speed),2)
if finite(pred.sta_speed) then data[20]=trim(round(pred.sta_speed)) else data[20]=strtrim(string(pred.sta_speed),2)
if finite(pred.stb_speed) then data[22]=trim(round(pred.stb_speed)) else data[22]=strtrim(string(pred.stb_speed),2)

for l=0, 4 do begin

  if finite(dt_all[l]) then data[l+23]=trim(round(dt_all[l]*100)/100.) else data[l+23]=strtrim(string(dt_all[l]),2)

endfor


form='(I3, " ", F3.1, " ", I3, " ", F5.2, " ", F6.2, " ", I2, " ", I3, " ", I4, " ", F5.2, " ", A17, " ", F6.2, " ", I4, " ", F4.1, " ", A17, " ", A4, " ", A17, " ", A4, " ", A17, " ", A4, " ", A17, " ", A4, " ", A17, " ", A4, " ", F6.2, " ", F6.2, " ", F6.2, " ", F6.2, " ", F6.2)'


OPENW, lun, dir+'eELEvoHI_results.txt', /GET_LUN, /APPEND
PRINTF, lun, data, format=form
CLOSE, lun
FREE_LUN, lun



end