;produce .txt input file for elevo.pro from ElCon/DBM-fitting output

PRO elevo_input, sc, hwidth, aspect_ratio, phi, tinit, rinit, vinit, swspeed, drag_parameter, dir, realtime=realtime

AU=149597871 ;km
r_sun=6.957d5; km


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
    if sc eq 'B' then begin
      sep=abs(pos_E[1]-pos_B[1])/!dtor
      dir_E=sep-phi
    endif
    if sc eq 'A' then begin
      sep=abs(pos_E[1]-pos_A[1])/!dtor
      dir_E=-(sep-phi)
    endif
endif

print, dir_E

rinit=strtrim(string(rinit),2)
vinit=strtrim(string(vinit),2)
swspeed=strtrim(string(swspeed),2)
drag=strtrim(string(drag_parameter*1e7), 2)
title=strmid(string(anytim(tinit,/vms)),0, 17)
epsfile=strmid(tinit,0, 11)+'_elevo.eps'
direction=strtrim(string(float(dir_E)),2)

d=strarr(30)

d[0]='Ellipse parameters'
d[1]='halfwidth:'
d[2]=hwidth
d[3]='aspectratio:'
d[4]=aspect_ratio
d[5]='direction from Earth in HEE longitude:'
d[6]=direction
d[7]='--------------------------------'
d[8]='Drag parameters'
d[9]='R0 initial distance:'
d[10]=rinit
d[11]='tinit:'
d[12]=tinit
d[13]='vinit:'
d[14]=vinit
d[15]='background wind:'
d[16]=swspeed
d[17]='gamma:'
d[18]=drag
d[19]='----------------------------------'
d[20]='plot controls'
d[21]='plot title:'
d[22]='ELEvoHI model of CME shock ('+title+'UT)'
d[23]=''
d[24]='----------------------------------'
d[25]='output control'
d[26]='eps plot filename:'
d[27]=epsfile
d[28]='steps in frames (1 frame = 10 min) in eps plot:'
d[29]='60'

;write to file

fnam=dir+'elevo_input.txt'


for o=0, n_elements(d)-1 do begin
  if o eq 0 then app=0 else app=1
  openw, lun, fnam, /get_lun, append=app
  printf, lun, d[o]
  close, lun
  free_lun, lun
endfor


end
