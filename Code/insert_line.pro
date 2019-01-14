
PRO insert_line, file, linenumber, data_insert



OPENR, unit, file, /GET_LUN
str = ''
count = 0ll
WHILE ~ EOF(unit) DO BEGIN
  READF, unit, str
  count = count + 1
  ENDWHILE
CLOSE, unit
FREE_LUN, unit

data=strarr(count)
OPENR, lun, file, /get_lun
READF, lun, data
CLOSE, lun
FREE_LUN, lun



data_new1=strarr(linenumber)
data_new1=data[0:linenumber-1]

data_new2=strarr(count-linenumber+1)
data_new2=data[linenumber+1:*]

new_data=[data_new1, data_insert, data_new2]

for o=0, n_elements(new_data)-1 do begin

  if o eq 0 then app=0 else app=1
OPENW, lun, file, /GET_LUN, append=app
PRINTF, lun, new_data[o]
CLOSE, lun
FREE_LUN, lun

endfor



end