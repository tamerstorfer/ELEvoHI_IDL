;+
;
; Name:       fitdbm
; 
; Purpose:    function to fit the ElCon outcome to forecast CME arrival times
; 
; Parameters: X
;
; Keywords:  	-
;   	    	
;    	    	
; Called by dbm_fit.pro
;
;                  
; History:    June 2015
; 
; Author:     Tanja Rollett 
;             Space Research Institute, Austrian Academy of Sciences
;-




FUNCTION fitdbm, A

common distfit, X, Y, r_init, v_init, sw_speed


fit = (1/A[0]) * alog(1 + A[0]*(v_init - sw_speed) * X) + sw_speed*X + r_init

;print, 'Gamma from fitdbm:'
;print, A[0]

residue=0
for i=0, n_elements(Y)-1 do begin
 residue=residue+abs(Y[i]-fit[i])^2
endfor
;print, residue
return, residue



END
