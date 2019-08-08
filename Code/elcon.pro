; +
;
; Name:       elcon.pro
;
; Purpose:    To use after tracking a CME to convert measured HI elongation into AU
;
; Parameters:
;             elon.....elongation of the CME-track (in degrees from the Sun) (fltarr)
;             d........heliocentric distance of STEREO (remote observer) for the event (in AU)
;             Phi......seperation angle between STEREO and the CME-track (in degrees)
;             lambda...angular half width (90 deg for harmonic mean geometry)
;   	        f .......ellipse aspect ratio = b/a (1 for SSE geometry)
;
; Keywords:   none
;
; History:    March 2015
;
; Author:     Tanja Amerstorfer and Chris Moestl
;             (Space Research Institute, Austrian Academy of Sciences)
;
;(c) 2018 T. Amerstorfer, The software is provided "as is", without warranty of any kind.
;         When using ELEvoHI for a publication, please cite Rollett et al. (2016, ApJ) and Amerstorfer et al. (2018, Space Weather).
;		  Please add in the acknowledgements section of your article, where the ELEvoHI package can be obtained (figshare doi).
;         We are happy if you could send a copy of the article to tanja.amerstorfer@oeaw.ac.at.
;
; -

PRO elcon, elon, d, phi, lambda, f, R_ell

p=!dtor*phi
l=!dtor*lambda
e=!dtor*elon
R_ell = fltarr(n_elements(e))


for i=0, n_elements(R_ell)-1 do begin

  beta1 = !dpi-e[i]-p

  if beta1 gt !dpi/2. then begin
      beta1=e[i]+p
  endif



   theta=atan(f^2*tan(beta1))
   w=sqrt((cos(theta)^2)*(f^2-1)+1)

   thetas=atan(f^2*tan(l))
   ws=sqrt((cos(thetas)^2)*(f^2-1)+1)



   X=((cos(l-thetas)/sin(l))+ws)^(-1) * ((cos(e[i]+p+theta)/(w*sin(e[i]+p)))+1)


   if !dpi-e[i]-p gt !dpi/2 then begin
      X=((cos(l-thetas)/sin(l))+ws)^(-1) * ((-sin(!dpi/2+theta-e[i]-p)/(w*sin(e[i]+p)))+1)
   endif


   Y=(d*sin(e[i]))/(sin(e[i]+p))

   R_ell[i]=Y/(1-ws*X)


endfor


end





