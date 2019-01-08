; Name:       elevo_analytic
; 
; Purpose:    Calculating the distance from the Sun to a point along an ellipse which 
;             describes a solar transients front shape. 
;             This point is specified by the angle delta between the ellipse apex and 
;             the in situ spacecraft. The ellipse is specified by the distance of the apex, 
;             the half width in degrees in heliospheric longitude, and the aspect ratio.
;             One of the ellipse main axes is oriented along the propagation direction, the 
;             other perpendicular to it. The aspect ratio can take any values, but physical 
;             ones suitable for CMEs cluster most likely around 1.3 +/ 0.2. 
;
; Parameters: R (apex distance from Sun in AU)
;             ellipse aspect ratio a/b
;             lambda (half width of ellipse, degree heliospheric longitude)
;             delta (angle in degree between delta and in situ spacecraft, must be < lambda)
;
; Keywords:   optional:
;             speed ... set to speed of the apex, speed will be provided as output;    	             
;             out .... displays text output of ellipse parameters, default=0
; 	    	     
; Calling sequence: results=elevo_analytic(R, aspectratio, lambda, delta)       
;
; Example: d=elevo_analytic(0.5, 1.6, 80, 60)
;          d=elevo_analytic(0.5, 1.6, 80, 60, speed=500, out=1)       	
;          d=elevo_analytic(0.5, 1.6, 80, 60, speed=800, out=0)       	
;
;
; Output: the result d is the distance from Sun to the point on the ellipse front
;         in direction of delta; if the "speed" keyword is given a vector [d,v] is returned
;         with v being the speed of the point at d
;  
; Side effects: none
;                  
; History:    17 Sept 2014  version 1.0 numerical solution
;             18 Sept 2014  added keywords speed and out
;             8 Oct 2014 replaced numerical with analytic procedure
;
; Authors:    Christian Moestl, IWF/OEAW and University of Graz, Austria
;             inputs from Tanja Rollett, Pascal Demoulin and Miho Janvier
;
;-

function elevo_analytic, R, aspectratio, lambda, delta, speed=speed, out=out

f=1d/aspectratio


if abs(delta) ge lambda then begin 
  ;Half width lambda must be greater than delta!
  n=r*!Values.F_NAN
  return, n
end  

;*********************** construct ellipse



;two angles
theta=atan(f^2*tan(lambda*!dtor)) 
omega=sqrt(cos(theta)^2*(f^2-1)+1)
;value for semi-minor axis in AU
b=R*omega*sin(lambda*!dtor)/(cos(lambda*!dtor-theta)+omega*sin(lambda*!dtor))
;semi-major axis in AU
a=b/f
;center of ellipse in AU
c=R-b

;*********** if function elevo_shape is used do it like this
;res=elevo_shape(R,1/f,lambda)
;a=res[0]
;b=res[1]
;c=res[2]
;*************



;************ Speed and distance from Sun of given point along ellipse front

; analytic solution
root=sqrt(sin(delta*!dtor)^2*f^2*(b^2-c^2)+cos(delta*!dtor)^2*b^2)

dvalue_analytic_front=(c*cos(delta*!dtor)+root)/(sin(delta*!dtor)^2*f^2+cos(delta*!dtor)^2)
dvalue_analytic_rear=(c*cos(delta*!dtor)-root)/(sin(delta*!dtor)^2*f^2+cos(delta*!dtor)^2)

;taking the root positive is the front solution, the negative root is the rear solution
dvalue=dvalue_analytic_front


IF KEYWORD_SET(speed) THEN BEGIN
  deltaspeed=dvalue/R*speed
END    

IF KEYWORD_SET(out) THEN BEGIN

 print, 'a', a
 print, 'b', b
 print, 'c', c
 print, 'R', R
 print, 'f', f

 print, 'distance d in AU', dvalue
 print, 'distance of apex in AU ', R
 print, 'ratio of apex to delta point', dvalue/R

 IF KEYWORD_SET(speed) THEN BEGIN
    print, 'Apex speed  = ', speed, ' km/s'
    print, 'Delta speed = ', deltaspeed, ' km/s'
 END

END
; return dvalue and corrected speed; if speed is not set return just dvalue

IF KEYWORD_SET(speed) THEN return, [dvalue, deltaspeed] ELSE return, dvalue


end
