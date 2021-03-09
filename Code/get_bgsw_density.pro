; Name:		get_bgsw_density
;
; Purpose: 	Calculates the ambient solar wind density with respect to the radial distance and ambient solar wind speed
;			based on the Eyni & Steinitz (1980) paper: https://link.springer.com/chapter/10.1007/978-94-009-9100-2_21
;
; Calling sequence: rhoSW = double(get_bgsw_density(fitend, swspeed, mass_dens=1)) ; [g/km^3]
;
; Parameters (input):
;			r: radial distance from the sun [AU]
;			v: background solar wind speed  [km/s]
;
; Parameters (output):
;			n: ambient solar wind density in [protons/cm^3]
;				if keyword mass_dens is set, the density is given in gramm per kubic kilometers
;				n: ambient solar wind density [g/km^3]
;
; History:    2021/03: created (Juergen Hinterreiter)
;
; Authors:    Tanja Amerstorfer & Christian Moestl & Juergen Hinterreiter
;             Space Research Institute, Austrian Academy of Sciences
;			  Graz, Austria
; -
function get_bgsw_density, r, v, mass_dens=mass_dens

	r = double(r)
	v = double(v)
	n = 1.3e6 * 1/(r*r) * 1/(v*v)
	if keyword_set(mass_dens) then n=n*1.67d-24*1d15
	
	return, n

end

