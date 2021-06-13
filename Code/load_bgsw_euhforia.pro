;+
;
; Name:       laod_bgsw_euhforia
;
; Purpose:    loads the ambient solar wind data from the EUHFORIA model for the given event
;
; Calling sequence: load_bgsw_euhforia, bgsw_file
;
; Parameters (input):
;			  bgsw_file: Path to the ambient solar wind file
;
; Parameters (output):
;			  bgswData: structure for the ambient solar wind speed and density and the given initial time
;
;
; History:    2021/05: created (Juergen Hinterreiter)
;
; Authors:    Tanja Amerstorfer & Christian Moestl & Juergen Hinterreiter
;             Space Research Institute, Austrian Academy of Sciences
;			  Graz, Austria
; -

function load_bgsw_euhforia, bgsw_file
    tinitFile = strmid(bgsw_file,0, strpos(bgsw_file, '.')) + '_tinit.txt'
    array = ''
    line = ''
    OPENR, lun, tinitfile, /GET_LUN
    WHILE NOT EOF(lun) DO BEGIN & $
    READF, lun, line & $
        array = [array, line] & $
    ENDWHILE
    ; Close the file and free the file unit
    FREE_LUN, lun
    tinit = array[1]

    r_sun=double(695700000.) ;in meter
    tinitnum = anytim(tinit)

    file_id = H5F_OPEN(bgsw_file)
    
    
    ; Open the image dataset within the file.
    ; We could also have used H5G_OPEN to open up the group first.
    dataset_lon = H5D_OPEN(file_id, 'lon_edges')
    dataset_r = H5D_OPEN(file_id, 'r_edges')
    dataset_vr = H5D_OPEN(file_id, 'vr_lat_m6')
    dataset_n = H5D_OPEN(file_id, 'n_lat_m6')

    lon_edge = H5D_READ(dataset_lon)
    r_edge = H5D_READ(dataset_r)
    vr = H5D_READ(dataset_vr)
    n = H5D_READ(dataset_n)

    H5D_CLOSE, dataset_lon
    H5D_CLOSE, dataset_r
    H5D_CLOSE, dataset_vr
    H5D_CLOSE, dataset_n
    H5F_CLOSE, file_id
        
    r1 = r_edge[1:n_elements(r_edge)-1]
    r2 = r_edge[0:n_elements(r_edge)-2]
    r_center = 0.5*(r1 + r2)/r_sun
    
    lon1 = lon_edge[1:n_elements(lon_edge)-1]
    lon2 = lon_edge[0:n_elements(lon_edge)-2]
    lon_center = fix(0.5*(lon1 + lon2)/!dtor)

    bgswTime = tinitnum
    bgswSpeed = vr
    bgswDensity = n
    
    bgswData = {tinitStr:tinit, $
                tinit:tinitnum, $
				r:r_center, $
				lon:lon_center, $
				varr:vr, $
				narr:n}
								
	return, bgswData
end

