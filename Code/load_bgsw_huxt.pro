;+
;
; Name:       laod_bgsw_huxt
;
; Purpose:    loads the ambient solar wind data from the HUXt model for the given event
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


function load_bgsw_huxt, file

    tinitFile = strmid(file,0, strpos(file, '.')) + '_tinit.txt'
    array = ''
    line = ''
    OPENR, lun, tinitfile, /GET_LUN
    WHILE NOT EOF(lun) DO BEGIN & $
    READF, lun, line & $
        array = [array, line] & $
    ENDWHILE
    ; Close the file and free the file unit
    FREE_LUN, lun
    tinitial = array[1]

    file_id = H5F_OPEN(file)

    ; Open the image dataset within the file.
    ; This is located within thhtml/images group.
    ; We could also have used H5G_OPEN to open up the group first.
    dataset_lon = H5D_OPEN(file_id, 'lon')
    dataset_dlon = H5D_OPEN(file_id, 'dlon')
    dataset_r = H5D_OPEN(file_id, 'r')
    dataset_cr_lon_init = H5D_OPEN(file_id, 'cr_lon_init')
    dataset_dt_out = H5D_OPEN(file_id, 'dt_out')
    dataset_v_grid = H5D_OPEN(file_id, 'v_grid')
    dataset_time = H5D_OPEN(file_id, 'time_out')

    lon = H5D_READ(dataset_lon)
    dlon = H5D_READ(dataset_dlon)
    r = H5D_READ(dataset_r)
    cr_lon_init = H5D_READ(dataset_cr_lon_init)
    dt = H5D_READ(dataset_dt_out)
    v_grid = H5D_READ(dataset_v_grid)
    timeHUXT = H5D_READ(dataset_time)

    H5D_CLOSE, dataset_lon
    H5D_CLOSE, dataset_dlon
    H5D_CLOSE, dataset_r
    H5D_CLOSE, dataset_cr_lon_init
    H5D_CLOSE, dataset_dt_out
    H5D_CLOSE, dataset_v_grid
    H5D_CLOSE, dataset_time

    H5F_CLOSE, file_id

    ; TODO: remove initial time to be set fixed

    bgswData = {time:timeHUXT,  $
                lon:lon,    $
                dlon:dlon,  $
                r:r,        $
                cr_lon_init:cr_lon_init,    $
                dt:dt,      $
                v_grid:v_grid,  $
                tinit:tinitial}
                
    return, bgswData
end