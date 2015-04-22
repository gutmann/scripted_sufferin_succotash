; reads in CDP data, and plots CDP statistics on a global hexagonal grid.
; Ming Ge March 2015
; edited to run on hydro-c1 : Ethan Gutmann April 2015

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

function load_track_data(filename, x, y, grid_size)
begin
	print("Loading: " + filename + " " + grid_size)
    track_count            = ispan(1, grid_size, 1)
    track_count@units      = ""
    track_count@_FillValue = -99
    track_count            = 0

    track_cdp            = new((/dimsizes(track_count), 8000/), "float")
    track_cdp@units      = ""
    track_cdp@_FillValue = -999.

    track_cdp_med = new((/dimsizes(track_count)/), "float")
    track_cdp_med = 0.
    track_cdp_max = track_cdp_med
    track_cdp_mean= track_cdp_med

  ; read CDP 
    flnm_cdp = filename ;dir_i + "cdp_" + basin + ".txt"
    n_row = numAsciiRow(flnm_cdp)
    n_col = numAsciiCol(flnm_cdp)
    print(n_row + " " + n_col)
    ; n_col = 7
    data  = asciiread(flnm_cdp, (/n_row, n_col/), "float")
    lat   = data(:,1)
    lon   = data(:,2)
    cdp   = data(:,3) 

    do tc_point = 0, n_row -1
        gcdist  = gc_latlon(lat(tc_point), lon(tc_point), y, x, 30, 4) ; return distance in km
        ind_min = ind(gcdist.eq.min(gcdist)) ; where is the closest hexagon center

        track_cdp(ind_min, track_count(ind_min) ) = (/cdp(tc_point)/)
        track_count(ind_min) = track_count(ind_min) + 1
        delete([/gcdist, ind_min/])
    end do	; end of tc_point-loop

    do ttt = 0, grid_size - 1
      if(.not.all(ismissing(track_cdp(ttt,:)))) then
        track_cdp_max (ttt) =        max(track_cdp (ttt,:))
        track_cdp_mean(ttt) =    dim_avg(track_cdp (ttt,:))
        track_cdp_med (ttt) = dim_median(track_cdp (ttt,:))
      end if
    end do

    track_cdp_max = where(ismissing(track_cdp_max), 0, track_cdp_max)
    track_cdp_mean= where(ismissing(track_cdp_mean),0, track_cdp_mean)
    track_cdp_med = where(ismissing(track_cdp_med), 0, track_cdp_med)
	
	return(track_cdp_mean)
	
end

begin
  dir_i = "./"
  dir_o = "./"

  basin = "ebr64_atlc"

  r2d = 57.2957795           

  f   = addfile(dir_i + "C02562.global.nc","r")

  x  = f->grid_center_lon * r2d
  y  = f->grid_center_lat * r2d
  cx = f->grid_corner_lon * r2d
  cy = f->grid_corner_lat * r2d

  x@units = "degrees_east"
  y@units = "degrees_north"

  grid_size = dimsizes(x)

  track_cdp_mean = load_track_data(dir_i + "cdp_" + basin + ".txt",     x, y, grid_size)
  track_cdp_pgw  = load_track_data(dir_i + "cdp_" + basin + "_pgw.txt", x, y, grid_size)
  
  wks = gsn_open_wks("pdf", dir_o + "geodesic_cdp_pgw")
  ; gsn_define_colormap (wks, "wh-bl-gr-ye-re")
  gsn_define_colormap (wks, "WhiteYellowOrangeRed")
  
  map = new (2 , "graphic")          ; 1d array to hold plots

  res                      = True
  res@gsnDraw 	           = False
  res@gsnFrame	           = False
  res@gsnMaximize          = True
  res@gsnSpreadColors      = True

  res@sfXArray             = x
  res@sfYArray             = y
  res@sfXCellBounds        = cx
  res@sfYCellBounds        = cy

  res@lbLabelAutoStride    = True
  res@lbBoxLinesOn         = False

  lat_center_basin  =  30.
  lon_center_basin  = -60.
  res@mpMinLatF     = lat_center_basin - 35.
  res@mpMaxLatF     = lat_center_basin + 35.

  res@mpProjection      = "Orthographic"
  res@mpDataBaseVersion = "MediumRes"
  res@mpGridAndLimbOn   = True
  res@mpCenterLatF      = lat_center_basin
  res@mpCenterLonF      = lon_center_basin
  res@mpFillOn          = False
  res@mpLimitMode          = "LatLon"
  res@mpMinLonF            = lon_center_basin - 50.
  res@mpMaxLonF            = lon_center_basin + 50.
  
; CDP
  mnmxint                  = nice_mnmxintvl(0.000001, max(track_cdp_mean), 21, False)
  res@cnFillOn             = True
  res@cnFillMode           = "RasterFill"
  res@cnLinesOn            = False
  res@cnLineLabelsOn       = False
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF       = 0.4 ; mnmxint(0)
  res@cnMaxLevelValF       = max(track_cdp_mean) ;10. ; mnmxint(1)
  res@cnLevelSpacingF      = 0.4 ; mnmxint(2)

  res@tiMainString          = "Current CDP"
  map(0) = gsn_csm_contour_map(wks, track_cdp_mean, res)
  
  gsn_define_colormap (wks, "BlWhRe")
  res@cnMinLevelValF       = -5.1
  res@cnMaxLevelValF       = 5.1
  res@cnLevelSpacingF      = 0.2
  res@tiMainString          = "(Future - Current) CDP"
  map(1) = gsn_csm_contour_map(wks, track_cdp_pgw - track_cdp_mean, res)

  ; res@tiMainString          = "Max CDP"
  ; map(1) = gsn_csm_contour_map(wks, track_cdp_max, res)

  resP                = True
  resP@gsnMaximize    = True    
  gsn_panel(wks, map, (/1,2/), resP)             

end
