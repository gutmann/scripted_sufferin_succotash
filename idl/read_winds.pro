pro read_winds
	files=file_search('*.nc')
	x=117
	y=103
	z=10
	print, files[0]
	for i=0l,n_elements(files)-1 do begin
		ncid=ncdf_open(files[i])
		ncdf_varget, ncid,0,u
		ncdf_varget, ncid,1,v
		print,files[i], u[x,y,z],v[x,y,z]
		ncdf_close, ncid
		
	endfor
	
end
