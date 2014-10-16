from matplotlib.colors import hsv_to_rgb, rgb_to_hsv
import numpy as np

class LightSource(object):
    """
    Create a light source coming from the specified azimuth and elevation.
    Angles are in degrees, with the azimuth measured
    clockwise from north and elevation up from the zero plane of the surface.
    The :meth:`shade` is used to produce rgb values for a shaded relief image
    given a data array.
    
    Original in matplotlib.colors, modified 1/18/2012 - Ethan Gutmann to provide more accurate/realistic shading
    In particular: 
        added deltax,y keywords to allow grid to be on a dx/dy increment other than 1
        calculate cos(theta) term and just multiply HSV-Value by that term
        removed modifications to Saturation, and Hue as well as lots of normalization code
    """
    def __init__(self,azdeg=315,altdeg=45,\
                 hsv_min_val=0,hsv_max_val=1,hsv_min_sat=1,hsv_max_sat=0):
       """
       Specify the azimuth (measured clockwise from south) and altitude
       (measured up from the plane of the surface) of the light source
       in degrees.

       The color of the resulting image will be darkened
       by moving the (s,v) values (in hsv colorspace) toward
       (hsv_min_sat, hsv_min_val) in the shaded regions, or
       lightened by sliding (s,v) toward
       (hsv_max_sat hsv_max_val) in regions that are illuminated.
       The default extremes are chose so that completely shaded points
       are nearly black (s = 1, v = 0) and completely illuminated points
       are nearly white (s = 0, v = 1).
       """
       self.azdeg = azdeg
       self.altdeg = altdeg
       self.hsv_min_val = hsv_min_val
       self.hsv_max_val = hsv_max_val
       self.hsv_min_sat = hsv_min_sat
       self.hsv_max_sat = hsv_max_sat

    def shade(self,data,cmap,deltax=1.0,deltay=1.0,minval=None):
        """
        Take the input data array, convert to HSV values in the
        given colormap, then adjust those color values
        to given the impression of a shaded relief map with a
        specified light source.
        RGBA values are returned, which can then be used to
        plot the shaded image with imshow.
        
        Actually, conversion to and from HSV is done in shade_rgb
        """

        if minval==None:minval=data.min()
        normdata=(data-minval)/(data.max()-minval)
        # not necessary, cmap will cutoff values below 0?
        # normdata=np.choose(normdata<0,(normdata,0))
        rgb0 = cmap(normdata)
        rgb1 = self.shade_rgb(rgb0, elevation=data,deltax=deltax,deltay=deltay)
        rgb0[:,:,:3] = rgb1[:,:,:3]
        return rgb0

    def shade_rgb(self,rgb, elevation, deltax=1.0,deltay=1.0,fraction=1.0):
        """
        Take the input RGB array (ny*nx*3) adjust their color values
        to given the impression of a shaded relief map with a
        specified light source using the elevation (ny*nx).
        A new RGB array ((ny*nx*3)) is returned.
        """
        # imagine an artificial sun placed at infinity in
        # some azimuth and elevation position illuminating our surface. The parts of
        # the surface that slope toward the sun should brighten while those sides
        # facing away should become darker.
        # convert alt, az to radians
        az = (self.azdeg-90)*np.pi/180.0
        alt = self.altdeg*np.pi/180.0
        # gradient in x and y directions
        dzx, dzy = np.gradient(elevation)
        d0=np.zeros(elevation.shape)
        dx=np.zeros(elevation.shape)+deltax
        dy=np.zeros(elevation.shape)+deltay
        
        #  setup surface vectors in X and Y directions
        vecshape=(d0.shape[0],d0.shape[1],1)
        vec1=np.concatenate((dx.reshape(vecshape),d0.reshape(vecshape),dzx.reshape(vecshape)),axis=2)
        vec2=np.concatenate((d0.reshape(vecshape),dy.reshape(vecshape),dzy.reshape(vecshape)),axis=2)
        # compute normal vectors at each point on the surface (normal = cross product of two vectors on the surface)
        normals=np.cross(vec1,vec2)

        # compute a vector pointing towards the sun
        cosalt=np.cos(alt)
        sunvec=np.array([cosalt*np.sin(az),cosalt*np.cos(az),np.sin(alt)])
        #compute cos(theta)  (=dot product of normal and sun vectors divided by the magnitude of each)
        shade=np.dot(normals,sunvec)/(np.sqrt((normals**2).sum(axis=2))*np.sqrt((sunvec**2).sum()))
        # set a minimum shadeing value so we don't set anything completely black
        intensity=np.choose(shade<0.05,(shade,0.05))
        # normalize so that the entire map never ends up a blah shade of grey, might as well have the brightest part still be bright
        intensity/=np.max(intensity)

        hsv = rgb_to_hsv(rgb[:,:,0:3])
        # modify hsv values to simulate illumination.
        hsv[:,:,2]*=intensity

        # convert modified hsv back to rgb.
        outputrgb=hsv_to_rgb(hsv)
        if rgb.shape[2]>3:
            tmp=rgb.copy()
            tmp[:,:,0:3]=outputrgb
            return tmp
        else:
            return outputrgb
