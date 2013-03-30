# 
# Downloaded from http://www.ariel.com.au/a/python-point-int-poly.html
# 
# points_in_poly implemented by edg Dec. 17 2012
# 

def point_in_poly(x,y,poly):
    """Determine if a point is inside a given polygon or not
         Polygon is a list of (x,y) pairs. This fuction returns True or False.  
         The algorithm is called the Ray Casting Method.
         
       from: http://www.ariel.com.au/a/python-point-int-poly.html
       
       Edge and Vertex testing from: 
           http://geospatialpython.com/2011/08/point-in-polygon-2-on-line.html
    """
    # check if point is a vertex
    if (x,y) in poly: return True

    # check if point is on a boundary
    n = len(poly)
    for i in range(n):
        p1 = poly[i-1]
        p2 = poly[i]
        if p1[1] == p2[1] and p1[1] == y and x > min(p1[0], p2[0]) and x < max(p1[0], p2[0]):
            return True
      
    inside = False

    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x,p1y = p2x,p2y

    return inside

def points_in_poly(x,y,poly):
    """Same algorithm as above but requires numpy arrays for x and y
        all points are processed simultaneously (fast)
        returns a boolean array of len(x.size) 
        NOTE: this does not currently include the edge cases at the beginning
        """
    # only require numpy as runtime so point_in_poly is available even without numpy installed
    import numpy as np
    # maxlength=0.3
    
    n = len(poly)
    crossings = np.zeros(x.size)

    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        # note if a side length is greater than a max, 
        # it is probably a connection to e.g. an island and should be ignored
        # if (np.sqrt((p2x-p1x)**2+(p2y-p1y)**2))<maxlength:
        # find all points that match the first three if statements in point_in_poly
        # curpts=np.where((y>min(p1y,p2y)) & (y<=max(p1y,p2y)) & (x<=max(p1x,p2x)))
        ytest1=np.where(y>min(p1y,p2y))
        curpts=[[]]
        if len(ytest1[0])>0:
            ytest2=np.where(y[ytest1]<=max(p1y,p2y))
            if len(ytest2[0])>0:
                curpts=np.where(x[ytest1[0][ytest2[0]]]<=max(p1x,p2x))
        # if there are matching points...
        if len(curpts[0])>0:
            curpts=ytest1[0][ytest2[0]][curpts[0]]
            # inner code from point_in_poly
            if p1y != p2y:
                xinters = (y[curpts]-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
            # find all matching points (from among the previous matching points)
            test=np.where((p1x==p2x)|(x[curpts]<=xinters))
            # if there are matching points...
            if len(test[0])>0:
                # these points just crossed an edge so add on to their edge crossing count
                try:
                    crossings[curpts[test]]+=1
                except Exception as e:
                    print(e)
                    print(curpts)
                    print(test)
        p1x,p1y = p2x,p2y
    
    # a point is inside the polygon there were an odd number of edge crossings. 
    inside= (crossings%2)==1

    # check if point is a vertex
    for i in range(n):
        p1x,p1y = poly[i]
        onxvertex=np.where(x==p1x)
        if len(onxvertex[0])>0:
            onxyvertex=np.where(y[onxvertex]==p1y)
            if len(onxyvertex[0])>0:
                inside[tmp]=1
    
    # check if point is on a boundary
    for i in range(n):
        p1 = None
        p2 = None
        if i==0:
            p1 = poly[-1]
            p2 = poly[0]
        else:
            p1 = poly[i-1]
            p2 = poly[i]
        if p1[1] == p2[1]:
            # if np.abs(p1[0]-p2[0])<maxlength:
            onborder=np.where((p1[1] == y) & (x > min(p1[0], p2[0])) & (x < max(p1[0], p2[0])))
            if len(onborder[0])>0:
                inside[onborder]=1
                
    return inside


# ## Test
# 
# polygon = [(0,10),(10,10),(10,0),(0,0)]
# 
# point_x = 5
# point_y = 5
# 
# ## Call the fuction with the points and the polygon
# print point_in_poly(point_x,point_y,polygon)
# 
# print points_in_poly(np.arange(10)*2-3,np.arange(10)*2-3,polygon)


# see http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html for an explanation of the c code below:
# int pnpoly(int nvert, float *vertx, float *verty, float testx, float testy)
# {
#   int i, j, c = 0;
#   for (i = 0, j = nvert-1; i < nvert; j = i++) {
#     if ( ((verty[i]>testy) != (verty[j]>testy)) &&
#      (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
#        c = !c;
#   }
#   return c;
# }
# LICENSE for C-code commented above. 
# Copyright (c) 1970-2003, Wm. Randolph Franklin
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of
#  this software and associated documentation files (the "Software"), to deal in 
#  the Software without restriction, including without limitation the rights to 
#  use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
#  of the Software, and to permit persons to whom the Software is furnished to do
#  so, subject to the following conditions:
# 
# 1) Redistributions of source code must retain the above copyright notice, this list
#  of conditions and the following disclaimers.
# 2) Redistributions in binary form must reproduce the above copyright notice in the 
# documentation and/or other materials provided with the distribution.
# 3) The name of W. Randolph Franklin may not be used to endorse or promote products 
# derived from this Software without specific prior written permission.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS 
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER 
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
