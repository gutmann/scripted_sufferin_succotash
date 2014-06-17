import glob,os
import copy

class RunCollector(object):
    """Simple object to collect data from the driver. 
    This is particularly useful if the file combinations must be collected 
    prior to processing, for example if correlations between variables are
    to be examined. 
    
    Note that the extra argument passed around gives you another option to deal
    with this problem by storing information about the previous call(s) in the extra argument. 
    
    Usage:
        collection=RunCollector()
        drive(collection.collect)
        for current_run in collection:
            do something (possibly getting more ellements from collection)
    """
    _allmodels=[]
    def collect(self,files,v,output_base,combinations,extra):
        self._allmodels.append([files,v,output_base,combinations,extra])
    def __iter__():
        return self
    def next(self):
        return self._allmodels.pop(0)
    __next__=next
    def getCollection(self):
        return copy.deepcopy(self._allmodels)
        
def drive(foo,yearsearch="200*",obs=True,stat=True,runforce=False,extra=[None],
          methods=["SDmon","SDe","CAe","CAe0","CAe1","SDe0","SDe1","CA","SD","SARe0","SARe1"],
          BCs=["BC",""],
          models=["ncep","narr"],
          resolutions=["6km","12km"],
          variables=["pr","tasmax","tasmin"],
          extendedmethods=False):
    """Base driver for many statistical downscaling comparison routines
    
    Usage:
        def processing_function(files,variable,output_base,combinations,extra):
            do something with files
            
        def main():
            drive(processing_function)
    
    Finds all valid method,forcing,resolution,variable combinations that have files
    and calls an external method with the files, variable name, outputbase, and any extra information
    
    foo: function to be called for each combination
    yearsearch: search string for the years to be collected
    obs: boolean If true, return observed datasets
    stat: boolean If true, return statistical datasets
    
    methods, BC status, models, resolutions, and variables can all be supplied as lists
    """
    
    if extendedmethods:
        methods.extend(["CAcold","CAdry","SDcold","SDdry","CAhot","CAwet","SDhot","SDwet"])
    if stat:
        for m in methods:
            for b in BCs:
                for forc in models:
                    for r in resolutions:
                        for v in variables:
                            # print(m,forc,r,v,b)
                            if (forc=="ncep"):
                                gridtype="gauss"
                            else:
                                gridtype=""
                            # note SAR and SDmon do not have non-bias corrected versions
                            if ((m[:3]=="SAR") or (m[:3]=="SDm")) and b=="":
                                pass
                            else:
                                if type(yearsearch)==list:
                                    files=[]
                                    for ys in yearsearch:
                                        filesearch=m+"/"+forc+"/"+v+"/"+b+m[:2]+"*"+r+"*"+gridtype+"*"+ys+"*.nc"
                                        files.extend(glob.glob(filesearch))
                                else:
                                    filesearch=m+"/"+forc+"/"+v+"/"+b+m[:2]+"*"+r+"*"+gridtype+"*"+yearsearch+"*.nc"
                                    files=glob.glob(filesearch)
                                    
                                output_base="-".join([m,forc,v,b+r])
                                print(len(files),output_base)
                                if len(files)>1:
                                    foo(files,v,output_base,[m,b,forc,r,v],extra)
                                else:
                                    print("--------Not enough files in "+ output_base.replace("-"," "))
                                    print("     "+filesearch)

    if obs:
        obs_base="../obs/"
        for r in resolutions:
            if r=="12km":
                o="maurer.125"
            elif r=="6km":
                o="uw.0625"
            if r=="4km":
                o="uw.4km"
                
            for v in variables:
                filesearch=obs_base+o+"/"+v+"/*"+v+"."+yearsearch+".nc"
                files=glob.glob(filesearch)
                output_base="-".join(["obs",o,v])
                print(len(files),output_base)
                if len(files)>1:
                    foo(files,v,output_base,[o,"BC","obs",r,v],extra)
                else:
                    print("--------Not enough files in "+ output_base.replace("-"," "))
     
    if runforce:
        force_base="../../forcing/"
        for forc in models:
                
            for v in variables:
                if forc=="ncep":
                    additional="gauss*"
                else:
                    additional="*"
                filesearch=force_base+forc+"/"+v+"/*"+v+"*"+additional+yearsearch+".daily.nc"
                files=glob.glob(filesearch)
                output_base="-".join(["forcing",forc,v])
                print(len(files),output_base)
                if len(files)>1:
                    foo(files,v,output_base,[forc,"BC","forcing","12km",v],extra)
                else:
                    print(filesearch)
                    print("--------Not enough files in "+ output_base.replace("-"," "))
        