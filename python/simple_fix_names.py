#!/usr/bin/env python
import numpy as np
import glob
import os

def fix_names():
    oldfiles=glob.glob("BCSAR*")
    # oldfiles=glob.glob("added*")
    # if len(oldfiles)==len(newfiles):
    for o in oldfiles:
        fileparts=o.split("_")
        fileparts.insert(2,"4km_gauss")
        newfile="_".join(fileparts)
        # print(o,newfile)
        os.rename(o,newfile)
    # else:
    #     print(str(os.getcwd().split("/")[-2:])+" : \n       "+str(len(oldfiles))+"  "+str(len(newfiles)))

def main():
    topdirs=["SARe0","SARe1"]
    moddirs=["ncep"]# ,"narr"]
    vardirs=["pr","tasmax","tasmin"]
    for t in topdirs:
        for m in moddirs:
            os.chdir(t+"/"+m)
            for v in vardirs:
                os.chdir(v)
                print(os.getcwd())
                fix_names()
                os.chdir("../")
            os.chdir("../../")

if __name__ == '__main__':
    main()