import numpy as np
from glob import glob
import matplotlib.pyplot as plt

filesearch="20*/*/*.png"

def main():
    files=glob(filesearch)
    files.sort()
    
    albedos=np.zeros(len(files))
    colors=np.zeros(len(files))
    lengths=np.zeros(len(files))
    for i,f in enumerate(files):
        d=plt.imread(f)
        print(f, i, len(files))
        curdata=d[215:345,400:480,:]
        r=curdata[:,:,0]
        g=curdata[:,:,1]
        good=np.where(r>0.3)
        colors[i]=np.mean(g[good]/r[good])
        albedos[i]=np.mean(g[good]+r[good])/2
        lengths[i]=len(good[0])
            
    
    plt.plot(colors,'x')
    plt.savefig("colors.png")
    return (colors,albedos,lengths)

if __name__ == '__main__':
    main()