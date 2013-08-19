import numpy as np

def load_intercep():
    f=open("INTERCEP.TXT");i=0;current_length=1000000
    data=np.zeros((1000000,15))
    dates=[]
    for l in f:
        try:
            n=len(str.split(l,",")[1:])
            data[i,:n]=np.array(l.split(",")[1:]).astype("f")
            dates.append(l.split(",")[0])
            i+=1
        except:
            print(i,l)
            i+=1
        if i>=current_length:
            new_data=np.zeros((current_length*2,data.shape[1]))
            new_data[:data.shape[0],:]=data
            data=new_data
            current_length*=2
    f.close()
    return dates,data[:i,...]
