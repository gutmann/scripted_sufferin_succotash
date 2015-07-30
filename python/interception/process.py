import datetime
import numpy as np
from numpy.fft import fft
import matplotlib.pyplot as plt

# ## simple CLI code to load, process sway and environmental paramters and make some plots
# import mygis,glob,date_fun, process
# 
### LOAD INTERCEPTOMETER DATA (and calc sway)
# d=mygis.read_nc("newinterception.nc").data
# mjd=mygis.read_nc("newinterception.nc","mjd").data
# tmpmjd=mjd[8000::10000]
# ms=(d[:,0].copy()+32768).astype(np.ulonglong)
# start=0
# end=0
# for i in range(1,len(ms)):
#     if ms[i]<ms[i-1]:
#         print(i)
#         if start!=0:
#             end=i
#             ms[start:end]+=ms[start-1]
#             start=end
#         else:
#             start=i
# f3=process.calc_scargle_frequencies(d[:,7],time=ms)
# newoutput=np.zeros(3099)
# for i in range(len(newoutput)):
#     curdata=f3[i*100:(i+1)*100,0][f3[i*100:(i+1)*100,1]>15]
#     curdata=curdata[(1/curdata)<2500]
#     if len(curdata)>1:
#         newoutput[i]=np.mean(curdata)
### LOAD ENVIRONMENTAL DATA (assumes location of files)
# snotel=load_data.cols("snotel_precip.txt")
# swe=load_data.cols("snotel_swe.txt")
# snodates=date_fun.datearr2datetime(snotel[:,1:6])
# temperature=np.concatenate(mygis.read_files("nwt_data/*.nc","T12"))
# speed=np.concatenate(mygis.read_files("nwt_data/nwt.131*","w_spd_25m"))
# wdir=np.concatenate(mygis.read_files("nwt_data/nwt.131*","w_dir_25m"))
# wdates=[]
# files=glob.glob("nwt_data/*")
# for f in files:
#      wdates.append(datetime.datetime(2000+int(f.split(".")[1][:2]),
#                  int(f.split(".")[1][2:4]),int(f.split(".")[1][4:])))
# for i in range(len(files)):
#     for j in range(1,24*12):
#         wdates.insert(i*24*12+j,wdates[i*24*12]+datetime.timedelta(0,5*60*j))
# 
# ## PLOTS
# clf();plot(tmpdates,1/newoutput,'x',label="Sway Period [ms]")
# plot(snodates,snotel[:,-1]*200+1200,label="SNOTEL precip [in*200]")
# plot(snodates,snotel[:,-1]*200+1200,label="SNOTEL precip [in*200]")
# plot(wdates,wdir+2100,'x',label="Wind Dir [mv+2100]")
# plot(wdates,speed*10+2000,'x',label="Wind Speed [m/s*10+2000]")
# plot(wdates,temperature*10+1000,'x',label="Temperature [C*10+1000]")
# legend(ncol=3)
# ylim(900,2800)
# plot([snodates[0],snodates[-1]],[1000,1000],color="yellow",linewidth=2)
# plot([snodates[0],snodates[-1]],[2000,2000],color="magenta",linewidth=2)


def std_over_dates(data1,date1,date2):
    """Useful for calculating wind speed [speed = f(sqrt(stddev)) ]"""
    start=0
    current=0
    outputdata=np.zeros(len(date2))
    
    for i in range(len(data1)):
        if date1[i]>date2[current]:
            if start==0:
                start=i
            else:
                end=i
                outputdata[current]=np.std(data1[start:end])
                start=end
                current+=1
    
    return outputdata
    


def dtfromstring(datestr):
    """Convert a date string (e.g. 2010/10/2 13:3:30) to a datetime object"""
    year=datestr[:4]
    month=datestr.split("/")[1]
    day=datestr.split("/")[-1].split()[0]
    hour=datestr.split(":")[0].split()[-1]
    minute=datestr.split(":")[1]
    second=datestr.split(":")[-1].strip()
    return datetime.datetime(int(year),int(month),int(day),int(hour),int(minute),int(second))

def load_nc_data(filename):
    """docstring for load_nc_data"""
    import mygis,date_fun
    data=mygis.read_nc(filename).data
    dates=date_fun.mjd2datetime(mygis.read_nc(filename,"mjd").data)
    return dates,data

def load_data(filename,ncols=None):
    """Read an interceptometer data file
    
    Assumes a comma delimited text file with no header lines
    Reads each line separating the first column out as a date and converting using dtfromstrin
    The rest of the line is converted to 32bit integers (milliseconds need this)
    
    Returns a tuple of (date, data)
    
    Each line is converted within it's own try/except block, and lines that fail are printed to the screen
    """
    f=open(filename)
    d=[]
    dates=[]
    d=None
    ntimes=25920000 # 30days of 10hz data
    additional_time=ntimes
    
    i=long(0)
    for l in f:
        try:
            curdate=dtfromstring(l.split(",")[0])
            curdata=np.array(l.split(",")[1:],dtype=np.int32)
            if ncols==None:
                ncols=len(curdata)
            if len(curdata)==ncols:
                dates.append(curdate)
                if d==None:
                    d=np.zeros((ntimes,ncols),dtype=np.int32)
                
                if (i>=ntimes):
                    print(str(ntimes)+" "+str(curdate))
                    ntimes+=additional_time # add ~one month of memory ~1.5GB at a time (avoids doubling near the last call and requireing 128GB of RAM)
                    
                    newd=np.zeros((ntimes,ncols),dtype=np.int32)
                    newd[:i,:]=d
                    del d # help the gc out, this could be a LOT of data ~20GB
                    d=newd
                    del newd
                if ((i%864000)==0):
                    print(dates[i])
                d[i,:]=curdata
                i+=1
        except:
            print(l[:min(len("YYYY/MM/DD 00:00"),len(l))])
    f.close()
    
    return (dates,d[:i,:])

def calc_scargle_frequencies(data,time,verbose=True):
    import scargle
    
    window=1024.0
    window=128.0
    stepsize=100
    fy_range=[10,300]
    
    output=np.zeros((len(data)/stepsize,2))
    count=0
    for i in np.arange(0,len(data)-window,stepsize):
        if verbose:
            if (count%500)==0:
                print("cur={} of {} = {}%".format(i,len(data),int(100*i/len(data))) )
        count+=1
        
        fx,fy, nout, jmax, prob = scargle.fasper(time[i:i+window],data[i:i+window], 6., 6.)
        output[i/stepsize,0]=fx[np.argmax(fy[fy_range[0]:fy_range[1]])+fy_range[0]]
        output[i/stepsize,1]=np.max(fy[fy_range[0]:fy_range[1]])

    return output[:i,:]



def calc_frequencies(data,time=None):
    """ Calculate frequency response for a time series of data

    Uses a moving 1024 sample window to calculate peak FFT response within a window. 
    Moves 100 samples at a time so the output is has n/100 samples
    Converted to period with 1024/fft/sample_hz
    If time is given (as array of elapsed milli-seconds), uses that to calculate sample_hz, else assumes 10hz
    Ideally this should use a Lomb-Scargle type transform to incorporate the time stamp into the FFT itself. 
    """
    window=1024.0
    stepsize=100
    fft_range=[40,500]
    
    output=np.zeros((len(data)/stepsize,2))
    for i in np.arange(0,len(data)-window,stepsize):
        ft=np.abs(fft(data[i:i+window]))
        output[i/stepsize,0]=np.argmax(ft[fft_range[0]:fft_range[1]])+fft_range[0]
        output[i/stepsize,1]=np.max(ft[fft_range[0]:fft_range[1]])
        
        if time !=None:
            curdt=np.median(time[i+1:i+window]-time[i:i+window-1])
            sample_hz=1000.0/curdt #convert time step in milliseconds into samples/second [hz]
        else:
            sample_hz=10.0
        output[i/stepsize,0]=window/output[i/stepsize,0]/sample_hz
        
    return output[:i,:]
    
def smooth_frequencies(freq,threshold=7,fthresh=0.0005):
    fs=np.zeros(freq.shape[0]/100)
    for j in range(fs.size):
        cur=freq[j*100:j*100+100,0][freq[j*100:j*100+100,1]>threshold]
        if len(cur)>1:
            cur=cur[cur>fthresh]
        if len(cur)>1:
            fs[j]=np.mean(cur)
    return fs

def sway_frequencies(data):
    """Calculate sway frequencies for a given time series of data
    
    Data should be as formatted by load_data (i.e. a tuple with the second element being the data)
    """
    output=[]
    # for i in range(data[1].shape[1]):
    for i in range(3,4):
        print(i)
        freq=calc_scargle_frequencies(data[1][:,i],time=data[1][:,0])
        # freq[:,0]=1024/freq[:,0]/10.5
        fs=np.zeros(freq.shape[0]/100)
        for j in range(fs.size):
            cur=freq[j*100:j*100+100,0][freq[j*100:j*100+100,1]>500]
            if len(cur)>1:
                fs[j]=np.mean(cur)
        sdates=data[0][5000::10000]
        sdates=sdates[:fs.size]
        output.append([freq,fs,sdates])
    return output

def write_sway_to_textfile(swaydata,textfile):
    sites=[1,2,3,4,5,6,7,9,10,11,12,13]
    with open(textfile,'w') as f:
        f.write("Date,Tree1W,Tree1E,Tree1S,Tree1N,Tree2E,Tree2W,Tree2S,Tree2N,Tree3E,Tree3W,Tree3S,Tree3N\n")
        for line in range(len(swaydata[0][1])):
            f.write(str(swaydata[0][2][line]))
            for site in range(len(swaydata)):
                f.write(","+str(swaydata[site][1][line]))
            f.write("\n")
    

def plot_sway_frequency(dates,data):
    # simple plot
    plt.plot(dates,data,'x')
    plt.ylim(1,2)
    plt.ylabel("Period [s]")
    plt.title("Sway Frequency vs Time")

        

def convert_to_netcdf(filename="INTERCEP.TXT",outputfile="interception_data",data=None):
    """Convert an interceptometer text datafile (or data tuple) into netcdf
    
    This is separated out a bit from the other routines because it requires a netcdf library be installed
    Along with my helper module (mygis) to write the netcdf file and I'm guessing most won't use this
    Also uses my date_fun module to convert the datetime objects into modified julian days (easier for netcdf)
    """
    import mygis
    from bunch import Bunch
    import date_fun
    
    if data==None:
        data=load_data(filename)
    
    mjd=date_fun.datetime2mjd(data[0])
    extravars=[Bunch(data=mjd,name="mjd",dims=("y",),dtype="d",
                attributes=Bunch(units="days",description="Modified Julian Day"))]
    
    desc="raw voltage read from interceptometer poteniometers and datalogger voltage levels"
    cols="0=milliseconds, 1-4=Tree1:W,E,S,N, 5-7,9=Tree2,8=Tree temperature, 10-13=Tree3, 14,15,16=junk,17=logger temp,18=G,19=Ex(~4.95v)"
    note="""
    To convert data to voltage, divide by 32768 and multiply by 6.144
    Also, normalize by the excitation voltage by dividing by the data in column 19 (0 based numbering)
    I'll figure out conversions for Temperature sensors later"""
    attributes=Bunch(description=desc,columns=cols,note=note)
    
    mygis.write(outputfile,
                    data[1],dtype="h",units="volts*32768/6.144",
                    extravars=extravars,attributes=attributes)