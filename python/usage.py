#!/usr/bin/env python
import time
import subprocess

INTERVAL = 0.3

def getTimeList():
    """
    Fetches a list of time units the cpu has spent in various modes
    Detailed explanation at http://www.linuxhowtos.org/System/procstat.htm
    """
    cpuStats = file("/proc/stat", "r").readline()
    columns = cpuStats.replace("cpu", "").split(" ")
    return map(int, filter(None, columns))

def deltaTime(interval):
    """
    Returns the difference of the cpu statistics returned by getTimeList
    that occurred in the given time delta
    """
    timeList1 = getTimeList()
    time.sleep(interval)
    timeList2 = getTimeList()
    return [(t2-t1) for t1, t2 in zip(timeList1, timeList2)]

def get_usage():
    """
    Returns the cpu load as a value from the interval [0.0, 1.0]
    """
    dt = list(deltaTime(INTERVAL))
    idle_time = float(dt[3])
    total_time = sum(dt)
    load = 1-(idle_time/total_time)
    return load

def get_load():
    output = subprocess.Popen(['uptime'], stdout=subprocess.PIPE).communicate()[0]
    load=" ".join(output.split()[-3:])
    return load

#while True:
print "CPU usage=%.2f%%   %s" % (get_usage()*100.0,get_load())
#    time.sleep(0.1)
