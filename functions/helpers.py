import time

def nice_time():
    return time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime())