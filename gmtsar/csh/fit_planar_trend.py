#!/usr/bin/python

import subprocess
import os.path
import numpy as np
import sys

def get_trend(x,y,z):
    G=np.column_stack((np.ones(np.shape(x)),x,y))
    #p=(G'*G)\G'*z
    p=np.dot(np.dot(np.linalg.inv(np.dot(np.transpose(G),G)),np.transpose(G)),z)
    return p

if __name__ == '__main__':

    # load data from the unwrapped subset
    unwrap=sys.argv[1] # eg. unwrap.dat
    data=np.loadtxt(unwrap)
    x=data[:,0]
    y=data[:,1]
    z=data[:,2]

    # find the trend
    trend_params=get_trend(x,y,z)
    print('%.15f %.15f %.15f') %(trend_params[0], trend_params[1],trend_params[2])


    
## END ##
