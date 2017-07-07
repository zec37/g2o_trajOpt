#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 11:10:09 2017

@author: yj
"""
import matplotlib.pyplot as plt

file_num=80
for i in range(file_num):
    file_name='./data/traj'+str(i+1)+'.txt'
    f=open(file_name)
    x=[]
    y=[]
    for line in f:
        data=line.split()
        x.append(float(data[1]))
        y.append(float(data[2]))
	#print(data[0])
    plt.plot(x,y)
    plt.show()
