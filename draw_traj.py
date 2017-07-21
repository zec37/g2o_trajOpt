#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 11:10:09 2017

@author: yj
"""
import matplotlib.pyplot as plt
import math

while(1):
    num_recv = raw_input("Which traj? ")
    num = num_recv.split(" ")
    if(int(num[1]) == 1):
        file_name='./data_ori/traj'+str(int(num[0])+1)+'.txt'
        f=open(file_name)
        x=[]
        y=[]
        for line in f:
            data=line.split()
            x.append(float(data[1]))
            y.append(float(data[2]))
        plt.plot(y,x)
        plt.show()
    else:
        n = int(num[1])
        size_x = round(math.sqrt(n))
        size_y = math.ceil(float(n) / float(size_x))
        plt.figure(1)
        for i in range(n):
            plt.figure(1)
            plt.subplot(size_x, size_y, i + 1)
            file_name = './data_ori/traj' + str(int(num[0]) + i + 1) + '.txt'
            f = open(file_name)
            x = []
            y = []
            for line in f:
                data = line.split()
                x.append( float(data[1]))
                y.append( float(data[2]))
            plt.plot(y,x)
        plt.show()
