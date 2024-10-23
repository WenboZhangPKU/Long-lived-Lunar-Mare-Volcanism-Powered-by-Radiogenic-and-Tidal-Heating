#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 14:21:05 2020

@author: p300
"""

import os, sys
import tool_readfiles as RF
import numpy as NP

'''
This is for CitcomS-3.3.1_new4/5
recalculating after radial nodes'redistribution
'''
#generating radial array of length of (2*noz-1)#
def Generation_radius(nr, rr):
    r =NP.zeros(nr[-1])
    j =1
    k =1
    for i in nr[:-1]:
        r[i-1]=rr[int(k)-1]
        s =(rr[k]-rr[k-1])/(nr[k]-nr[k-1])
        while j!= nr[int(k)]:
            r[j] =r[i-1] +(j-i+1)*s
            j +=1
        k +=1
    s =(rr[-1]-rr[-2])/(nr[-1]-nr[-2])
    while j!=nr[-1]:
        r[j] =r[-2] +(j-i+1)*s
        j +=1
    r[-1] =rr[-1]
    return r



def LinearInterpolate(src_r, src_temp, velo_r):
    velo_temp =[]
    for i in velo_r:
        for j in range(len(src_r)-1):
            if i>=src_r[j] and i<=src_r[j+1]:
                if i==src_r[j]:
                    velo_temp.append(src_temp[j])
                elif i==src_r[j+1]:
                    velo_temp.append(src_temp[j+1])
                else:
                    tmp =(src_temp[j+1]*(i-src_r[j]) +src_temp[j]*(src_r[j+1]-i))/(src_r[j+1]-src_r[j])
                    velo_temp.append(tmp)
                break
            else:
                continue
    return velo_temp

def BinaryInterpolate(src_r, src_field, dest_r):
    for i in dest_r:
        j = NP.argmin(NP.array([abs(k-i) for k in src_r]))
        yield src_field[j]

def BinaryInterpolate_Index(src_r, dest_r):
    for i in dest_r:
        j = NP.argmin(NP.array([abs(k-i) for k in src_r]))
        yield j

def Write_velo(cycle,nno,time,nproc,noz,vtk_down,vtk_up,velo_down,velo_up,src_r,velo_r):
    #read in temperature from upper and lower blocks#
    data_down, tmp =RF.read_vtk_file(vtk_down)
    data_up, tmp =RF.read_vtk_file(vtk_up)
    temp_down_total =data_down['temperature']
    temp_up_total =data_up['temperature']
    #head#
    with open(velo_down, 'a') as down:
        down.write("%d %d %.5e\n" %(cycle, nno, time))
        down.write("%3d %7d\n" %(nproc, nno))
    with open(velo_up, 'a') as up:
        up.write("%d %d %.5e\n" %(cycle, nno, time))
        up.write("%3d %7d\n" %(nproc, nno))
    #interpolate#
    for i in range(33*33):
        temp_down =[]
        temp_up =[]
        for j in range(33):
            temp_down.append(temp_down_total[i*33+j])
            temp_up.append(temp_up_total[i*33+j])
        #merge and delete duplicate nodes#
        temp_down.pop(-1)
        temp_down =NP.array(temp_down)
        temp_up =NP.array(temp_up)
        src_temp =NP.vstack((temp_down,temp_up))
        velo_temp=LinearInterpolate(src_r, src_temp, velo_r)
        #output data to velo files#
        velo_temp_down =velo_temp[:noz]
        velo_temp_up =velo_temp[(noz-1):]
        with open(velo_down, 'a') as down:
            for i in velo_temp_down:
                down.write("%.6e %.6e %.6e %.6e\n" %(0,0,0,i))
        with open(velo_up, 'a') as up:
            for i in velo_temp_up:
                up.write("%.6e %.6e %.6e %.6e\n" %(0,0,0,i))


def main_writing(solution_cycle_init):
    #mesh information#
    nproc =96
    noz =33
    nno =33*33*33,
    elapsed_time =1.1974e-04 #useless, indeed#
    #generating source and target radial arrays#
    nr=[1,3,45,65]
    rr=[1.9540e-01,2.2410e-01,8.8793e-01,1.0000e+00]
    src_r =Generation_radius(nr, rr)
    nr=[1,33,55,65]    #cannot be too arbitary otherwise the error will be Exit 10
    rr=[1.9540e-01,0.57502,0.95118,1.0000e+00]
    velo_r =Generation_radius(nr, rr)
    #find lower and upper source temperature files#
    #create corresponding upper and lower velo files#
    vtk_path ='/home/p300/newdisk/back/input/D150v1e-1V1e21H50_20201011'
    velo_path ='/home/p300/newdisk/back/input/D150v1e-1V1e21H50_202010151/temp0'
    if os.path.exists(velo_path)==0:
        os.mkdir(velo_path)
    for i in range(int(nproc/2)):
        vtk_down ='a.proc'+str(2*i)+'.'+str(solution_cycle_init)+'.vts'
        vtk_down =os.path.join(vtk_path,vtk_down)
        vtk_up ='a.proc'+str(2*i+1)+'.'+str(solution_cycle_init)+'.vts'
        vtk_up =os.path.join(vtk_path,vtk_up)
        velo_down ='a.velo'+'.'+str(2*i)+'.'+str(solution_cycle_init)
        velo_down =os.path.join(velo_path,velo_down)
        velo_up ='a.velo'+'.'+str(2*i+1)+'.'+str(solution_cycle_init)
        velo_up =os.path.join(velo_path,velo_up)
        #output velo files#
        Write_velo(solution_cycle_init,nno,elapsed_time,nproc,noz,vtk_down,vtk_up,velo_down,velo_up,src_r,velo_r)
