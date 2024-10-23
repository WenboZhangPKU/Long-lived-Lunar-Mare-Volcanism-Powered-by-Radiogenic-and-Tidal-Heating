#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Jan. 15th, 2024
"""

import os, sys
import tool_readfiles as RF
import numpy as NP

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

def Find_MaterialLayer(radius, rr):
    for MatLayer, R_interface in enumerate(rr):
        if R_interface < radius:
            continue
        else:
            break
    
    return MatLayer
