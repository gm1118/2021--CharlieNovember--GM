#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 13 14:00:31 2021

@author: georgemercieca
"""

import numpy as np
import random as ran
import matplotlib.pyplot as plt
import logbin_2020 as L
import time


class OsloModel:
    """
    A Class for describing an iteration of the Oslo model with its own
    variables to describe the system's properties such as L, h and t.
    Has functions for driving the model and measuring its properties.
    """
    
    def __init__(self, L=10, p = 0.5):
        self._L = L           #Set system length to L
        self._p = p
        if self._L <= 0:
            raise Exception("Length is negative or 0")
        self._h = np.zeros(self._L + 1)    #Array containing site heights
        self._z = np.zeros(self._L)    #Array for site threshold slopes
        for i in range(L):
            r = ran.random()
            if r <= self._p:
                self._z[i] = 1
            else:
                self._z[i] = 2
        self._T = 0

        
        
    def drive(self, n):
        """
        Function for driving the Oslo Model 'n' times returns the 
        avalanche size 's' for the avalanches that have occured as an array
        """
        if n <= 0:
            raise Exception("Number of iterations is negative or 0")
        if isinstance(n,int) != True:
            raise Exception("'n' is not an integer")
        s = np.zeros(n)      # array for no. of avalanches per iteration
        while n > 0:
            self._h[0] += 1       #Drive the system
            i = 0                 #Variable for checking for relaxations
            if self._h[i] - self._h[i+1] <= self._z[i]:
                self._T += 1
                n -= 1            #If there is no relaxation at site 1 following drive, there will be none in this iteration
            else:
                s[len(s)-n] += 1             #add Avalanche
                self._h[i] -= 1   
                self._h[i+1] += 1
                r = ran.random()
                if r <= self._p:     #Create new threshold slope
                    self._z[i] = 1
                else:
                    self._z[i] = 2
                while i < self._L:  #Check entire system for avalanches
                    if self._h[i] - self._h[i+1] <= self._z[i]:
                        i += 1       #If no avalanche possible proceed to check next site
                    else:
                        s[len(s)-n] += 1
                        self._h[i] -= 1
                        self._h[i+1] += 1
                        r = ran.random()
                        if r <= self._p:
                            self._z[i] = 1
                        else:
                            self._z[i] = 2
                        if i == self._L-1:
                            self._h[i+1] = 0   #For final site 'h' will always be 0
                            i -= 1
                        elif i == 0:           
                            pass               #For initial site, there is no previous site to check, so check current site again before proceeding
                        else:
                            i -= 1       #First check to see if previous site has been affected by avalanche and needs to relax, then proceed forward
                self._T += 1
                n -= 1  #Once all possible avalanches have occured proceed to next iteration

        return s         #Return total number of avalanches
    
    def h(self, num):
        return self._h[num]
    
    def z(self, num):
        return (self._h[num]-self._h[num+1])
    
    def T(self):
        return self._T
                                       
    
    def Reinit(self):
        """
        Function for reinitialising an Oslo Model 
        """
        self._h = np.zeros(self._L + 1)   
        self._z = np.zeros(self._L)
        for i in range(self._L):
            r = ran.random()
            if r <= self._p:
                self._z[i] = 1
            else:
                self._z[i] = 2
        self._T = 0  
        
    def hAve(self, n):
        """
        Function for testing the OsloModel Class along with the drive function
        to see if it is working as intended using the test suggested in the
        project notes, n is the no. of values used for average height value.
        Returns the average height value at the 1st site
        """
        while self._T == sum(self._h):
            self.drive(1)
        h1 = np.zeros(n)
        for i in range (n):          #Drive system 'n' times and record height at 1st site each time
            self.drive(1)
            h1[i] = self.h(0)        #Add to array
        return (sum(h1)/n)           #Sum array and return average 
        
        
        
    def Test2(self):
        """
        Second Function for Testing the OsloModel Class and drive function
        by seeing if it can recreate the 1D BTW model by setting p = 1
        """
        self.Reinit()
        self._p = 1                       #Recreate 1D BTW model
        n = int(0.5*(self._L**2+self._L))      #Drive until expected recurrent configs
        self.drive(n)
        print (self._h)
        self.drive(1)
        if self._T == sum(self._h)+1:
            print("Only 1 grain has left the model as expected")
        else:
            raise Exception("Test Unsuccessful")
        self.drive(99)
        if self._T == sum(self._h)+100:
            print("Test also successful following 100 drives into Recurrent Configuration")
        else:
            raise Exception("Test not successful when attempting 100 drives - check model variable values")
        

    def htdata(self,T):
        """
        Function for obatining the data to plot h vs. t for a model of size L
        Returns an array for h, and the point where the model
        moves from transient to recurrent configurations. End time T can be set.
        """
        h = np.zeros(T+1)
        Tc = 0
        for i in range(T):
            h[i] = self.h(0)
            self.drive(1)
            if Tc == 0:
                if self._T != sum(self._h):
                    Tc = i
        h[T] = self.h(0)
        return h,Tc
    
    def TcAve(self,n=100):
        """
        Function for obtaining an average of values for Tc for a given system's parameters determined
        from 'n' values
        """
        Tcs = np.zeros(n)
        for i in range(n):
            self.Reinit()
            while self._T == sum(self._h):
                self.drive(1)
            Tcs[i] = self._T
        return sum(Tcs)/n
    
    def htdataAve(self, T, n):
        """
        Function for obtaining an average of h vs t data from multiple realisations of
        the same system
        """
        H = np.zeros(T+1)
        for i in range(n):
            self.Reinit()
            h, Tc = self.htdata(T)
            H += h
        H /= n
        return H
    

    def hVarAve(self, n):
        """
        Function for obtaining the variance for the average height
        """
        while self._T == sum(self._h):
            self.drive(1)
        h1Arr = np.zeros(n)
        h2Arr = np.zeros(n)
        for i in range (n):          #Drive system 'n' times and record height at 1st site each time
            self.drive(1)
            h2Arr[i] = (self.h(0))**2        #Add to array
            h1Arr[i] = self.h(0)
        h2 = sum(h2Arr)/n
        h1 = sum(h1Arr)/n
        hVar = h2 - h1**2
        return hVar
        
    
    def hProbdata(self, n):
        """
        Function for creating a dataset made up of two arrays, one for the various
        heights for the model, and the other for the frequency of said height, hence
        producing a model for the probability of each height
        """
        while self._T == sum(self._h):
            self.drive(1)
        hArr = np.zeros(0)
        hFreqArr = np.zeros(0)
        for i in range (n):          #Drive system 'n' times and record height at 1st site each time
            repeat = False    
            self.drive(1)
            h = self.h(0)
            for j in range(len(hArr)):
                if h == hArr[j]:
                    hFreqArr[j] += 1/n
                    repeat = True
                else:
                    pass
            if repeat != True:
               hArr = np.append(hArr, [h])
               hFreqArr = np.append(hFreqArr, [1/n])
        return hArr, hFreqArr
    
    def zdata(self, n):
        """
        Function for obtaining probability data for each possible value of z from
        the z values in this Oslo Model
        """
        while self._T == sum(self._h):
            self.drive(1)
        zArr = np.array([0,1,2])
        zFreqArr = np.zeros(3)
        for i in range(n):
            self.drive(1)
            for j in range(self._L):
                z = self.z(j)
                if z == 0:
                    zFreqArr[0] += 1/(n*self._L)
                elif z == 1:
                    zFreqArr[1] += 1/(n*self._L)
                else:
                    zFreqArr[2] += 1/(n*self._L)
        return zArr, zFreqArr
    
    def sdata(self, N):
        """
        Function for obtaining a data set of avalanche sizes for the Oslo Model
        Returns an array conatining all of the avalanches unordered
        """
        while self._T == sum(self._h):
            self.drive(1)
        s = self.drive(1)    
        sArr = np.array(int(s))
        for i in range(N-1):
            s = self.drive(1)
            sArr = np.append(sArr, int(s))
        return sArr
        

def InterpolatePoly(x, yin, xin):
    """
    Function for a lagrange-polynomial interpolation on xy data performed for
    a specific x-point
    """
    np.seterr(divide='ignore', invalid='ignore')
    if len(xin) != len(yin):
        raise Exception("Size of data sets must be equal")
    l = len(xin)

    nom = 1
    Intersum = 0
    Inter = np.zeros(l)      #create array for terms in the sum
    
    for i in range(l):
        nom *= (x-xin[i])      #Create variable for nominator
        Inter[i] = yin[i]

    for i in range(l):
        prod = 1               #Variable for the product of the product terms
        for j in range(l):
            if i >= j:
                pass
            else:
                p = (xin[i]-xin[j])
                prod /= p
                Inter[j] /= -p   #ensure each 'p' can be reused in a future sum where possible to reduce calculation time
        Intersum += prod*Inter[i]*nom/(x-xin[i])   #Add term to sum by multiplying product and y-value

    return Intersum


#%%
#1
O = OsloModel(L = 32)
result = O.hAve(1000)
print(result)
O_2 = OsloModel(L = 16)
result2 = O.hAve(1000)
print(result2)
O_3 = OsloModel(L = 32, p = 1)
O_3.Test2()


#%%
#2a
O4 = OsloModel(L = 4)
h4,Tc4 = O4.htdata(1000)
O8 = OsloModel(L = 8)
h8,Tc8 = O8.htdata(1000)
O16 = OsloModel(L = 16)
h16,Tc16 = O16.htdata(1000)
O32 = OsloModel(L = 32)
h32,Tc32 = O32.htdata(1000)
O64 = OsloModel(L = 64)
h64,Tc64 = O64.htdata(1000)
O128 = OsloModel(L = 128)
h128,Tc128 = O128.htdata(1000)
O256 = OsloModel(L = 256)
h256,Tc256 = O256.htdata(1000)

t =np.arange(0,1001,1)
#%%
#2a


f,ax = plt.subplots()
ax.plot(t, h4, 'k-', label = 'L = 4')
ax.plot([Tc4,Tc4],[0,60], 'k--', label = 'Start of R for L = 4')
ax.plot(t, h8, 'b-', label = 'L = 8')
ax.plot([Tc8,Tc8],[0,60], 'b--', label = 'Start of R for L = 8')
ax.plot(t, h16, 'g-', label = 'L = 16')
ax.plot([Tc16,Tc16],[0,60], 'g--', label = 'Start of R for L = 16')
ax.plot(t, h32, 'r-', label = 'L = 32')
ax.plot([Tc32,Tc32],[0,60], 'r--', label = 'Start of R for L = 32')
ax.plot(t, h64, '-', label = 'L=64')
ax.plot(t, h128, 'm-', label = 'L = 128')
ax.plot(t, h256, '-', label = 'L = 256')
ax.set(xlabel='Time', ylabel='Height', title = 'Height versus Time for various system sizes')
ax.grid()
ax.legend(loc='upper left', ncol=2, prop={'size': 7})
plt.show

#%%
#2B


Tc4Ave = O4.TcAve(100)
Tc8Ave = O8.TcAve(100)
Tc16Ave = O16.TcAve(100)
Tc32Ave = O32.TcAve(100)
Tc64Ave = O64.TcAve(100)
Tc128Ave = O128.TcAve(100)
Tc256Ave = O256.TcAve(100)

#%%

O4 = OsloModel(L = 4)
H4 = O4.htdataAve(1000, 1000)
O8 = OsloModel(L = 8)
H8 = O8.htdataAve(1000, 1000)
O16 = OsloModel(L = 16)
H16 = O16.htdataAve(1000, 1000)
O32 = OsloModel(L = 32)
H32 = O32.htdataAve(1000, 1000)
O64 = OsloModel(L = 64)
H64 = O64.htdataAve(1000, 1000)
O128 = OsloModel(L = 128)
H128 = O128.htdataAve(1000, 1000)
O256 = OsloModel(L = 256)
H256 = O256.htdataAve(1000, 1000)



#%%


H4ave = O4.hAve(100000)
H8ave = O8.hAve(100000)
H16ave = O16.hAve(100000)
H32ave = O32.hAve(100000)
H64ave = O64.hAve(100000)
H128ave = O128.hAve(100000)
H256ave = O256.hAve(100000)
#%%

f,ax = plt.subplots()
ax.plot(t/Tc4Ave, H4/H4ave, 'k-', label = 'L = 4')
ax.plot(t/Tc8Ave, H8/H8ave, 'b-', label = 'L = 8')
ax.plot(t/Tc16Ave, H16/H16ave, 'g-', label = 'L = 16')
ax.plot(t/Tc32Ave, H32/H32ave, 'r-', label = 'L = 32')
ax.plot(t/Tc64Ave, H64/H64ave, '-', label = 'L=64')
ax.plot(t/Tc128Ave, H128/H128ave, 'm-', label = 'L = 128')
ax.plot(t/Tc256Ave, H256/H256ave, '-', label = 'L = 256')
ax.set(xlabel=r'$t/{\langle t_c \rangle}$', ylabel=r'$h/\langle h\rangle$', title = 'Rescaled Height versus Time Data Collapse', xlim = (-1,20))
ax.grid()
ax.legend(loc='lower right', ncol=2, prop={'size': 10})
plt.show


#%%

O4 = OsloModel(L = 4)
O8 = OsloModel(L = 8)
O16 = OsloModel(L = 16)
O32 = OsloModel(L = 32)
O64 = OsloModel(L = 64)
O128 = OsloModel(L = 128)
O256 = OsloModel(L = 256)


#%%
print(H4ave)
print(H8ave)
print(H16ave)
print(H32ave)
print(H64ave)
print(H128ave)
print(H256ave)
#%%

H4s = (O4.hVarAve(10000))**0.5
H8s = (O8.hVarAve(10000))**0.5
H16s = (O16.hVarAve(10000))**0.5
H32s = (O32.hVarAve(10000))**0.5
H64s = (O64.hVarAve(10000))**0.5
H128s = (O128.hVarAve(10000))**0.5
H256s = (O256.hVarAve(10000))**0.5
#%%
print(H4s)
print(H8s)
print(H16s)
print(H32s)
print(H64s)
print(H128s)
print(H256s)


#%%
import sys
sys.path.append("/Users/georgemercieca/OneDrive - Imperial College London/Year 3 Physics - Imperial College London/Year 3 Comp_Phys/Assignment")
import CodeModule_copy as C
Hs = [H4s,H8s,H16s,H32s,H64s,H128s,H256s]
LArr = [4,8,16,32,64,128,256]
sqrtLArr = [2, 8**0.5 ,4 , 32**0.5, 8, 128**0.5, 16]
xt = np.arange(0,16,0.1)

fix,ax1 = plt.subplots()
ax1.plot(sqrtLArr, Hs, 'kx', label = 'Recorded Stadard Deviations')
#ax1.plot(xt, C.InterpolatePoly(xt, Hs, L), 'g', label = 'Lagrange')
ax1.plot(xt, C.CubicSplineArr(xt, Hs, sqrtLArr),'b-', label='Interpolation')

ax1.set(xlabel=r'$\sqrt{L}$', ylabel= (r'$\sigma$'), title = 'Interpolation of Standard Deviations')
ax1.grid()
ax1.legend(loc='upper right', ncol=2)
plt.show


#%%
h4data, ph4data = O4.hProbdata(100000)
h8data, ph8data = O8.hProbdata(100000)
h16data, ph16data = O16.hProbdata(100000)
h32data, ph32data = O32.hProbdata(100000)
h64data, ph64data = O64.hProbdata(100000)
h128data, ph128data = O128.hProbdata(100000)
h256data, ph256data = O256.hProbdata(100000)
#%%
f,ax2 = plt.subplots()
ax2.plot(h4data, ph4data, 'kx', label = 'L = 4')
ax2.plot(h8data, ph8data, 'bx', label = 'L = 8')
ax2.plot(h16data, ph16data, 'gx', label = 'L = 16')
ax2.plot(h32data, ph32data, 'rx', label = 'L = 32')
ax2.plot(h64data, ph64data, 'x', label = 'L=64')
ax2.plot(h128data, ph128data, 'mx', label = 'L = 128')
ax2.plot(h256data, ph256data, 'x', label = 'L = 256')
ax2.set(xlabel='h', ylabel='P(h;L)', title = 'Probabilities of h for each system')
ax2.grid()
ax2.legend(loc='upper right', ncol=2, prop={'size': 10})
plt.show
#%%
f,ax3 = plt.subplots()
ax3.plot((h4data-H4ave)/H4s, ph4data*H4s, 'kx', label = 'L = 4')
ax3.plot((h8data-H8ave)/H8s, ph8data*H8s, 'bx', label = 'L = 8')
ax3.plot((h32data-H32ave)/H32s, ph32data*H32s, 'rx', label = 'L = 16')
ax3.plot((h16data-H16ave)/H16s, ph16data*H16s, 'gx', label = 'L = 32')
ax3.plot((h64data-H64ave)/H64s, ph64data*H64s, 'x', label = 'L=64')
ax3.plot((h128data-H128ave)/H128s, ph128data*H128s, 'mx', label = 'L = 128')
ax3.plot((h256data-H256ave)/H256s, ph256data*H256s, 'x', label = 'L = 256')
ax3.set(xlabel=r'$(h-\langle h\rangle)/\sigma_h$', ylabel=r'$\sigma_h P(h;L)$', title = 'Probabilities of h for each system data collapse')
ax3.grid(True)
ax3.legend(loc='upper right', ncol=2, prop={'size': 10})
plt.show
#%%

z4data, pz4data = O4.zdata(64000)
z8data, pz8data = O8.zdata(32000)
z16data, pz16data = O16.zdata(16000)
z32data, pz32data = O32.zdata(8000)
z64data, pz64data = O64.zdata(4000)
z128data, pz128data = O128.zdata(2000)
z256data, pz256data = O256.zdata(1000)

print(pz4data)
print(pz8data)
print(pz16data)
print(pz32data)
print(pz64data)
print(pz128data)
print(pz256data)

#%%
start = time.time()
s4data = O4.sdata(100000)
end = time.time()
runtime = end - start
print(runtime)
s8data = O8.sdata(100000)
end = time.time()
runtime = end - start
print(runtime)
s16data = O16.sdata(100000)
end = time.time()
runtime = end - start
print(runtime)
s32data = O32.sdata(100000)
end = time.time()
runtime = end - start
print(runtime)
s64data = O64.sdata(100000)
end = time.time()
runtime = end - start
print(runtime)
s128data = O128.sdata(100000)
end = time.time()
runtime = end - start
print(runtime)
s256data = O256.sdata(100000)
end = time.time()
runtime = end - start
print(runtime)

#%%
s4bins, s4binFreq = L.logbin(s4data, scale = 1.1, zeros = True)
s8bins, s8binFreq = L.logbin(s8data, scale = 1.1, zeros = True)
s16bins, s16binFreq = L.logbin(s16data, scale = 1.1, zeros = True)
s32bins, s32binFreq = L.logbin(s32data, scale = 1.1, zeros = True)
s64bins, s64binFreq = L.logbin(s64data, scale = 1.1, zeros = True)
s128bins, s128binFreq = L.logbin(s128data, scale = 1.1, zeros = True)
s256bins, s256binFreq = L.logbin(s256data, scale = 1.1, zeros = True)

#%%
f,ax4 = plt.subplots()
ax4.plot(s4bins, s4binFreq, 'k--', label = 'L = 4')
ax4.plot(s8bins, s8binFreq, 'b--', label = 'L = 8')
ax4.plot(s16bins, s16binFreq, 'r--', label = 'L = 16')
ax4.plot(s32bins, s32binFreq, 'g--', label = 'L = 32')
ax4.plot(s64bins, s64binFreq, '--', label = 'L=64')
ax4.plot(s128bins, s128binFreq, 'm--', label = 'L = 128')
ax4.plot(s256bins, s256binFreq, '--', label = 'L = 256')
ax4.set(xlabel='s', ylabel=r'$\tildeP_N (s;L)$', title = 'Probabilities of s for each system')
ax4.grid(True)
ax4.legend(loc='upper right', ncol=2, prop={'size': 10})
ax4.set_xscale('log')
ax4.set_yscale('log')
plt.show

#%%
D = 2.25
ts = 1.55
s4binsdc = s4bins/(O4._L**D)
s8binsdc = s8bins/(O8._L**D)
s16binsdc = s16bins/(O16._L**D)
s32binsdc = s32bins/(O32._L**D)
s64binsdc = s64bins/(O64._L**D)
s128binsdc = s128bins/(O128._L**D)
s256binsdc = s256bins/(O256._L**D)

def binFreqtobinFreqdc(binFreq, bins, ts):
    """
    Function to save code and make it easier
    to transform the binFreq arrays to perform a data collapse
    Needs the binFreq array, its respective bins array and the critcial
    exponent \tau_s
    """
    l = len(binFreq)
    binFreqdc = np.zeros(l)
    for i in range(l):
        binFreqdc[i] = binFreq[i]*(bins[i]**ts)
    return binFreqdc

s4binFreqdc = binFreqtobinFreqdc(s4binFreq, s4bins, ts)
s8binFreqdc = binFreqtobinFreqdc(s8binFreq, s8bins, ts)
s16binFreqdc = binFreqtobinFreqdc(s16binFreq, s16bins, ts)
s32binFreqdc = binFreqtobinFreqdc(s32binFreq, s32bins, ts)
s64binFreqdc = binFreqtobinFreqdc(s64binFreq, s64bins, ts)
s128binFreqdc = binFreqtobinFreqdc(s128binFreq, s128bins, ts)
s256binFreqdc = binFreqtobinFreqdc(s256binFreq, s256bins, ts)



f,ax4 = plt.subplots()
ax4.plot(s4binsdc, s4binFreqdc, 'k--', label = 'L = 4')
ax4.plot(s8binsdc, s8binFreqdc, 'b--', label = 'L = 8')
ax4.plot(s16binsdc, s16binFreqdc, 'r--', label = 'L = 16')
ax4.plot(s32binsdc, s32binFreqdc, 'g--', label = 'L = 32')
ax4.plot(s64binsdc, s64binFreqdc, '--', label = 'L=64')
ax4.plot(s128binsdc, s128binFreqdc, 'm--', label = 'L = 128')
ax4.plot(s256binsdc, s256binFreqdc, '--', label = 'L = 256')
ax4.set(xlabel=r'$s/L^D$', ylabel=r'$s^{\tau_s} \tildeP_N (s;L)$', title = 'Probabilities of s for each system data collapse')
ax4.grid(True)
ax4.legend(loc='upper right', ncol=2, prop={'size': 10})
ax4.set_xscale('log')
ax4.set_yscale('log')
plt.show

#%%

def kthmoment(sdata, k):
    """
    Function for determining the kth moment of an avalanche size
    manually if provided with data on T avalanches in array form
    during the steady state. Returns 
    """
    T = len(sdata)
    sk = 0
    if k == 1:
        return sum(sdata)/T
    elif k > 1:
        for i in range(T):
            s = sdata[i]**k
            if s < 0:
                s = float(sdata[i])
                s = s**4
            sk += s
        sk /= T
        return sk
#%%

s4k1 = kthmoment(s4data, 1)
s4k2 = kthmoment(s4data, 2)
s4k3 = kthmoment(s4data, 3)
s4k4 = kthmoment(s4data, 4)

s8k1 = kthmoment(s8data, 1)
s8k2 = kthmoment(s8data, 2)
s8k3 = kthmoment(s8data, 3)
s8k4 = kthmoment(s8data, 4)

s16k1 = kthmoment(s16data, 1)
s16k2 = kthmoment(s16data, 2)
s16k3 = kthmoment(s16data, 3)
s16k4 = kthmoment(s16data, 4)

s32k1 = kthmoment(s32data, 1)
s32k2 = kthmoment(s32data, 2)
s32k3 = kthmoment(s32data, 3)
s32k4 = kthmoment(s32data, 4)

s64k1 = kthmoment(s64data, 1)
s64k2 = kthmoment(s64data, 2)
s64k3 = kthmoment(s64data, 3)
s64k4 = kthmoment(s64data, 4)

s128k1 = kthmoment(s128data, 1)    
s128k2 = kthmoment(s128data, 2)    
s128k3 = kthmoment(s128data, 3)    
s128k4 = kthmoment(s128data, 4)    

s256k1 = kthmoment(s256data, 1)    
s256k2 = kthmoment(s256data, 2)    
s256k3 = kthmoment(s256data, 3)    
s256k4 = kthmoment(s256data, 4)   
#%%

sk1arr = [s4k1,s8k1,s16k1,s32k1,s64k1,s128k1,s256k1]
sk2arr = [s4k2,s8k2,s16k2,s32k2,s64k2,s128k2,s256k2]
sk3arr = [s4k3,s8k3,s16k3,s32k3,s64k3,s128k3,s256k3]
sk4arr = [s4k4,s8k4,s16k4,s32k4,s64k4,s128k4,s256k4]
   
#%%


f,ax5 = plt.subplots()
ax5.set_xscale('log')
ax5.set_yscale('log')
ax5.plot(LArr, sk1arr, 'k-', label = 'k = 1')
ax5.plot(LArr, sk2arr, 'b-', label = 'k = 2')
ax5.plot(LArr, sk3arr, 'r-', label = 'k = 3')
ax5.plot(LArr, sk4arr, 'g-', label = 'k = 4')
ax5.plot(LArr, sk1arr, 'kx')
ax5.plot(LArr, sk2arr, 'bx')
ax5.plot(LArr, sk3arr, 'rx')
ax5.plot(LArr, sk4arr, 'gx')

ax5.set(xlabel='L', ylabel=r'$\langle s^k \rangle$', title = 'Probabilities of s for each system data collapse')
ax5.grid(True)
ax5.legend(loc='upper left', ncol=2, prop={'size': 10})

plt.show

 

#%%
print(1)
