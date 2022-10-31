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
import Complexity_Code_Module as C


#%%
#1
O = C.OsloModel(L = 32)
result = O.hAve(1000)
print(result)
O_2 = C.OsloModel(L = 16)
result2 = O.hAve(1000)
print(result2)
O_3 = C.OsloModel(L = 32, p = 1)
O_3.Test2()


#%%
#2a
O4 = C.OsloModel(L = 4)
h4,Tc4 = O4.htdata(1000)
O8 = C.OsloModel(L = 8)
h8,Tc8 = O8.htdata(1000)
O16 = C.OsloModel(L = 16)
h16,Tc16 = O16.htdata(1000)
O32 = C.OsloModel(L = 32)
h32,Tc32 = O32.htdata(1000)
O64 = C.OsloModel(L = 64)
h64,Tc64 = O64.htdata(1000)
O128 = C.OsloModel(L = 128)
h128,Tc128 = O128.htdata(1000)
O256 = C.OsloModel(L = 256)
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

O4 = C.OsloModel(L = 4)
H4 = O4.htdataAve(1000, 1000)
O8 = C.OsloModel(L = 8)
H8 = O8.htdataAve(1000, 1000)
O16 = C.OsloModel(L = 16)
H16 = O16.htdataAve(1000, 1000)
O32 = C.OsloModel(L = 32)
H32 = O32.htdataAve(1000, 1000)
O64 = C.OsloModel(L = 64)
H64 = O64.htdataAve(1000, 1000)
O128 = C.OsloModel(L = 128)
H128 = O128.htdataAve(1000, 1000)
O256 = C.OsloModel(L = 256)
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

O4 = C.OsloModel(L = 4)
O8 = C.OsloModel(L = 8)
O16 = C.OsloModel(L = 16)
O32 = C.OsloModel(L = 32)
O64 = C.OsloModel(L = 64)
O128 = C.OsloModel(L = 128)
O256 = C.OsloModel(L = 256)


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

s4binFreqdc = C.binFreqtobinFreqdc(s4binFreq, s4bins, ts)
s8binFreqdc = C.binFreqtobinFreqdc(s8binFreq, s8bins, ts)
s16binFreqdc = C.binFreqtobinFreqdc(s16binFreq, s16bins, ts)
s32binFreqdc = C.binFreqtobinFreqdc(s32binFreq, s32bins, ts)
s64binFreqdc = C.binFreqtobinFreqdc(s64binFreq, s64bins, ts)
s128binFreqdc = C.binFreqtobinFreqdc(s128binFreq, s128bins, ts)
s256binFreqdc = C.binFreqtobinFreqdc(s256binFreq, s256bins, ts)


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

s4k1 = C.kthmoment(s4data, 1)
s4k2 = C.kthmoment(s4data, 2)
s4k3 = C.kthmoment(s4data, 3)
s4k4 = C.kthmoment(s4data, 4)

s8k1 = C.kthmoment(s8data, 1)
s8k2 = C.kthmoment(s8data, 2)
s8k3 = C.kthmoment(s8data, 3)
s8k4 = C.kthmoment(s8data, 4)

s16k1 = C.kthmoment(s16data, 1)
s16k2 = C.kthmoment(s16data, 2)
s16k3 = C.kthmoment(s16data, 3)
s16k4 = C.kthmoment(s16data, 4)

s32k1 = C.kthmoment(s32data, 1)
s32k2 = C.kthmoment(s32data, 2)
s32k3 = C.kthmoment(s32data, 3)
s32k4 = C.kthmoment(s32data, 4)

s64k1 = C.kthmoment(s64data, 1)
s64k2 = C.kthmoment(s64data, 2)
s64k3 = C.kthmoment(s64data, 3)
s64k4 = C.kthmoment(s64data, 4)

s128k1 = C.kthmoment(s128data, 1)    
s128k2 = C.kthmoment(s128data, 2)    
s128k3 = C.kthmoment(s128data, 3)    
s128k4 = C.kthmoment(s128data, 4)    

s256k1 = C.kthmoment(s256data, 1)    
s256k2 = C.kthmoment(s256data, 2)    
s256k3 = C.kthmoment(s256data, 3)    
s256k4 = C.kthmoment(s256data, 4)   
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

 


