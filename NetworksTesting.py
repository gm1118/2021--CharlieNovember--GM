
import numpy as np
import random as ran
import matplotlib.pyplot as plt
import logbin_2020 as L
import time
import networkx as nx
from scipy.stats import chisquare
import Networks_Code_Module as N

#%%
#PA Test

Freq = np.array([0,0,0,0,0,0])
Exp = np.array([13/54,5/54,11/54,9/54,9/54,7/54])
for i in range(10**6):
    G0 = N.TestG0()    
    G1 = N.BAModel(7, 1, G0 = G0, q = 2/3)
    Freq[(list(nx.edges(G0,nbunch = 6))[0][1])]+=1
Prob = Freq/(10**6)
print(chisquare(Prob,Exp))
#%%
#Figure 2
G0 = N.TestG0()
plt.subplot(121)
nx.draw(G0, with_labels=True, font_weight='bold')
plt.show()
#%%
#Growth Test
for i in range(10**6):
    G1 = N.BAModel(10**6,3)
    if nx.number_of_nodes(G1) != 10**6:
        print(i)
        raise Exception('Fail')
print('Success')

   
#%%        

#start = time.time()
m1data = N.KdistData(10000, 1, 3200)       
m1bins, m1freq = L.logbin(m1data, scale = 1.1, zeros = False)
m1theo = N.TheoDist(1, m1bins)
m3data = N.KdistData(10000, 3, 1600)       
m3bins, m3freq = L.logbin(m3data, scale = 1.1, zeros = False)
m3theo = N.TheoDist(3, m3bins)
m9data = KdistData(10000, 9, 800)       
m9bins, m9freq = L.logbin(m9data, scale = 1.1, zeros = False)
m9theo = N.TheoDist(9, m9bins)
m27data = N.KdistData(10000, 27, 400)       
m27bins, m27freq = L.logbin(m27data, scale = 1.1, zeros = False)
m27theo = N.TheoDist(27, m27bins)
m81data = N.KdistData(10000, 81, 200)       
m81bins, m81freq = L.logbin(m81data, scale = 1.1, zeros = False)
m81theo = N.TheoDist(81, m81bins)
m243data = N.KdistData(10000, 243, 100)       
m243bins, m243freq = L.logbin(m243data, scale = 1.1, zeros = False)
m243theo = N.TheoDist(243, m243bins)
#end = time.time()
#runtime = end-start
#print(runtime)
#%%
#First Figure - No logbinning
m3ks, m3probs = L.logbin(m3data, scale = 1, zeros = True)

f,ax = plt.subplots()
ax.plot(m3bins, m3freq, 'bo', label = 'Log-Binned')
ax.plot(m3ks, m3probs, 'kx', label = 'Raw Data')
ax.set(xlabel=r'$k$', ylabel=r'$p(k)$')
ax.legend(loc='upper right', ncol=2, prop={'size': 10})
ax.set_xscale('log')
ax.set_yscale('log')
plt.show    

#%%
#Second Figure   - Kdist vs. Theory
f,ax2 = plt.subplots()
ax2.plot(m1bins, m1freq, 'bo', label = 'm = 1', markeredgecolor = 'k')
ax2.plot(m1bins, m1theo, 'b--')
ax2.plot(m3bins, m3freq, 'go', label = 'm = 3', markeredgecolor = 'k')
ax2.plot(m3bins, m3theo, 'g--')
ax2.plot(m9bins, m9freq, 'ko', label = 'm = 9', markeredgecolor = 'k')
ax2.plot(m9bins, m9theo, 'k--')
ax2.plot(m27bins, m27freq, 'ro', label = 'm = 27', markeredgecolor = 'k')
ax2.plot(m27bins, m27theo, 'r--')
ax2.plot(m81bins, m81freq, 'co', label = 'm = 81', markeredgecolor = 'k')
ax2.plot(m81bins, m81theo, 'c--')
ax2.plot(m243bins, m243freq, 'mo', label = 'm = 243', markeredgecolor = 'k')
ax2.plot(m243bins, m243theo, 'm--')
ax2.set(xlabel=r'$k$', ylabel=r'$p(k)$', xlim = (1,10**4))
ax2.legend(loc='lower left', ncol=2, prop={'size': 10})
ax2.set_xscale('log')
ax2.set_yscale('log')
plt.show   
#%%
print(N.R2Test(m1theo, m1freq))
print(N.R2Test(m3theo, m3freq))
print(N.R2Test(m9theo, m9freq))
print(N.R2Test(m27theo, m27freq))
print(N.R2Test(m81theo, m81freq))
print(N.R2Test(m243theo, m243freq))
#%%
#Third Figure - K1 vs. Theory
#start = time.time()
N2data, N2err = N.largestkdata(10**2, 3, 10**6)
N3data, N3err = N.largestkdata(10**3, 3, 10**5)
N4data, N4err = N.largestkdata(10**4, 3, 10**4)
N5data, N5err = N.largestkdata(10**5, 3, 10**3)
N6data, N6err = N.largestkdata(10**6, 3, 10**2)
#end = time.time()
#runtime = end - start
#print(runtime)
#%%
Ns = np.array([10**2,10**3,10**4,10**5,10**6])
fit, fitdata = np.polynomial.polynomial.polyfit(np.log(Ns), np.log(np.array([N2data, N3data, N4data, N5data, N6data])), 1, full = True)
poly1d_fn = np.poly1d(fit)
fiterr = fitdata[0][0]
#%%
f,ax3 = plt.subplots()
ax3.errorbar(Ns, np.array([N2data, N3data, N4data, N5data, N6data]),yerr = np.array([N2err,N3err,N4err,N5err,N6err]), fmt = 'kx', label = 'Numerical Data')
ax3.plot(Ns, [Theok(10**2,3),Theok(10**3,3),Theok(10**4,3),Theok(10**5,3),Theok(10**6,3)], 'k--', label = 'Theoretical Data')
ax3.plot(Ns, np.exp(fit[1]*np.log(Ns)+fit[0]), 'k-', label = 'Fit')
ax3.set(xlabel=r'$N$', ylabel=r'$k_1$')
ax3.legend(loc='upper left', ncol=2, prop={'size': 10})
ax3.set_xscale('log')
ax3.set_yscale('log')
plt.show   

#%%
#start = time.time()
m1kdata, m1err = N.largestkdata(10**4, 1, 3200)
m3kdata, m3err = N.largestkdata(10**4, 3, 1600)
m9kdata, m9err = N.largestkdata(10**4, 9, 800)
m27kdata, m27err = N.largestkdata(10**4, 27, 400)
m81kdata, m81err = N.largestkdata(10**4, 81, 200)
m243kdata, m243err = N.largestkdata(10**4, 243, 100)
#end = time.time()
#runtime = end - start
#print(runtime)
#%%
m1ploty = [x/y for x,y in zip(m1freq,m1theo)]
m3ploty = [x/y for x,y in zip(m3freq,m3theo)]
m9ploty = [x/y for x,y in zip(m9freq,m9theo)]
m27ploty = [x/y for x,y in zip(m27freq,m27theo)]
m81ploty = [x/y for x,y in zip(m81freq,m81theo)]
m243ploty = [x/y for x,y in zip(m243freq,m243theo)]
m1plotx = [x/m1kdata for x in m1bins]
m3plotx = [x/m3kdata for x in m3bins]
m9plotx = [x/m9kdata for x in m9bins]
m27plotx = [x/m27kdata for x in m27bins]
m81plotx = [x/m81kdata for x in m81bins]
m243plotx = [x/m243kdata for x in m243bins]
#%%
#Figure 4 - Data Collapse (various m, unused in report)
f,ax4 = plt.subplots()
ax4.plot(m1plotx, m1ploty, 'bo', label = 'm = 1', markeredgecolor = 'k')
ax4.plot(m3plotx, m3ploty, 'go', label = 'm = 3', markeredgecolor = 'k')
ax4.plot(m9plotx, m9ploty, 'ko', label = 'm = 9', markeredgecolor = 'k')
ax4.plot(m27plotx, m27ploty, 'ro', label = 'm = 27', markeredgecolor = 'k')
ax4.plot(m81plotx, m81ploty, 'co', label = 'm = 81', markeredgecolor = 'k')
ax4.plot(m243plotx, m243ploty, 'mo', label = 'm = 243', markeredgecolor = 'k')
ax4.set(xlabel=r'$k/k_1$', ylabel=r'$p(k)/p_\infty (k)$')
ax4.legend(loc='lower left', ncol=2, prop={'size': 10})
ax4.set_xscale('log')
ax4.set_yscale('log')
plt.show   
#%%
N2distdata = N.KdistData(10**2, 3, 10**5)       
N2bins, N2freq = L.logbin(N2distdata, scale = 1.1, zeros = False)
N2theo = N.TheoDist(3, N2bins)
N3distdata = N.KdistData(10**3, 3, 10**4)       
N3bins, N3freq = L.logbin(N3distdata, scale = 1.1, zeros = False)
N3theo = N.TheoDist(3, N3bins)
N4distdata = N.KdistData(10**4, 3, 10**3)       
N4bins, N4freq = L.logbin(N4distdata, scale = 1.1, zeros = False)
N4theo = N.TheoDist(3, N4bins)
N5distdata = N.KdistData(10**5, 3, 10**2)       
N5bins, N5freq = L.logbin(N5distdata, scale = 1.1, zeros = False)
N5theo = N.TheoDist(3, N5bins)
N6distdata = N.KdistData(10**6, 3, 10**1)       
N6bins, N6freq = L.logbin(N6distdata, scale = 1.1, zeros = False)
N6theo = N.TheoDist(3, N6bins)

#%%
N2ploty = [x/y for x,y in zip(N2freq,N2theo)]
N3ploty = [x/y for x,y in zip(N3freq,N3theo)]
N4ploty = [x/y for x,y in zip(N4freq,N4theo)]
N5ploty = [x/y for x,y in zip(N5freq,N5theo)]
N6ploty = [x/y for x,y in zip(N6freq,N6theo)]

N2plotx = [x/N2data for x in N2bins]
N3plotx = [x/N3data for x in N3bins]
N4plotx = [x/N4data for x in N4bins]
N5plotx = [x/N5data for x in N5bins]
N6plotx = [x/N6data for x in N6bins]





#%%
#Figure 4 - Data Collapse (various N, used in report)
f,ax4_5 = plt.subplots()
ax4_5.plot(N2plotx, N2ploty, 'bo', label = 'N = 2', markeredgecolor = 'k')
ax4_5.plot(N3plotx, N3ploty, 'go', label = 'N = 3', markeredgecolor = 'k')
ax4_5.plot(N4plotx, N4ploty, 'ko', label = 'N = 4', markeredgecolor = 'k')
ax4_5.plot(N5plotx, N5ploty, 'ro', label = 'N = 5', markeredgecolor = 'k')
ax4_5.plot(N6plotx, N6ploty, 'co', label = 'N = 6', markeredgecolor = 'k')
ax4_5.set(xlabel=r'$k/k_1$', ylabel=r'$p(k)/p_\infty (k)$')
ax4_5.legend(loc='lower left', ncol=2, prop={'size': 10})
ax4_5.set_xscale('log')
ax4_5.set_yscale('log')
plt.show   
#%%
#Figure 5 - RA K-Distributions
#start = time.time()
ram1data = N.KdistData(10000, 1, 3200, q = 0)       
ram1bins, ram1freq = L.logbin(ram1data, scale = 1.1, zeros = False)
ram1theo = N.raTheoDist(1, ram1bins)
ram3data = N.KdistData(10000, 3, 1600, q = 0)       
ram3bins, ram3freq = L.logbin(ram3data, scale = 1.1, zeros = False)
ram3theo = N.raTheoDist(3, ram3bins)
ram9data = N.KdistData(10000, 9, 800, q = 0)       
ram9bins, ram9freq = L.logbin(ram9data, scale = 1.1, zeros = False)
ram9theo = N.raTheoDist(9, ram9bins)
ram27data = N.KdistData(10000, 27, 400, q = 0)       
ram27bins, ram27freq = L.logbin(ram27data, scale = 1.1, zeros = False)
ram27theo = N.raTheoDist(27, ram27bins)
ram81data = N.KdistData(10000, 81, 200, q = 0)       
ram81bins, ram81freq = L.logbin(ram81data, scale = 1.1, zeros = False)
#ram81theo = N.raTheoDist(81, ram81bins)   Cause Error - Cannot be calculated
ram243data = N.KdistData(10000, 243, 100, q = 0)       
ram243bins, ram243freq = L.logbin(ram243data, scale = 1.1, zeros = False)
#ram243theo = N.raTheoDist(243, ram243bins)  #Cause Error - Cannot be calculated
#end = time.time()
#runtime = end - start
#print(runtime)

#%%
f,ax5 = plt.subplots()
ax5.plot(ram1bins, ram1freq, 'bo', label = 'm = 1', markeredgecolor = 'k')
ax5.plot(ram1bins, ram1theo, 'b--')
ax5.plot(ram3bins, ram3freq, 'go', label = 'm = 3', markeredgecolor = 'k')
ax5.plot(ram3bins, ram3theo, 'g--')
ax5.plot(ram9bins, ram9freq, 'ko', label = 'm = 9', markeredgecolor = 'k')
ax5.plot(ram9bins, ram9theo, 'k--')
ax5.plot(ram27bins, ram27freq, 'ro', label = 'm = 27', markeredgecolor = 'k')
ax5.plot(ram27bins, ram27theo, 'r--')
ax5.plot(ram81bins, ram81freq, 'co', label = 'm = 81', markeredgecolor = 'k')
#ax5.plot(ram81bins, ram81theo, 'c--')   #Produces Nothing
ax5.plot(ram243bins, ram243freq, 'mo', label = 'm = 243', markeredgecolor = 'k')
#ax5.plot(ram243bins, ram243theo, 'm--')  #Produces Nothing
ax5.set(xlabel=r'$k$', ylabel=r'$p(k)$', xlim = (1,10**4))
ax5.legend(loc='lower left', ncol=2, prop={'size': 10})
ax5.set_xscale('log')
ax5.set_yscale('log')
plt.show

#%%
# Figure 6 - RA K1 values
start = time.time()
raN2data, raN2err = N.largestkdata(10**2, 3, 10**6, q = 0)
raN3data, raN3err = N.largestkdata(10**3, 3, 10**5, q = 0)
raN4data, raN4err = N.largestkdata(10**4, 3, 10**4, q = 0)
raN5data, raN5err = N.largestkdata(10**5, 3, 10**3, q = 0)
raN6data, raN6err = N.largestkdata(10**6, 3, 10**2, q = 0)
#end = time.time()
#runtime - end-start
#print(runtime)
#%%
Ns = np.array([10**2,10**3,10**4,10**5,10**6])
#rafit, rafitdata = np.polynomial.polynomial.polyfit(np.log(Ns), np.log(np.array([raN2data, raN3data, raN4data, raN5data, raN6data])), 1, full = True)
#rapoly1d_fn = np.poly1d(rafit)
#rafiterr = rafitdata[0][0]
Ratheo_x = np.arange(1.8,6.2,0.1)
ratheo_x = [10**x for x in Ratheo_x]
ratheo_y = [N.raTheok(x,3) for x in ratheo_x ]
#%%
f,ax6 = plt.subplots()
ax6.errorbar(Ns, np.array([raN2data, raN3data, raN4data, raN5data, raN6data]),yerr = np.array([raN2err,raN3err,raN4err,raN5err,raN6err]), fmt = 'kx', label = 'Numerical Data')
ax6.plot(ratheo_x, ratheo_y, 'k--', label = 'Theoretical Data')
#ax6.plot(Ns, np.exp(rafit[1]*np.log(Ns)+rafit[0]), 'k-', label = 'Fit')
ax6.set(xlabel=r'$N$', ylabel=r'$k_1$')
ax6.legend(loc='upper left', ncol=2, prop={'size': 10})
ax6.set_xscale('log')
ax6.set_yscale('log')
plt.show  
#%%
# Figure 7 - MA - K Distribution
#start = time.time()
mam1data = N.KdistData(10000, 1, 3200, q = 2/3)       
mam1bins, mam1freq = L.logbin(mam1data, scale = 1.1, zeros = False)
mam1theo = N.maTheoDist(1, mam1bins)
mam3data = N.KdistData(10000, 3, 1600, q = 2/3)       
mam3bins, mam3freq = L.logbin(mam3data, scale = 1.1, zeros = False)
mam3theo = N.maTheoDist(3, mam3bins)
mam9data = N.KdistData(10000, 9, 800, q = 2/3)       
mam9bins, mam9freq = L.logbin(mam9data, scale = 1.1, zeros = False)
mam9theo = N.maTheoDist(9, mam9bins)
mam27data = N.KdistData(10000, 27, 400, q = 2/3)       
mam27bins, mam27freq = L.logbin(mam27data, scale = 1.1, zeros = False)
mam27theo = N.maTheoDist(27, mam27bins)
mam81data = N.KdistData(10000, 81, 200, q = 2/3)       
mam81bins, mam81freq = L.logbin(mam81data, scale = 1.1, zeros = False)
mam81theo = N.maTheoDist(81, mam81bins)   #Calculate Manually
mam243data = N.KdistData(10000, 243, 100, q = 2/3)       
mam243bins, mam243freq = L.logbin(mam243data, scale = 1.1, zeros = False)
mam243theo = N.maTheoDist(243, mam243bins)  #Calculate manually
#end = time.time()
#runtime = end-start
#print(runtime)
#%%
f,ax7 = plt.subplots()
ax7.plot(mam1bins, mam1freq, 'bo', label = 'm = 1', markeredgecolor = 'k')
ax7.plot(mam1bins, mam1theo, 'b--')
ax7.plot(mam3bins, mam3freq, 'go', label = 'm = 3', markeredgecolor = 'k')
ax7.plot(mam3bins, mam3theo, 'g--')
ax7.plot(mam9bins, mam9freq, 'ko', label = 'm = 9', markeredgecolor = 'k')
ax7.plot(mam9bins, mam9theo, 'k--')
ax7.plot(mam27bins, mam27freq, 'ro', label = 'm = 27', markeredgecolor = 'k')
ax7.plot(mam27bins, mam27theo, 'r--')
ax7.plot(mam81bins, mam81freq, 'co', label = 'm = 81', markeredgecolor = 'k')
ax7.plot(mam81bins, mam81theo, 'c--')
ax7.plot(mam243bins, mam243freq, 'mo', label = 'm = 243', markeredgecolor = 'k')
ax7.plot(mam243bins, mam243theo, 'm--')
ax7.set(xlabel=r'$k$', ylabel=r'$p(k)$', xlim = (1,10**4))
ax7.legend(loc='lower left', ncol=2, prop={'size': 10})
ax7.set_xscale('log')
ax7.set_yscale('log')
plt.show

#%%
# Figure 8 - MA K1 values
#start = time.time()
maN2data, maN2err = N.largestkdata(10**2, 3, 10**6, q = 2/3)
maN3data, maN3err = N.largestkdata(10**3, 3, 10**5, q = 2/3)
maN4data, maN4err = N.largestkdata(10**4, 3, 10**4, q = 2/3)
maN5data, maN5err = N.largestkdata(10**5, 3, 10**3, q = 2/3)
maN6data, maN6err = N.largestkdata(10**6, 3, 10**2, q = 2/3)
#end = time.time()
#runtime = end-start
#print(runtime)
#%%
Ns = np.array([10**2,10**3,10**4,10**5,10**6])
mafit, mafitdata = np.polynomial.polynomial.polyfit(np.log(Ns), np.log(np.array([maN2data, maN3data, maN4data, maN5data, maN6data])), 1, full = True)
mapoly1d_fn = np.poly1d(mafit)
mafiterr = mafitdata[0][0]
#Ratheo_x = np.arange(1.8,6.2,0.1)
#ratheo_x = [10**x for x in Ratheo_x]
#ratheo_y = [N.raTheok(x,3) for x in ratheo_x ]
#%%
f,ax8 = plt.subplots()
ax8.errorbar(Ns, np.array([maN2data, maN3data, maN4data, maN5data, maN6data]),yerr = np.array([maN2err,maN3err,maN4err,maN5err,maN6err]), fmt = 'kx', label = 'Numerical Data')
#ax8.plot(ratheo_x, ratheo_y, 'k--', label = 'Theoretical Data')
ax8.plot(Ns, np.exp(mafit[1]*np.log(Ns)+mafit[0]), 'k-', label = 'Fit')
ax8.set(xlabel=r'$N$', ylabel=r'$k_1$')
ax8.legend(loc='upper left', ncol=2, prop={'size': 10})
ax8.set_xscale('log')
ax8.set_yscale('log')
plt.show  

