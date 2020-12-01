import numpy as np
import scipy.stats as sci
import matplotlib.pyplot as plt
import time
import math
import os

# 1. Start from scratch and recode
    # General format and rewrite
    # Hide all code except comments and recode from there

# 2. Document every single step (make a figure/visualize variable)
    # Make a flowchart for the program
    # Then go back and recode

def driftmodeling(flynum, numberofbins, numberofdays, prefmean, prefvariance, envimean, envivariance, driftvariance, gain, per, maxsurvivalrate, birthrate, matureage, percentbh, showgraphs, figuresavepath):
    x=np.linspace(-1,1,numberofbins) # number of bins between -1 and 1
    maxage=30 #maximum age
    prefvariance=np.array([prefvariance])
    prefmean=np.array([prefmean])
    driftvariance=np.array([driftvariance])
    numconditions=max(prefvariance.shape) # Number of conditions based on the total conditions we're running
    finalpop=np.zeros((numconditions))
    #print(np.floor(numberofbins/2))

    for q in range(numconditions): # for loop for each condition\
        # print(numconditions)
        # print(prefvariance[q])
        # print(driftvariance[q])
        pref=np.zeros((numberofbins,numberofdays,maxage,2)) # Matrix, [bins, days, maxage, bh vs. reducebh]
        #reducebethedge=np.zeros((numberofbins,numberofdays,maxage,2))

        # Set pref[:,0,0,0], which is the "reduced bet hedge version"

        if prefvariance[q]>=0.015:  #Check if variance is so small to just eliminate bet-hedging
            pref[:,0,0,0]=sci.norm.pdf(x,prefmean[q],prefvariance[q]) # A fly's first day preference gaussian of preference with center around 0
        else: # Make the bin in the middle have all the flies
            #print('Zero bet-hedging')
            #pref[50,0,0,0]=flynum
            pref[math.floor(numberofbins/2),0,0,0]=flynum
            #print(pref[:,0,0,0])

        pref[:,0,0,0]=pref[:,0,0,0]/np.sum(pref[:,0,0,0])*flynum # total # of flies=flynum
        pref[:,1,1,0]=pref[:,0,0,0] # Fly ages to 1, day changes to 1, set the same as initial


        # Now set pref[:,0,0,1], which is the "reduced bet hedge version"
        if prefvariance[q]*percentbh>=0.015: #Check if variance is so small to just eliminate bet-hedging
            pref[:,0,0,1]=sci.norm.pdf(x,prefmean[q],np.multiply(prefvariance[q],percentbh))
        else:
            #print('Also Zero bet-hedging')
            #pref[50,0,0,1]=flynum
            pref[math.floor(numberofbins/2),0,0,1]=flynum
        pref[:,0,0,1]=pref[:,0,0,1]/np.sum(pref[:,0,0,1])*flynum # total # of flies=flynum
        #print(reducebethedge[:,0,0])
        pref[:,1,1,1]=pref[:,0,0,1] # Fly ages to 1, day changes to 1, set the same as initial

        envi=np.zeros((numberofbins,numberofdays))
        envi[:,0]=sci.norm.pdf(x,envimean,envivariance) # A gaussian of environment with center around 0
        envi=envi/(np.max(envi))*maxsurvivalrate # Normalizing the maximum envi value and factoring in deathrate
        driftadvantage=np.zeros((numberofdays))
        betadvantage=np.zeros((numberofdays))
        blur=np.zeros((numberofbins,numberofbins,2)) # [which bin profile it's for, what the distribution is between that bin and all other bins, 2]

        for b in range(numberofbins):
            blur[b,:,0]=sci.norm.pdf(x,x[b],driftvariance[q])
            blur[b,:,1]=sci.norm.pdf(x,x[b],driftvariance[q]) #Made both to compare reducedrift with drift, but never implemented

        for t in range(1,numberofdays):
            # print('t is: '+str(t))
            for w in range(2):
                # print(w)
                # print('hi1')
                pref[:,t,0,w]=pref[:,0,0,w]*birthrate/flynum*np.sum(pref[:,t-1,matureage:,0]) # Calculate newborn flies

                # print('hi2')
                envi[:,t]=sci.norm.pdf(x,(envimean+gain*np.sin(t*np.pi*2/per)),envivariance) # Making envi a sin wave that changes over time
                envi[:,t]=envi[:,t]/np.max(envi[:,t])*maxsurvivalrate # Normalizing the envi and multiplying by maxsurvival rate
                # toc = time.perf_counter()
                # print(toc-tic)

                #hi = time.perf_counter()
                # print('hi3')
                for a in range(maxage):
                    if w==0:
                        driftadvantage[t]+=np.sum(np.multiply(pref[:,t-1,a,0], envi[:,t])) # Calculating the number of flies that survive without drift #Should extend to include BH #FInd bug
                    #print(driftadvantage[t])
                # bye = time.perf_counter()
                # print(bye-hi)

                #one = time.perf_counter()
                #TAKES ABOUT 10X LONGER
                # print('hi4')
                for a in range(maxage):
                    #pref[:,t,a]+=pref[:,t-1,a]
                    #not sure why this line was here!
                    if a>0:
                        # print(pref[b,t-1,a,0])
                        for b in range(numberofbins):
                            if driftvariance[q]<.05/numberofbins: #check if blur is too low to be worth blurring. NOTE: This number should be based on the limit of sci.norm.pdf
                                pref[b,t,a,w]+=pref[b,t-1,a,0] #FInd bug
                            else:
                                pref[:,t,a,w]+=pref[b,t-1,a,0]*blur[b,:,w]/np.sum(blur[b,:,w])
                            #pref[:,t,a]+=pref[b,t-1,a]*sci.norm.pdf(x,x[b],driftvariance)/np.sum(sci.norm.pdf(x,x[b],driftvariance))
                    pref[:,t,a,w]=np.multiply(pref[:,t,a,w], envi[:,t]) # Multiplying the preference to the environment #FInd bug
                # two = time.perf_counter()
                # print(two-one)

            # print('hi5')
            driftadvantage[t]=np.sum(pref[:,t,:,0])-driftadvantage[t] #FInd bug
            betadvantage[t]=np.sum(pref[:,t,0,0]-pref[:,t,0,1])
            pref[:,t,1:,0]=pref[:,t,:-1,0]
        if showgraphs:
            #before = time.perf_counter()
            fig, (ax0, ax1, ax2, ax3, ax4) = plt.subplots(5, 1)
            fig.set_figwidth(10)
            fig.set_figheight(12)
            fig.tight_layout()
            plt.subplots_adjust(hspace=.6)
            c=ax0.pcolormesh(envi)
            fig.colorbar(c,ax=ax0)
            ax0.set_title('Environment (color is fraction of flies of given pref survive)')
            ax0.set_ylabel('Preference')
            ax0.set_xlabel('Day')

            c=ax1.pcolormesh(np.sum(pref[:,:,:,0],axis=2))
            fig.colorbar(c,ax=ax1)
            ax1.set_title('Fly Preference (color is log(num) flies each day)')
            ax1.set_ylabel('Preference')
            ax1.set_xlabel('Day')

            ax2.plot(np.log(np.sum(pref[:,:,:,0],axis=(0,2))))
            ax2.set_title('total log(num) flies)') # lowest value is 0.0001 (prefvariance = 0.01 with percent bh=0.01)
            ax2.set_ylabel('log(num) flies)')
            ax2.set_xlabel('Day')
            ax2.set_xlim(0,numberofdays)

            ax3.plot(driftadvantage/np.sum(pref[:,:,:,0],axis=(0,2)))
            ax3.set_title('Change in death rate due to last day\'s drift ')
            ax3.set_ylabel('∆surviving flies/total flies')
            ax3.set_xlabel('Day')
            ax3.set_xlim(0,numberofdays)

            ax4.plot(betadvantage/np.sum(pref[:,:,:,0],axis=(0,2)))
            ax4.set_title('Change in death rate due to last day\'s bethedging ')
            ax4.set_ylabel('∆surviving flies/total flies')
            ax4.set_xlabel('Day')
            ax4.set_xlim(0,numberofdays)

            fig.colorbar(c,ax=ax2)
            fig.colorbar(c,ax=ax3)
            fig.colorbar(c,ax=ax4)

            fig.suptitle('Bet-hedge variance: '+str(prefvariance[q])+', Drift variance: '+str(driftvariance[q]), y=-.05, fontsize=16)

            plt.show()

            if os.path.exists(figuresavepath):
                fig.savefig(os.path.join(figuresavepath,'bh'+str(prefvariance[q])+'dv'+str(driftvariance[q])+'.png'),bbox_inches='tight', pad_inches=.3)
                print(figuresavepath)
            else:
                print('not saving, no valid path')
        #after = time.perf_counter()
        #print(after-before)

        finalpop[q]=np.sum(pref[:,-1,:,0])

    return finalpop
