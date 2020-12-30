import numpy as np
import scipy.stats as sci
import matplotlib.pyplot as plt
import time
import math
import os

def driftmodeling(flynum, numberofbins, numberofdays, prefmean, prefvariance, envimean, envivariance, driftvariance, adaptivetracking, gain, per, maxsurvivalrate, birthrate, matureage, percentbh, showgraphs, figuresavepath):
    # adaptivetracking=0
    x=np.linspace(-1,1,numberofbins) # Make a certain number of bins for numbers between -1 and 1
    maxage=30 # Set maximum age
    prefvariance=np.array([prefvariance]) # Make an array for prefvariance
    prefmean=np.array([prefmean]) # Make an array for prefmean
    driftvariance=np.array([driftvariance]) # Make an array for driftvariance
    adaptivetracking=np.array([adaptivetracking]) # Make an array for adaptivetracking
    numconditions=max(prefvariance.shape) # Set the number of conditions based on the total conditions we're running
    finalpop=np.zeros((numconditions)) # Set a zeros array for the finalpop

    for q in range(numconditions): # Run a for loop for each condition
        pref=np.zeros((numberofbins,numberofdays,maxage)) # Make a zeros matrix for [bins, days, maxage]
        if prefvariance[q]>=0.015:  # Check if variance is too small so that we can just eliminate bet-hedging
            pref[:,0,0]=sci.norm.pdf(x,prefmean[q],prefvariance[q]) # If variance is large enough, set a fly's first day preference gaussian with center around 0
        else: # If not, make the bin in the middle have all the flies
            pref[math.floor(numberofbins/2),0,0]=flynum
        pref[:,0,0]=pref[:,0,0]/np.sum(pref[:,0,0])*flynum # Set the total number of flies to equal flynum

        if prefvariance[q]*percentbh>=0.015: # Check if variance is too small so that we can just eliminate bet-hedging
            reducedbethedgeinitial=sci.norm.pdf(x,prefmean[q],np.multiply(prefvariance[q],percentbh)) # If it is large enough, set the reduced bh fly's first day preference gaussian with center around 0 and variance reduced
        else: # If not, make the bin in the middle have all the flies
            reducedbethedgeinitial=np.zeros(numberofbins)
            reducedbethedgeinitial[math.floor(numberofbins/2)]=flynum
        reducedbethedgeinitial[:]=reducedbethedgeinitial[:]/np.sum(reducedbethedgeinitial[:])*flynum # Set the total number of reduced bh flies to equal flynum

        envi=np.zeros((numberofbins,numberofdays)) # Make a zeros array for environment
        envi[:,0]=sci.norm.pdf(x,envimean,envivariance) # Make a gaussian of environment with center around 0
        envi=envi/(np.max(envi))*maxsurvivalrate # Normalizing the maximum envi value and factoring in deathrate
        driftadvantage=np.zeros((numberofdays)) # Make a zeros array for driftadvantage
        betadvantage=np.zeros((numberofdays)) # Make a zeros array for betadvantage
        blur=np.zeros((numberofbins,numberofbins)) # Make a zeros array for blur [which bin profile it's for, what the distribution is between that bin and all other bins, 2]

        for b in range(numberofbins): # Calculate a value to blur the bins
            blur[b,:]=sci.norm.pdf(x,x[b],driftvariance[q])

        for t in range(1,numberofdays): # For each day...

            envi[:,t]=sci.norm.pdf(x,(envimean+gain*np.sin(t*np.pi*2/per)),envivariance) # Making envi a sin wave that changes over time
            envi[:,t]=envi[:,t]/np.max(envi[:,t])*maxsurvivalrate # Normalizing the envi and multiplying by maxsurvival rate
            pref[:,t,0]=pref[:,0,0]*birthrate/flynum*np.sum(pref[:,t-1,matureage:]) # Calculate newborn flies
            if adaptivetracking[q]>0:
                pref[:,t,0]=pref[:,t,0]*(1-adaptivetracking[q])+adaptivetracking[q]*birthrate*np.sum(pref[:,t-1,matureage:],1)

                #maybe we should consider putting in some amount of variation on adaptivetracking (shift mean but keep bet hedging?)
            numfliesborntoday=np.sum(pref[:,t,0]) # Calculate the new number of flies born on that day

            betadvantage[t]=np.sum(np.multiply(pref[:,t,0], envi[:,t]))-numfliesborntoday*envi[math.floor(numberofbins/2),t] #
            
            print(np.sum(blur[b,:]))

            for a in range(maxage): # For the flies age...

                driftadvantage[t]+=np.sum(np.multiply(pref[:,t-1,a-1], envi[:,t])) # Calculating the number of flies that survive without drift #Should extend to include BH

                if a>0:
                    for b in range(numberofbins):
                        if driftvariance[q]<.05/numberofbins: # Check if blur is too low to be worth blurring. NOTE: This number should be based on the limit of sci.norm.pdf
                            pref[b,t,a]+=pref[b,t-1,a-1] # If it is too low, just carry it forward
                        else: # If not, multiply by the blur value and carry it forward
                            pref[:,t,a]+=pref[b,t-1,a-1]*blur[b,:]/np.sum(blur[b,:])
                    pref[:,t,a]=np.multiply(pref[:,t,a], envi[:,t]) # Multiplying the preference to the environment

            driftadvantage[t]=np.sum(pref[:,t,:])-driftadvantage[t]-numfliesborntoday


        # if showgraphs:
        #     fig, (ax0, ax1, ax2, ax3, ax4) = plt.subplots(5, 1)
        #     fig.set_figwidth(10)
        #     fig.set_figheight(12)
        #     fig.tight_layout()
        #     plt.subplots_adjust(hspace=.6)
        #     c=ax0.pcolormesh(envi)
        #     fig.colorbar(c,ax=ax0)
        #     ax0.set_title('Environment (color is fraction of flies of given pref survive)')
        #     ax0.set_ylabel('Preference')
        #     ax0.set_xlabel('Day')

        #     c=ax1.pcolormesh(np.sum(pref[:,:,:],axis=2))
        #     fig.colorbar(c,ax=ax1)
        #     ax1.set_title('Fly Preference (color is log(num) flies each day)')
        #     ax1.set_ylabel('Preference')
        #     ax1.set_xlabel('Day')

        #     ax2.plot(np.log(np.sum(pref[:,:,:],axis=(0,2))))
        #     ax2.set_title('total log(num) flies)') # lowest value is 0.0001 (prefvariance = 0.01 with percent bh=0.01)
        #     ax2.set_ylabel('log(num) flies)')
        #     ax2.set_xlabel('Day')
        #     ax2.set_xlim(0,numberofdays)

        #     # ax3.plot(driftadvantage)
        #     ax3.plot(driftadvantage/np.sum(pref[:,:,:],axis=(0,2)))
        #     ax3.set_title('Change in death rate due to last day\'s drift ')
        #     ax3.set_ylabel('∆surviving flies/total flies')
        #     ax3.set_xlabel('Day')
        #     ax3.set_xlim(0,numberofdays)

        #     ax4.plot(betadvantage/np.sum(pref[:,:,:],axis=(0,2)))
        #     ax4.set_title('Change in death rate due to last day\'s bethedging ')
        #     ax4.set_ylabel('∆surviving flies/total flies')
        #     ax4.set_xlabel('Day')
        #     ax4.set_xlim(0,numberofdays)

        #     fig.colorbar(c,ax=ax2)
        #     fig.colorbar(c,ax=ax3)
        #     fig.colorbar(c,ax=ax4)

        #     fig.suptitle('Bet-hedge variance: '+str(prefvariance[q])+', Drift variance: '+str(driftvariance[q])+', Adaptive Tracking: '+str(adaptivetracking[q]), y=-.05, fontsize=16)

        #     plt.show()

        #     if os.path.exists(figuresavepath):
        #         fig.savefig(os.path.join(figuresavepath,'bh'+str(prefvariance[q])+'dv'+str(driftvariance[q])+'.png'),bbox_inches='tight', pad_inches=.3)
        #         print(figuresavepath)
        #     else:
        #         print('not saving, no valid path')

        finalpop[q]=np.sum(pref[:,-1,:])

    return finalpop
