#pip install colorednoise
import numpy as np
import scipy.stats as sci
import matplotlib.pyplot as plt
import time
import math
import os
import scipy.stats as stat
import colorednoise as cn



# 1. Start from scratch and recode
    # General format and rewrite
    # Hide all code except comments and recode from there

# 2. Document every single step (make a figure/visualize variable)
    # Make a flowchart for the program
    # Then go back and recode
def sinwaveinput(numberofbins, numberofdays, envimean, envivariance, maxsurvivalrate, gain, per):
    envi=np.zeros((numberofbins,numberofdays))
    x=np.linspace(-1,1,numberofbins)
    envi[:,0]=sci.norm.pdf(x,envimean,envivariance) # A gaussian of environment with center around 0
    envi=envi/(np.max(envi))*maxsurvivalrate # Normalizing the maximum envi value and factoring in deathrate

    for t in range(1,numberofdays):
        envi[:,t]=sci.norm.pdf(x,(envimean+gain*np.sin(t*np.pi*2/per)),envivariance) # Making envi a sin wave that changes over time
        envi[:,t]=envi[:,t]/np.max(envi[:,t])*maxsurvivalrate
    return(envi)

def diffusion(numberofbins, numberofdays, envimean,envivariance, maxsurvivalrate, dailydrift):
    envi=np.zeros((numberofbins,numberofdays))
    x=np.linspace(-1,1,numberofbins)
    envi[:,0]=sci.norm.pdf(x,envimean,envivariance) # A gaussian of environment with center around 0
    envi=envi/(np.max(envi))*maxsurvivalrate # Normalizing the maximum envi value and factoring in deathrate
    blur=np.zeros((numberofbins,numberofbins)) # [which bin profile it's for, what the distribution is between that bin and all other bins, 2]

    for b in range(numberofbins):
        blur[b,:]=sci.norm.pdf(x,x[b],dailydrift)

    for t in range(1,numberofdays):
        for b in range (numberofbins):
            # blur=sci.norm.pdf(x,x[b],dailydrift)
            envi[:,t]+=envi[b,t-1]*blur[b,:]/np.sum(blur[b,:])
        # envi[:,t]=sci.norm.pdf(x,(envimean,envivariance) # Making envi a sin wave that changes over time
        # envi[:,t]=envi[:,t]/np.max(envi[:,t])*maxsurvivalrate
    return(envi)

def randomwalk(numberofbins, numberofdays, envimean,envivariance, maxsurvivalrate, dailydrift):
    envi=np.zeros((numberofbins,numberofdays))
    x=np.linspace(-1,1,numberofbins)
    envi[:,0]=sci.norm.pdf(x,envimean,envivariance) # A gaussian of environment with center around 0
    envi=envi/(np.max(envi))*maxsurvivalrate # Normalizing the maximum envi value and factoring in deathrate
    blur=np.zeros((numberofbins,numberofbins)) # [which bin profile it's for, what the distribution is between that bin and all other bins, 2]

    for t in range(1,numberofdays):
        envimean+=np.random.normal(0,dailydrift)
        envi[:,t]=sci.norm.pdf(x,envimean,envivariance) # Making envi a sin wave that changes over time
        envi[:,t]=envi[:,t]/np.max(envi[:,t])*maxsurvivalrate
    return(envi)

def metropolishastingsdrift(numberofbins, numberofdays, envimeanvariance, envivariance, maxsurvivalrate, dailydrift):
    envi=np.zeros((numberofbins,numberofdays))
    x=np.linspace(-1,1,numberofbins)
    envimean=np.random.normal(0,envimeanvariance)
    envi[:,0]=sci.norm.pdf(x, envimean, envivariance) # A gaussian of environment with center around 0

    for t in range(1, numberofdays):
        pcurrentvalue=stat.norm.pdf(envimean,0,envimeanvariance)
        # proposedvalue=np.random.normal(0,variability)
        proposedvalue=np.random.normal(0,envimeanvariance)*dailydrift+envimean*(1-dailydrift)

        pproposedvalue=stat.norm.pdf(proposedvalue,0,envimeanvariance)
        if pproposedvalue/pcurrentvalue>(1-np.random.rand()):
            # return proposedvalue*percentdrift+currentvalue*(1-percentdrift)
            envimean=proposedvalue
        envi[:,t]=sci.norm.pdf(x, envimean, envivariance) # A gaussian of environment with center around 0
        envi[:,t]=envi[:,t]/np.max(envi[:,t])*maxsurvivalrate
    return(envi)

def whitedrift(numberofbins, numberofdays, envimeanvariance, envivariance, maxsurvivalrate, dailydrift):
    envi=np.zeros((numberofbins,numberofdays))
    x=np.linspace(-1,1,numberofbins)
    # envimean=np.random.normal(0,envimeanvariance)

    for t in range(1, numberofdays):
        envimean=np.random.normal(0,envimeanvariance)
        envi[:,t]=sci.norm.pdf(x, envimean, envivariance) # A gaussian of environment with center around 0
        envi[:,t]=envi[:,t]/np.max(envi[:,t])*maxsurvivalrate

    return(envi)

def powerdrift(numberofbins, numberofdays, envimeanvariance, envivariance, maxsurvivalrate, power):
    # beta = 1 # the exponent
    samples = 100 # number of samples to generate
    y = cn.powerlaw_psd_gaussian(power, numberofdays)*envimeanvariance
    envi=np.zeros((numberofbins,numberofdays))
    x=np.linspace(-1,1,numberofbins)
    # envimean=np.random.normal(0,envimeanvariance)
    # stat.norm.ppf()
    for t in range(numberofdays):
        envi[:,t]=sci.norm.pdf(x, y[t], envivariance) # A gaussian of environment with center around 0
        envi[:,t]=envi[:,t]/np.max(envi[:,t])*maxsurvivalrate    
    return(envi)

# def driftmodeling(flynum, numberofbins, numberofdays, prefmean, prefvariance, envimean, envivariance, driftvariance, adaptivetracking, gain, per, maxsurvivalrate, birthrate, matureage, percentbh, showgraphs, figuresavepath):
    # adaptivetracking=0

def driftmodeling(envi, prefmean, prefvariance, driftvariance, adaptivetracking, birthrate, matureage, percentbh, showgraphs, figuresavepath, driftmaxdistribution=0, saveloc='none'):
    flynum=1
    numberofbins=envi.shape[0]
    numberofdays=envi.shape[1]

    x=np.linspace(-1,1,numberofbins) # number of bins between -1 and 1
    maxage=30 #maximum age
    prefvariance=np.array([prefvariance])
    prefmean=np.array([prefmean])
    driftvariance=np.array([driftvariance])
    adaptivetracking=np.array([adaptivetracking])
    numconditions=max(prefvariance.shape) # Number of conditions based on the total conditions we're running
    finalpop=np.zeros((numconditions))
    numdays=np.array([numberofdays])

    for q in range(numconditions): # for loop for each condition\
        pref=np.zeros((numberofbins,numberofdays,maxage)) # Matrix, [bins, days, maxage, bh vs. reducebh]
        # Set pref[:,0,0,0], which is the "reduced bet hedge version"
        if prefvariance[q]>=0.015:  #Check if variance is so small to just eliminate bet-hedging
            pref[:,0,0]=sci.norm.pdf(x,prefmean[q],prefvariance[q]) # A fly's first day preference gaussian of preference with center around 0
        else: # Make the bin in the middle have all the flies
            pref[math.floor(numberofbins/2),0,0]=flynum

        pref[:,0,0]=pref[:,0,0]/np.sum(pref[:,0,0])*flynum # total # of flies=flynum

        # Now set pref[:,0,0], which is the "reduced bet hedge version"
        if prefvariance[q]*percentbh>=0.015: #Check if variance is so small to just eliminate bet-hedging
            reducedbethedgeinitial=sci.norm.pdf(x,prefmean[q],np.multiply(prefvariance[q],percentbh))
        else:
            reducedbethedgeinitial=np.zeros(numberofbins)
            reducedbethedgeinitial[math.floor(numberofbins/2)]=flynum
        reducedbethedgeinitial[:]=reducedbethedgeinitial[:]/np.sum(reducedbethedgeinitial[:])*flynum # total # of flies=flynum

        # envi=np.zeros((numberofbins,numberofdays))
        # envi[:,0]=sci.norm.pdf(x,envimean,envivariance) # A gaussian of environment with center around 0
        # envi=envi/(np.max(envi))*maxsurvivalrate # Normalizing the maximum envi value and factoring in deathrate
        driftadvantage=np.zeros((numberofdays))
        betadvantage=np.zeros((numberofdays))
        numdays=np.zeros((numberofdays))
        blur=np.zeros((numberofbins,numberofbins)) # [which bin profile it's for, what the distribution is between that bin and all other bins, 2]

        for b in range(numberofbins): # Calculate a value to blur the bins
            blur[b,:]=sci.norm.pdf(x,x[b],driftvariance[q])
            blur[b,:]/=np.sum(blur[b,:])
        if driftmaxdistribution!=0:
            for b in range(numberofbins):
                blur[b,:]=blur[b,:]*sci.norm.pdf(x,0,driftmaxdistribution)
                blur[b,:]/=np.sum(blur[b,:])

        for t in range(1,numberofdays):

            # envi[:,t]=sci.norm.pdf(x,(envimean+gain*np.sin(t*np.pi*2/per)),envivariance) # Making envi a sin wave that changes over time
            # envi[:,t]=envi[:,t]/np.max(envi[:,t])*maxsurvivalrate # Normalizing the envi and multiplying by maxsurvival rate
            # # print('t is: '+str(t))
            # for w in range(2):
            pref[:,t,0]=pref[:,0,0]*birthrate/flynum*np.sum(pref[:,t-1,matureage:]) # Calculate newborn flies
            if adaptivetracking[q]>0:
                pref[:,t,0]=pref[:,t,0]*(1-adaptivetracking[q])+adaptivetracking[q]*birthrate*np.sum(pref[:,t-1,matureage:],1)

                #maybe we should consider putting in some amount of variation on adaptivetracking (shift mean but keep bet hedging?)
            numfliesborntoday=np.sum(pref[:,t,0]) # Calculate the new number of flies born on that day

            betadvantage[t]=np.sum(np.multiply(pref[:,t,0], envi[:,t]))-numfliesborntoday*envi[math.floor(numberofbins/2),t] #

            for a in range(maxage): # For the flies age...

                driftadvantage[t]+=np.sum(np.multiply(pref[:,t-1,a-1], envi[:,t])) # Calculating the number of flies that survive without drift NOTE: Should extend to include BH?
            betadvantage[t]=np.sum(np.multiply(pref[:,t,0], envi[:,t]))-numfliesborntoday*envi[math.floor(numberofbins/2),t]
            
            for a in range(maxage):

                driftadvantage[t]+=np.sum(np.multiply(pref[:,t-1,a-1], envi[:,t])) # Calculating the number of flies that survive without drift #Should extend to include BH

                if a>0:
                    for b in range(numberofbins):
                        if driftvariance[q]<.05/numberofbins: # Check if blur is too low to be worth blurring. NOTE: This number should be based on the limit of sci.norm.pdf
                            pref[b,t,a]+=pref[b,t-1,a-1] # If it is too low, just carry it forward
                        else: # If not, multiply by the blur value and carry it forward
                            pref[:,t,a]+=pref[b,t-1,a-1]*blur[b,:]/np.sum(blur[b,:])
                    pref[:,t,a]=np.multiply(pref[:,t,a], envi[:,t]) # Multiplying the preference to the environment

                if os.path.exists(figuresavepath) and saveloc=='none': 
                    print('not saving, saveloc=none 1')
                if os.path.exists(figuresavepath) and saveloc=='csv':
                    np.savetxt(os.path.join(figuresavepath, 'da'+ str(driftvariance[q])+'ba'+ str(betadvantage[q])+'.csv'),np.concatenate((numdays[:,np.newaxis],driftadvantage[:,np.newaxis],betadvantage[:,np.newaxis]), axis=1), delimiter= ',' , fmt='%i, %.4e, %.4e') #header='Each row is one day\n Day, Drift Advantage, Bet Advantage')
                    print('saving as csv 1')
                if os.path.exists(figuresavepath) and saveloc=='npz': 
                    np.savetxt(os.path.join(figuresavepath, 'da'+ str(driftvariance[q])+'ba'+ str(betadvantage[q])+'.npz'),np.concatenate((numdays[:,np.newaxis],driftadvantage[:,np.newaxis],betadvantage[:,np.newaxis]), axis=1), delimiter= ',' , fmt='%i, %.4e, %.4e', header='Each row is one day\n Day, Drift Advantage, Bet Advantage')
                    print('saving as npz 1')
                if os.path.exists(figuresavepath) and saveloc=='both': 
                    np.savetxt(os.path.join(figuresavepath, 'da'+ str(driftvariance[q])+'ba'+ str(betadvantage[q])+'.npz'),np.concatenate((numdays[:,np.newaxis],driftadvantage[:,np.newaxis],betadvantage[:,np.newaxis]), axis=1), delimiter= ',' , fmt='%i, %.4e, %.4e', header='Each row is one day\n Day, Drift Advantage, Bet Advantage')
                    np.savetxt(os.path.join(figuresavepath, 'da'+ str(driftvariance[q])+'ba'+ str(betadvantage[q])+'.csv'),np.concatenate((numdays[:,np.newaxis],driftadvantage[:,np.newaxis],betadvantage[:,np.newaxis]), axis=1), delimiter= ',' , fmt='%i, %.4e, %.4e', header='Each row is one day\n Day, Drift Advantage, Bet Advantage')
                    print('saving as csv and npz 1')
                # else:
                #     print('not saving 1')

            driftadvantage[t]=np.sum(pref[:,t,:])-driftadvantage[t]-numfliesborntoday
            numdays[t]=int(numdays[t-1]+1)
            # print(np.sum(pref[:,t,0]))
            # print(np.sum(pref[:,t-1,matureage:]))

            # betadvantage[t]=np.sum(pref[:,t,0]-pref[:,t,0])
            # pref[:,t,1:,0]=pref[:,t,:-1,0] #replaced with a-1
        
        sumofpref=np.sum(pref[:,:,:],axis=(0,2))
        sumofprefax2=np.sum(pref[:,:,:],axis=2)
        print(sumofpref.shape)
        print(sumofprefax2.shape)

        if os.path.exists(figuresavepath) and saveloc=='none': 
            print('not saving, saveloc=none 2')
        if os.path.exists(figuresavepath) and saveloc=='csv':
            np.savetxt(os.path.join(figuresavepath, 'da'+ str(driftvariance[q])+'ba'+ str(betadvantage[q])+'envi.csv'),(envi), delimiter= ',') #header='prefmean, prefvariance, driftvariance, birthrate, matureage, percentbh, adaptivetracking, numberofdays')
            np.savetxt(os.path.join(figuresavepath, 'da'+ str(driftvariance[q])+'ba'+ str(betadvantage[q])+'ax2.csv'),(sumofprefax2), delimiter= ',') #header='prefmean, prefvariance, driftvariance, birthrate, matureage, percentbh, adaptivetracking, numberofdays')
            np.savetxt(os.path.join(figuresavepath, 'da'+ str(driftvariance[q])+'ba'+ str(betadvantage[q])+'pref.csv'),(sumofpref), delimiter= ',') #header='prefmean, prefvariance, driftvariance, birthrate, matureage, percentbh, adaptivetracking, numberofdays')
            np.savetxt(os.path.join(figuresavepath, 'da'+ str(driftvariance[q])+'ba'+ str(betadvantage[q])+'parameters.csv'),(np.column_stack([prefmean[0],prefvariance[0],driftvariance[0], birthrate, matureage,percentbh,adaptivetracking[0], numberofdays])), delimiter= ',') #header='prefmean, prefvariance, driftvariance, birthrate, matureage, percentbh, adaptivetracking, sumofpref, numberofdays')
            print('saving as csv 2')
        if os.path.exists(figuresavepath) and saveloc=='npz': 
            np.savetxt(os.path.join(figuresavepath, 'da'+ str(driftvariance[q])+'ba'+ str(betadvantage[q])+'parameters.npz'),(prefmean[0],prefvariance[0],driftvariance[0], birthrate,matureage,percentbh,adaptivetracking[0]), delimiter= ',')
            print('saving as npz 2')
        if os.path.exists(figuresavepath) and saveloc=='both': 
            np.savetxt(os.path.join(figuresavepath, 'da'+ str(driftvariance[q])+'ba'+ str(betadvantage[q])+'parameters.csv'),(prefmean[0],prefvariance[0],driftvariance[0], birthrate,matureage,percentbh,adaptivetracking[0]), delimiter= ',')
            np.savetxt(os.path.join(figuresavepath, 'da'+ str(driftvariance[q])+'ba'+ str(betadvantage[q])+'parameters.npz'),(prefmean[0],prefvariance[0],driftvariance[0], birthrate,matureage,percentbh,adaptivetracking[0]), delimiter= ',')
            print('saving as csv and npz 2')
        # else:
        #     print('not saving 2')


        
        if showgraphs:
            #before = time.perf_counter()
            fig, (ax0, ax1,  ax1d, ax2, ax3, ax4) = plt.subplots(6, 1)
            fig.set_figwidth(10)
            fig.set_figheight(12)
            fig.tight_layout()
            plt.subplots_adjust(hspace=.6)
            c=ax0.pcolormesh(envi)
            fig.colorbar(c,ax=ax0)
            ax0.set_title('Environment (color is fraction of flies of given pref survive)')
            ax0.set_ylabel('Preference')
            ax0.set_xlabel('Day')

            c=ax1.pcolormesh(np.sum(pref[:,:,:],axis=2))
            fig.colorbar(c,ax=ax1)
            ax1.set_title('Fly Preference (color is log(num) flies each day)')
            ax1.set_ylabel('Preference')
            ax1.set_xlabel('Day')

            c=ax1d.pcolormesh(np.sum(pref[:,:,:],axis=2)/np.max(np.sum(pref[:,:,:],axis=2),axis=0))
            fig.colorbar(c,ax=ax1d)
            ax1d.set_title('Fly Preference Distribution (color is %daily flies)')
            ax1d.set_ylabel('Preference')
            ax1d.set_xlabel('Day')

            ax2.plot(np.log(np.sum(pref[:,:,:],axis=(0,2))))
            ax2.set_title('total log(num) flies)') # lowest value is 0.0001 (prefvariance = 0.01 with percent bh=0.01)
            ax2.set_ylabel('log(num) flies)')
            ax2.set_xlabel('Day')
            ax2.set_xlim(0,numberofdays)

            # ax3.plot(driftadvantage)
            ax3.plot(driftadvantage/np.sum(pref[:,:,:],axis=(0,2)))
            ax3.set_title('Change in death rate due to last day\'s drift ')
            ax3.set_ylabel('∆surviving flies/total flies')
            ax3.set_xlabel('Day')
            ax3.set_xlim(0,numberofdays)

            ax4.plot(betadvantage/np.sum(pref[:,:,:],axis=(0,2)))
            ax4.set_title('Change in death rate due to last day\'s bethedging ')
            ax4.set_ylabel('∆surviving flies/total flies')
            ax4.set_xlabel('Day')
            ax4.set_xlim(0,numberofdays)

            fig.colorbar(c,ax=ax2)
            fig.colorbar(c,ax=ax3)
            fig.colorbar(c,ax=ax4)

            fig.suptitle('Bet-hedge variance: '+str(prefvariance[q])+', Drift variance: '+str(driftvariance[q])+', Adaptive Tracking: '+str(adaptivetracking[q]), y=-.05, fontsize=16)

            plt.show()

            if os.path.exists(figuresavepath):
                fig.savefig(os.path.join(figuresavepath,'bh'+str(prefvariance[q])+'dv'+str(driftvariance[q])+'.png'),bbox_inches='tight', pad_inches=.3)
                print(figuresavepath)
            else:
                print('not saving, no valid path')
        #after = time.perf_counter()
        #print(after-before)

        finalpop[q]=np.sum(pref[:,-1,:])
        anymatrix=pref[:,:,1]

        if os.path.exists(figuresavepath):
            np.savetxt(os.path.join(figuresavepath,'prefmatrix'),anymatrix)
            print(figuresavepath)
        else:
            print('not saving, no valid path')


    return finalpop

#Make a csv/npz/both/none (default=none) for the preference
#put all parameters in the npz file
#put the paramaters in header for csv
#make a function that loads a csv filea nd remakes a graph