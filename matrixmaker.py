# make the matrix
# call simpledrift to fill in the values for different drift and bh
# pass on to heatmap
# pay attention to axes so they give the correct values
import numpy as np
import matplotlib.pyplot as plt
import simpledrift as sd
import math
from joblib import Parallel, delayed
import os

def matrixmaker(envi, bhlower, bhupper, bhinterval, driftlower, driftupper, driftinterval, runindex=0, saveloc='none'):

    prefmean=0
    birthrate=10
    matureage=10
    percentbh=0.01
    adaptivetracking=0
    showgraphs=False

    figuresavepath='figs'
    prefvariance=np.linspace(bhlower, bhupper, bhinterval)
    driftvariance=np.linspace(driftlower, driftupper, driftinterval)

    numruns=bhinterval*driftinterval

    #testing values are right
    # for n in range(numruns):
    #     print('driftvariance'+str(driftvariance[math.floor(n/bhinterval)])+ 'prefvariance'+ str(prefvariance[n%bhinterval]))

    flatmatrix=np.zeros(numruns)

    flatmatrix[:]=Parallel(n_jobs=-1, verbose=10)(delayed(sd.driftmodeling)(envi, prefmean, [prefvariance[n%bhinterval]], [driftvariance[math.floor(n/bhinterval)]], adaptivetracking, birthrate,matureage, percentbh, showgraphs, figuresavepath) for n in range(numruns))
            # print(["Drift is: ", driftvariance[y], "Bet-hedging is: ", driftvariance[x]] )

    matrix=np.zeros((bhinterval,driftinterval))
    matrix[:,:]=flatmatrix.reshape((bhinterval, driftinterval), order='F')
    matrixlog=np.log(matrix)

    bhmargin=(bhupper-bhlower)/(bhinterval-1)/2
    driftmargin=(driftupper-driftlower)/(driftinterval-1)/2
    prefvariancemesh=np.linspace(bhlower-bhmargin, bhupper+bhmargin, bhinterval+1)
    driftvariancemesh=np.linspace(driftlower-driftmargin, driftupper+driftmargin, driftinterval+1)
    # driftvariancegrid2, prefvariancegrid2= np.meshgrid(driftvariance, prefvariance)

    fig,ax=plt.subplots()
    scale=max(-np.min(matrixlog),np.max(matrixlog))
    c=ax.pcolormesh(driftvariancemesh, prefvariancemesh, matrixlog, shading='flat', cmap='RdBu',  vmin=-scale, vmax=scale)


    # c=ax.pcolormesh(driftvariancemesh, prefvariancemesh, matrixlog, shading='flat', cmap='RdBu',  vmin=-20, vmax=20)
    ax.set_xlabel('Drift')
    ax.set_ylabel('Bet-Hedging')
    fig.colorbar(c, ax=ax)
    ax.set_title('Log of Final Population')
    plt.show

    if os.path.exists(figuresavepath) and saveloc=='none':
        print('none=not saving, no valid path')
    if os.path.exists(figuresavepath) and saveloc=='npz':
        heatmapname='R'+str(runindex)+'_heatmap.png'
        fig.savefig(os.path.join(figuresavepath,heatmapname),bbox_inches='tight', pad_inches=.3)
        filename='R'+str(runindex)+'_Env_FinalPopulations.npz'
        np.savez(os.path.join(figuresavepath,filename), finalpopulations=matrix, prefvariancemesh=prefvariancemesh, driftvariancemesh=driftvariancemesh, envi=envi)
        print(figuresavepath)
    if os.path.exists(figuresavepath) and saveloc=='csv':
        heatmapname='R'+str(runindex)+'_heatmap.png'
        fig.savefig(os.path.join(figuresavepath,heatmapname),bbox_inches='tight', pad_inches=.3)
        filename='R'+str(runindex)+'_Env_FinalPopulations.csv'
        np.savez(os.path.join(figuresavepath,filename), finalpopulations=matrix, prefvariancemesh=prefvariancemesh, driftvariancemesh=driftvariancemesh, envi=envi)
        print(figuresavepath)
    if os.path.exists(figuresavepath) and saveloc=='both':
        heatmapname='R'+str(runindex)+'_heatmap.png'
        fig.savefig(os.path.join(figuresavepath,heatmapname),bbox_inches='tight', pad_inches=.3)
        filename='R'+str(runindex)+'_Env_FinalPopulations.npz'
        np.savez(os.path.join(figuresavepath,filename), finalpopulations=matrix, prefvariancemesh=prefvariancemesh, driftvariancemesh=driftvariancemesh, envi=envi)
        print(figuresavepath)
        heatmapname='R'+str(runindex)+'_heatmap.png'
        fig.savefig(os.path.join(figuresavepath,heatmapname),bbox_inches='tight', pad_inches=.3)
        filename='R'+str(runindex)+'_Env_FinalPopulations.csv'
        np.savez(os.path.join(figuresavepath,filename), finalpopulations=matrix, prefvariancemesh=prefvariancemesh, driftvariancemesh=driftvariancemesh, envi=envi)
        print(figuresavepath)
    else:
        print('error, no valid path for specified saveloc')

    return matrix, driftvariance, prefvariance
