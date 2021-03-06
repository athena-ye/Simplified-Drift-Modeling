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
import csv

def csvgraph(driftvariance, betadvantage):

    matrix = csv.reader(open('Matrices/driftadvantage' + 'da'+ str(driftvariance)+'ba'+ str(betadvantage)+'.csv', 'r'))

    # figuresavepath='figs'
    # prefvariance=np.linspace(bhlower, bhupper, bhinterval)
    # driftvariance=np.linspace(driftlower, driftupper, driftinterval)

    numruns=bhinterval*driftinterval

    #testing values are right
    # for n in range(numruns):
    #     print('driftvariance'+str(driftvariance[math.floor(n/bhinterval)])+ 'prefvariance'+ str(prefvariance[n%bhinterval]))

    flatmatrix=np.zeros(numruns)

    # for x in range(bhinterval):
    flatmatrix[:]=Parallel(n_jobs=-1, verbose=10)(delayed(matrix for n in range(numruns)))
            # print(["Drift is: ", driftvariance[y], "Bet-hedging is: ", driftvariance[x]] )
            # print(matrix)
    print(flatmatrix)

    # flatmatrix=np.matrix(flatmatrix)
    # flatmatrix=
    matrix=np.zeros((bhinterval,driftinterval))
    matrix[:,:]=flatmatrix.reshape((bhinterval, driftinterval), order='F')
    matrixlog=np.log(matrix)
    print(matrixlog)
    print(matrix)



    bhmargin=(bhupper-bhlower)/(bhinterval-1)/2
    driftmargin=(driftupper-driftlower)/(driftinterval-1)/2
    prefvariancemesh=np.linspace(bhlower-bhmargin, bhupper+bhmargin, bhinterval+1)
    driftvariancemesh=np.linspace(driftlower-driftmargin, driftupper+driftmargin, driftinterval+1)
    print(prefvariancemesh)
    print(driftvariancemesh)

# driftvariancegrid2, prefvariancegrid2= np.meshgrid(driftvariance, prefvariance)

    # print(matrix)

    fig,ax=plt.subplots()
    scale=max(-np.min(matrixlog),np.max(matrixlog))
    c=ax.pcolormesh(driftvariancemesh, prefvariancemesh, matrixlog, shading='flat', cmap='RdBu',  vmin=-scale, vmax=scale)


    # c=ax.pcolormesh(driftvariancemesh, prefvariancemesh, matrixlog, shading='flat', cmap='RdBu',  vmin=-20, vmax=20)
    ax.set_xlabel('Drift')
    ax.set_ylabel('Bet-Hedging')
    fig.colorbar(c, ax=ax)
    ax.set_title('Log of Final Population')
    plt.show

    if os.path.exists(figuresavepath):
        fig.savefig(os.path.join(figuresavepath,'heatmap.png'),bbox_inches='tight', pad_inches=.3)
        print(figuresavepath)
    else:
        print('not saving, no valid path')

    return matrix, driftvariance, prefvariance