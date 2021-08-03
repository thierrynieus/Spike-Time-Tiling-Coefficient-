import numpy as np

def run_T(spiketrain, N, dt, t_start, t_stop):
    """
    Calculate the proportion of the total recording time 'tiled' by spikes.
    """
    time_A = 2 * N * dt  # maxium possible time

    if N == 1:  # for just one spike in train
        if spiketrain[0] - t_start < dt:
            time_A = time_A - dt + spiketrain[0] - t_start
        elif spiketrain[0] + dt > t_stop:
            time_A = time_A - dt - spiketrain[0] + t_stop

    else:  # if more than one spike in train
        '''
            This part of code speeds up calculation with respect to the original version
        '''
        diff = np.diff(spiketrain)
        idx = np.where(diff<(2*dt))[0]
        Lidx = len(idx)
        time_A = time_A - 2 * Lidx * dt + diff[idx].sum()

        if (spiketrain[0] - t_start) < dt:
            time_A = time_A + spiketrain[0] - dt - t_start

        if (t_stop - spiketrain[N - 1]) < dt:
            time_A = time_A - spiketrain[-1] - dt + t_stop

    T = (time_A / (t_stop - t_start)) #.item()
    return T

def run_P(spiketrain_1, spiketrain_2, N1, N2, dt):
    """
    Check every spike in train 1 to see if there's a spike in train 2
    within dt
    """
    Nab = 0
    j = 0
    for i in range(N1):
        L=0
        while j < N2:  # don't need to search all j each iteration
            if np.abs(spiketrain_1[i] - spiketrain_2[j]) <= dt:
                Nab = Nab + 1
                L+=1
                break
            elif spiketrain_2[j] > spiketrain_1[i]:
                break
            else:
                j = j + 1
    return Nab

def sttc_fast(spiketrain_1,spiketrain_2,N1,N2,dt,t_start=0,t_stop=3e5):
    '''
    '''
    TA = run_T(spiketrain_1, N1, dt, t_start=t_start, t_stop=t_stop)
    #print(TA)
    TB = run_T(spiketrain_2, N2, dt, t_start=t_start, t_stop=t_stop)
    #print(TB)
    PA = run_P(spiketrain_1, spiketrain_2, N1, N2, dt)
    PA = PA / float(N1)
    PB = run_P(spiketrain_2, spiketrain_1, N2, N1, dt)
    PB = PB / float(N2)
    index = 0.5 * (PA - TB) / (1 - PA * TB) + 0.5 * (PB - TA) / (1 - PB * TA)
    return index

 
