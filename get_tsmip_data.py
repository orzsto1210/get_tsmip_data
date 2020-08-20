from obspy import UTCDateTime, Stream, Trace
from obspy.signal.trigger import ar_pick
import numpy as np
import matplotlib.pyplot as plt
import joblib
import os

def unpack(inflie):
    with open(inflie) as f:
        station = f.readline().split()[1]
        instrument = f.readline().split()[1]
        starttime = f.readline().split()[1]
        starttime = UTCDateTime(int(starttime[0:4]), int(starttime[5:7]), int(starttime[8:10]), int(starttime[11:13]), int(starttime[14:16]), float(starttime[17:]))
        recordLength = float(f.readline().split()[1])
        sampleRate = float(f.readline().split()[1])
        npts = int(recordLength * sampleRate)
        # print(station, instrument, starttime, recordLength, sampleRate, npts)

        f.readline() #AmplitudeUnit
        pga_z_info = f.readline().split()
        pga_z, pga_z2 = float(pga_z_info[2][:-1]), float(pga_z_info[3])

        pga_n_info = f.readline().split()
        pga_n, pga_n2 = float(pga_n_info[2][:-1]), float(pga_n_info[3])

        pga_e_info = f.readline().split()
        pga_e, pga_e2 = float(pga_e_info[2][:-1]), float(pga_e_info[3])
        
        # print(pga_z, pga_z2, pga_n, pga_n2, pga_e, pga_e2)
        pga = max(np.abs([pga_z, pga_z2, pga_n, pga_n2, pga_e, pga_e2]))
        # print(pga)

        f.readline() #DataSequence
        f.readline() #Data

        z_data = []
        n_data = []
        e_data = []

        for line in f.readlines():
            line = line.split()
            z = float(line[1])
            n = float(line[2])
            e = float(line[3])
            z_data.append(z)
            n_data.append(n)
            e_data.append(e)

    zne_data = np.array([z_data, n_data, e_data])
        
    st = Stream()
    for ch in range(3):
        channel = 'Ch' + str(ch+1)    

        stats = {'station': station, 'npts': npts, 'sampling_rate': sampleRate, 
                'starttime': starttime, 'channel': channel}
        data = zne_data[ch]
        sttmp = Stream([Trace(data=data, header=stats)])
        st += sttmp

    st.resample(100.0)

    return st

def calc_pga_intensity(my_st):
    pga_z = max(abs(my_st[0].data))
    pga_n = max(abs(my_st[1].data))
    pga_e = max(abs(my_st[2].data))
    pga = max(pga_z, pga_n, pga_e)
    pga = round(pga, 1)
    # print(pga)
    if pga >= 0.0 and pga < 0.8:
        intensity = 0
    elif pga >= 0.8 and pga < 2.5:
        intensity = 1
    elif pga >= 2.5 and pga < 8.0:
        intensity = 2
    elif pga >= 8 and pga < 25:
        intensity = 3
    elif pga >= 25 and pga < 80:
        intensity = 4
    elif pga >= 80 and pga < 140:
        intensity = 5.1
    elif pga >= 140 and pga < 250:
        intensity = 5.5
    elif pga >= 250 and pga < 440:
        intensity = 6.1
    elif pga >= 440 and pga < 800:
        intensity = 6.5
    elif pga >= 800:
        intensity = 7
    else:
        intensity = -1

    return pga, intensity

def calc_pgv_intensity(my_st):
    st = my_st.copy()

    pgv_list = []
    for i in range(len(st[0].data)):
        tmp_pgv = ((st[0].data[i])**2+(st[1].data[i])**2+(st[2].data[i])**2)**(0.5)
        pgv_list.append(tmp_pgv)
    pgv = max(pgv_list)
    pgv = round(pgv, 2)
    
    if pgv < 0.2:
        intensity = 0
    elif pgv >= 0.2 and pgv < 0.7:
        intensity = 1
    elif pgv >= 0.7 and pgv < 1.9:
        intensity = 2
    elif pgv >= 1.9 and pgv < 5.7:
        intensity = 3
    elif pgv >= 5.7 and pgv < 15:
        intensity = 4
    elif pgv >= 15 and pgv < 30:
        intensity = 5.1
    elif pgv >= 30 and pgv < 50:
        intensity = 5.5
    elif pgv >= 50 and pgv < 80:
        intensity = 6.1
    elif pgv >= 80 and pgv < 140:
        intensity = 6.5
    elif pgv >= 140:
        intensity = 7
    else:
        intensity = -1
        
    return pgv, intensity

def plot_vel_waveform(my_st, p_arrival, outname):
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(3, 1, 1)

    ax.plot(my_st[0].data, "k-")
    ymin, ymax = ax.get_ylim()
    ax.vlines(p_arrival, ymin, ymax, color='r', linewidth=2, label="p_arrival")
    plt.legend()
    ax.set_ylabel("cm/sec")

    ax = fig.add_subplot(3, 1, 2)
    ax.plot(my_st[1].data, "k-")
    ymin, ymax = ax.get_ylim()
    ax.vlines(p_arrival, ymin, ymax, color='r', linewidth=2, label="p_arrival")
    plt.legend()
    ax.set_ylabel("cm/sec")

    ax = fig.add_subplot(3, 1, 3)
    ax.plot(my_st[2].data, "k-")
    ymin, ymax = ax.get_ylim()
    ax.vlines(p_arrival, ymin, ymax, color='r', linewidth=2, label="p_arrival")
    plt.legend()
    ax.set_ylabel("cm/sec")

    plt.suptitle(outname,size=14)

    if not os.path.exists("vel_png"):
        os.makedirs("vel_png")
    outname = "vel_png/" + outname+".png"
    
    plt.savefig(outname)
    # plt.show()

def plot_acc_waveform(my_st, p_arrival, outname):
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(3, 1, 1)

    ax.plot(my_st[0].data, "k-")
    ymin, ymax = ax.get_ylim()
    ax.vlines(p_arrival, ymin, ymax, color='r', linewidth=2, label="p_arrival")
    plt.legend()
    ax.set_ylabel("gal")

    ax = fig.add_subplot(3, 1, 2)
    ax.plot(my_st[1].data, "k-")
    ymin, ymax = ax.get_ylim()
    ax.vlines(p_arrival, ymin, ymax, color='r', linewidth=2, label="p_arrival")
    plt.legend()
    ax.set_ylabel("gal")

    ax = fig.add_subplot(3, 1, 3)
    ax.plot(my_st[2].data, "k-")
    ymin, ymax = ax.get_ylim()
    ax.vlines(p_arrival, ymin, ymax, color='r', linewidth=2, label="p_arrival")
    plt.legend()
    ax.set_ylabel("gal")

    plt.suptitle(outname,size=14)

    if not os.path.exists("acc_png"):
        os.makedirs("acc_png")
    outname = "acc_png/" + outname+".png"

    plt.savefig(outname)
    # plt.show()

def get_p_arrival_data(my_st, p_arrival, sample, outname):
    E = my_st[2].data[p_arrival:p_arrival+sample]
    N = my_st[1].data[p_arrival:p_arrival+sample]
    Z = my_st[0].data[p_arrival:p_arrival+sample]

    return E, N, Z

def calculateSNR(dd_before, dd_after):
    s = np.sum(np.square(dd_after))
    n = np.sum(np.square(dd_before))
    snr = s / n
    return snr

def get_snr(my_st, p_arrival):
    E = my_st[2].data
    N = my_st[1].data
    Z = my_st[0].data

    start = p_arrival - 200
    end = p_arrival + 200

    snr_e = round(calculateSNR(E[start:p_arrival], E[p_arrival:end+1]), 3)
    snr_n = round(calculateSNR(N[start:p_arrival], N[p_arrival:end+1]), 3)
    snr_z = round(calculateSNR(Z[start:p_arrival], Z[p_arrival:end+1]), 3)
    snr = round(np.mean([snr_e, snr_n, snr_z]), 3)
    return snr

def save_acceleration_data(inflie):
    st = unpack(inflie)
    my_st = st.copy()

    p_pick, s_pick = ar_pick(my_st[0].data, my_st[1].data, my_st[2].data, my_st[0].stats.sampling_rate,
                         1.0, 20.0, 1.0, 0.1, 4.0, 1.0, 2, 8, 0.1, 0.2)
    p_arrival = int(p_pick*100)
    print(p_arrival)

    snr = get_snr(my_st, p_arrival)
    pga, intensity = calc_pga_intensity(my_st)

    outname = "_".join([inflie[:-4], my_st[0].stats.station, str(intensity), str(pga), str(snr)]) 
    plot_acc_waveform(my_st, p_arrival, outname)
    sample = 100
    E, N, Z = get_p_arrival_data(my_st, p_arrival, sample, outname)

    if not os.path.exists("acc_joblib"):
        os.makedirs("acc_joblib")
    outname = "acc_joblib/" + outname+".joblib"
    joblib.dump([snr, pga, intensity, E, N, Z], outname)

def save_velocity_data(inflie):
    st = unpack(inflie)
    my_st = st.copy()
    my_st.integrate()
    my_st.filter("highpass", freq=0.075)

    p_pick, s_pick = ar_pick(my_st[0].data, my_st[1].data, my_st[2].data, my_st[0].stats.sampling_rate,
                         1.0, 20.0, 1.0, 0.1, 4.0, 1.0, 2, 8, 0.1, 0.2)
    p_arrival = int(p_pick*100)
    print(p_arrival)

    snr = get_snr(my_st, p_arrival)
    pgv, intensity = calc_pgv_intensity(my_st)

    outname = "_".join([inflie[:-4], my_st[0].stats.station, str(intensity), str(pgv), str(snr)]) 
    plot_vel_waveform(my_st, p_arrival, outname)
    sample = 100
    E, N, Z = get_p_arrival_data(my_st, p_arrival, sample, outname)

    if not os.path.exists("vel_joblib"):
        os.makedirs("vel_joblib")
    outname = "vel_joblib/" + outname+".joblib"
    joblib.dump([snr, pgv, intensity, E, N, Z], outname)

inflie = "10026302.SSX.txt"
save_acceleration_data(inflie)
save_velocity_data(inflie)
