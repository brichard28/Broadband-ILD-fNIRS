import os, sys
import getopt
#import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.signal
from pysndfile.sndio import read as sfread
from pysndfile.sndio import write as sfwrite
import psylab
import ir
#import spyral
#import ild

__doc__ = """
Generates stimuli for nirs/eeg study looking at ild-based informational masking release.

Word array should have 18 words, 2-4 of which should be colors. 

Each word is 300 ms long, and there should be 300 ms silence between, for 10.8 sec total.

4 stim per condition

Color words were added to the masker stream, as in zhang et al.

"""

def main(subjectID):

    exp = 'srm-nirs-eeg-1'

    n = 4              # # of runs, Zhang used 6. Here we use 4 since side is an IV (8 total, 4 per side)
    nwords = 18        # # of words in a list, Zhang used 25 (15 sec) (18 * (.3word + .3 silence) = 10.8s)
    ncolors = [2, 4]   # list of the range of possible numbers of colors, Zhang used 3-5
    silence = .3       # Amount of silence between words
    masker_delay = .15 # Delay between target start and masker start. Has been effectively zero <= nirs-im-8
                       # Won't change overall length. IE, will eat into silence var above; so if word dur == .3
                       # silence == .3, and masker_delay == .15, there will be .15 s of silence between offset
                       # of a masker and onset of next target
    targ_talker = 'Bob' # Target Talker (Bob or Mike)
    mask_talker = 'Bob' # Masker Talker

    target_first = -1  # if masker_delay > 0; 1 = target first; 0 = masker first; -1 = random

    # location cue; a tone to indicate which location to listen to.
    #loc_cue_f = 500
    #loc_cue_dur = 500
    loc_cue_isi = 1700 # ms; use a number for fixed, a list of size 2 for random
    loc_cue_word = 'bag'
    
    # Soundfile paths
    sf_path_out = os.path.expanduser(f'~/work/Python/stim/bobmike/{exp}/s_{subjectID}')
    data_path_out = os.path.expanduser(f'~/work/Python/data')
    sf_path_in = '~/work/stim/bobmike'
    sf_fs = 44100

    # Conditions
    # IMPORTANT: ind=0 should always be the practice condition
    sides =   ['l',    'l',      'r',      'l',     'r',     'l',      'r',      'l',     'r',     'l',      'r',      'l',      'r',     'l',      'r',      'l',      'r',     'l',      'r',      'l',     'r',    ]
    maskers = ['none', 'speech', 'speech', 'noise', 'noise', 'speech', 'speech', 'noise', 'noise', 'speech', 'speech', 'noise',  'noise', 'speech', 'speech', 'noise',  'noise', 'speech', 'speech', 'noise', 'noise', ]
    itds =    ['0',    '50',     '50',     '50',    '50',    '500',    '500',    '500',   '500',   '0',      '0',      '0',      '0',     '0',      '0',      '0',      '0',     '500',    '500',    '500',   '500',   ]
    ilds =    ['0',    '0',      '0',      '0',     '0',     '0',      '0',      '0',     '0',     '10',     '10',     '10',     '10',    '70n',    '70n',    '70n',    '70n',    '0',      '0',      '0',     '0',]
    controls =['0',    '0',      '0',      '0',     '0',     '0',      '0',      '0',     '0',     '0',      '0',      '0',      '0',     '0',      '0',      '0',      '0',     '1',      '1',      '1',     '1',]

    # Token info
    colors = ['blue', 'green', 'red',    'white']
    objects = ['bag', 'card',  'chairs', 'desks', 'glove', 'hat', 'pen', 'shoe', 'sock', 'spoons', 'tables', 'toy']

    snr = 0
    m_az = 30

    # ILD stuff not used, whole-waveform ilds applied instead
    ild_bands= psylab.signal.logspace(62.5,8000,22) # psylab.signal.frequency.logspace(70, 2240, 5+1) # 5 1-oct bands
    ild_max_itd_us = 1000.
    ild_max_ild_db = 15
    ild_wsize = 882
    ild_stepsize = 441
    ild_attack = 110
    fun1 = np.zeros(50)
    fun2 = np.linspace(0, 1, 250)
    fun3 = np.ones(np.int32(ild_max_itd_us) - 300 + 1)
    fun = np.concatenate((fun1, fun2, fun3))
    #fun = np.ones(np.int32(ild_max_itd_us)+1)
    nofun = np.zeros(np.int32(ild_max_itd_us))
    # Here, no positive gain, only attenuation
    ild_table_l = np.concatenate((nofun, -fun))
    ild_table_r = np.concatenate((-fun[::-1], nofun))


    # Vocoder parameters From Zhang, Alamatsaz & Ihlefeld 2021
    # Use filtfilt for 72 dB/oct and 0 phase shift
    v_channel_fcs = psylab.signal.logspace(300,10000,16)
    v_channel_cos = np.zeros((2, 16))

    # format: (b,a)/(:)/(hi,lo)/(channel)
    v_channel_coefs = np.zeros((2,7,2,16))
    for i in np.arange(v_channel_fcs.size):
        f1,f2 = psylab.signal.oct2f(v_channel_fcs[i], 1/20.)
        v_channel_cos[:,i] = f1,f2
        v_channel_coefs[:,:,0,i] = scipy.signal.butter(6, f1/sf_fs/2., btype='highpass')
        v_channel_coefs[:,:,1,i] = scipy.signal.butter(6, f2/sf_fs/2., btype='lowpass')

    hrtf_path = 'stim/hrtf_kemar_0el.npy'
    dat = np.load(hrtf_path)

    if not os.path.exists(sf_path_out):
        os.makedirs(sf_path_out)
    if not os.path.exists(data_path_out):
        os.makedirs(data_path_out)
#    fid_color_times_t = open(os.path.join(data_path_out, f'target_color_times__s_{subjectID}.txt'), 'w')
#    fid_color_times_m = open(os.path.join(data_path_out, f'masker_color_times__s_{subjectID}.txt'), 'w')
#    fid_target_words = open(os.path.join(data_path_out, f'target_words__s_{subjectID}.txt'), 'w')
    fid_all_times = open(os.path.join(data_path_out, f'{exp}__s_{subjectID}__Word_Times.csv'), 'w')

    #cue_data = psylab.signal.tone(loc_cue_f, sf_fs, loc_cue_dur, amp=.3)
    #cue_data = psylab.signal.ramps(cue_data,sf_fs)
    cue_data = os.path.join(sf_path_in, f'{targ_talker.lower()}_object', f"{loc_cue_word}_short.wav")
    rng = np.random.default_rng()

    header = "Condition,Run,Type"
    for k in np.arange(nwords):
        header += f",Word{k},Time{k}"
    header += "\n"
    fid_all_times.write(header)

    #for masker,itd,ild in zip(maskers, itds, ilds):
    pract = True
    for masker,side,ild,itd,control in zip(maskers,sides,ilds,itds,controls):

        if isinstance(loc_cue_isi,list):
            # list means random in specified range
            cue_isi_dur = rng.uniform(loc_cue_isi[0], loc_cue_isi[1], size=1)[0]
        else:
            cue_isi_dur = loc_cue_isi
        cue_isi_samp = psylab.signal.ms2samp(cue_isi_dur, sf_fs)
        cue_isi_data = np.zeros(cue_isi_samp)

        if itd != '0':
            itd_data = np.zeros(int(int(itd) / 1000000 * sf_fs))
        else:
            itd_data = np.array(())
        lead_silence_data_pre = np.zeros(0)
        lead_silence_data_post = np.zeros(np.int32(sf_fs*silence)-len(itd_data))
        lag_silence_data_pre = np.zeros(np.int32(sf_fs*masker_delay))
        lag_silence_data_post = np.zeros(np.int32(sf_fs*(silence-masker_delay)))

        sf_inds = np.random.randint(1, 100, size=n)
        if pract:
            this_sf_path_out_sub = 'practice'
            pract = False
        else:
            this_sf_path_out_sub = f"m_{masker}__ild_{ild}__itd_{itd}__targ_{side}__control_{control}" # For textfiles
        this_sf_path_out = os.path.join(sf_path_out, this_sf_path_out_sub)
        if not os.path.exists(this_sf_path_out):
            os.makedirs(this_sf_path_out)
        print(f"{this_sf_path_out_sub}: {itd_data}")
        # Generate lists of inds, ensure: 
        #   - each ind is unique across lists
        #   - no targ ind is same as previous
        for k in np.arange(n):
            # used for writing words and times to fid_all_times. In S
            curr_time = (loc_cue_dur + cue_isi_dur) / 1000

            targ_str = f"{this_sf_path_out_sub},{k},Target,"
            mask_str = f"{this_sf_path_out_sub},{k},Masker,"
            nobjectlists = nwords // len(objects)
            robjectlists = nwords % len(objects)
            objects_25 = []
            for j in range(nobjectlists):
                objects_25.extend(objects)
            for j in range(robjectlists):
                objects_25.append(objects[np.random.choice(len(objects))])
            #objects_25 = objects + objects + [objects[np.random.choice(len(objects))]]
            prev_targ = None
            targ_arr =  []
            mask1_arr = []
            mask2_arr = []
            sf_targ = []
            sf_mask1 = []
            sf_mask2 = []
            for j in np.arange(nwords):
                this_i = np.random.choice(nwords, 3, replace=False)
                if prev_targ:
                    while this_i[0] == prev_targ:
                        this_i = np.random.choice(nwords, 3, replace=False)
                prev_targ = this_i[0]
                targ_arr.append(this_i[0])
                mask1_arr.append(this_i[1])
                mask2_arr.append(this_i[2])

            targ_words =  []
            mask1_words = []
            mask2_words = []

            targ_times = []
            mask2_times = []

            targ_data =  np.zeros(sf_fs)
            mask1_data = np.zeros(sf_fs)
            mask2_data = np.zeros(sf_fs)

            colors_n_t = np.random.randint(ncolors[0],ncolors[1]+1) # between 3 - 5 colors
            colors_n_m = np.random.randint(ncolors[0],ncolors[1]+1)
            colors_ind = np.random.choice(nwords, colors_n_t+colors_n_m, replace=False)
            colors_i_t = colors_ind[:colors_n_t]
            colors_i_m = colors_ind[colors_n_t:]
            color_times_t = []
            color_times_m = []

            for i in np.arange(len(targ_arr)):

                # TODO: this_targ_i and targ_arr do not seem to be used, as 
                # we are choosing from objects or colors at random below
                this_targ_i =  targ_arr[i]
                this_mask1_i = mask1_arr[i]
                this_mask2_i = mask2_arr[i]
                if i in colors_i_t:
                    this_targ_word = colors[np.random.randint(len(colors))]
                    this_targ_kind = f'{targ_talker.lower()}_color'
                    color_times_t.append(str(np.round(i*(.3+silence), 3))) # 600 ms between words
                else:
                    #this_targ_word = objects[np.random.randint(len(objects))]
                    this_targ_word = objects_25[this_targ_i]
                    this_targ_kind = f'{targ_talker.lower()}_object'
                targ_words.append(this_targ_word)

                if masker == 'speech':
                    if i in colors_i_m:
                        this_mask2_word = colors[np.random.randint(len(colors))]
                        this_mask2_kind = f'{mask_talker.lower()}_color'
                        color_times_m.append(str(np.round(i*(.3+silence+masker_delay), 3))) # 600 ms between words
                    else:
                        this_mask1_word = objects_25[this_mask1_i]
                        this_mask1_kind = f'{mask_talker.lower()}_object'
                        this_mask2_word = objects_25[this_mask2_i]
                        this_mask2_kind = f'{mask_talker.lower()}_object'
                elif masker == 'noise':
                    this_mask1_word = "noise"
                    this_mask1_kind = ""
                    this_mask2_word = "noise"
                    this_mask2_kind = ""
                else:
                    this_mask1_word = "none"
                    this_mask1_kind = ""
                    this_mask2_word = "none"
                    this_mask2_kind = ""

                mask1_words.append(this_mask1_word)
                mask2_words.append(this_mask2_word)
                this_targ_sf = os.path.join(sf_path_in,  this_targ_kind,   f"{this_targ_word}_short.wav")
                sf_targ.append(this_targ_sf)
                this_targ_data,fs,_enc = sfread(this_targ_sf)
            #    target_first = -1  # if masker_delay > 0; 1 = target first; 0 = masker first; -1 = random
                if target_first in [0, 1]:
                    t_f = target_first
                else:
                    t_f = np.random.choice([0,1])
                if t_f == 1:
                    targ_data = np.concatenate((targ_data, lead_silence_data_pre, this_targ_data, lead_silence_data_post))
                    this_targ_time = np.round(curr_time, 3)
                else:
                    targ_data = np.concatenate((targ_data, lag_silence_data_pre, this_targ_data, lag_silence_data_post))
                    this_targ_time = np.round(curr_time + masker_delay)

                if this_mask1_word == 'noise':
                    this_mask1_data = psylab.signal.white(this_targ_data.shape[0])
                    this_mask2_data = psylab.signal.white(this_targ_data.shape[0])
                    this_mask1_data *= psylab.signal.rms(this_targ_data) / psylab.signal.rms(this_mask1_data)
                    this_mask2_data *= psylab.signal.rms(this_targ_data) / psylab.signal.rms(this_mask2_data)
                elif this_mask1_word == 'none':
                    this_mask1_data = np.zeros(this_targ_data.shape[0])
                    this_mask2_data = np.zeros(this_targ_data.shape[0])
                else:
                    #print(this_mask1_word)
                    this_mask1_sf = os.path.join(sf_path_in, this_mask1_kind, f"{this_mask1_word}_short.wav")
                    sf_mask1.append(this_mask1_sf)
                    this_mask1_data,fs,_enc = sfread(this_mask1_sf)

                    this_mask2_sf = os.path.join(sf_path_in, this_mask2_kind, f"{this_mask2_word}_short.wav")
                    sf_mask2.append(this_mask2_sf)
                    this_mask2_data,fs,_enc = sfread(this_mask2_sf)

                if t_f == 1:
                    mask1_data = np.concatenate((mask1_data, lag_silence_data_pre, this_mask1_data, lag_silence_data_post))
                    mask2_data = np.concatenate((mask2_data, lag_silence_data_pre, this_mask2_data, lag_silence_data_post))
                    this_mask_time = np.round(curr_time + masker_delay, 3)
                else:
                    mask1_data = np.concatenate((mask1_data, lead_silence_data_pre, this_mask1_data, lead_silence_data_post))
                    mask2_data = np.concatenate((mask2_data, lead_silence_data_pre, this_mask2_data, lead_silence_data_post))
                    this_mask_time = np.round(curr_time, 3)

                targ_times.append(this_targ_time)
                mask2_times.append(this_mask_time)
                #targ_times += f",{this_targ_word},{this_targ_time}"
                #mask_times += f",{this_mask2_word},{this_mask_time}"

                curr_time += .6 # TODO: Really shouldn't hardcode this val

            ## Do nothing to target, which is always at midline
            targ_data_l = targ_data.copy()
            targ_data_r = targ_data.copy()
            cue_data_l = cue_data.copy()
            cue_data_r = cue_data.copy()
            #if ild == 'inf':
            #    targ_data_r = np.zeros(targ_data_r.shape)
            #elif ild != '0':
            #    targ_data_r = psylab.signal.atten(targ_data_r, int(ild))

    #ilds =    ['0',    '0',      '10',     '20',     '30',     'inf',    '70',      '70'] # 70 here is az not ild
    #ftypes =  ['bb',   'bb',     'bb',     'bb',     'bb',     'bb',     'bb',      'nat'] # type is a py keyword
    #ild_70_bb = 16.66 # Use this ILD in dB when ild == 70 and ftype == bb
            if masker == 'none':
                mask1_data_l = np.zeros(targ_data.shape)
                mask1_data_r = np.zeros(targ_data.shape)
                mask2_data_l = np.zeros(targ_data.shape)
                mask2_data_r = np.zeros(targ_data.shape)
                targ_data_l = psylab.signal.atten(targ_data_l, -6) # W no maskers, the overall level is low
                targ_data_r = psylab.signal.atten(targ_data_r, -6) # Probably because there are only 7 bands total
            else:
                mask1_data_l = psylab.signal.atten(mask1_data.copy(),snr) #+3) # +3 to account for 2 maskers
                mask1_data_r = psylab.signal.atten(mask1_data.copy(),snr) #+3)
                mask2_data_l = psylab.signal.atten(mask2_data.copy(),snr) #+3)
                mask2_data_r = psylab.signal.atten(mask2_data.copy(),snr) #+3)
#    maskers = ['none', 'speech', 'speech', 'speech', 'speech', 'speech', 'speech', 'noise']
#    itds =    ['500',  '0',      '0',      '0',      '0',      '0',      '500',    '500']
#    ilds =    ['0',    'inf',    '0',      '10',     '20',     '30',     '0',      '0']

                # 0 is handled correctly
                # We only use 1 masker for nirs-im-8, so ignore mask1 and apply spatial cues to targ
                targ_data_l = np.concatenate((targ_data_l, itd_data))
                targ_data_r = np.concatenate((itd_data, targ_data_r))
                mask1_data_l = np.concatenate((itd_data, mask1_data_l))
                mask1_data_r = np.concatenate((mask1_data_r, itd_data))
                mask2_data_l = np.concatenate((mask2_data_l, itd_data))
                mask2_data_r = np.concatenate((itd_data, mask2_data_r))
                cue_data_l = np.concatenate((cue_data_l, itd_data))
                cue_data_r = np.concatenate((itd_data, cue_data_r))

                if ild == 'inf':
                    targ_data_r = np.zeros(len(targ_data_r))
                    mask1_data_r = np.zeros(len(mask1_data_r))
                    mask2_data_l = np.zeros(len(mask2_data_l))
                    cue_data_r = np.zeros(len(cue_data_r))
                    
                    #mask1_data_r = psylab.signal.atten(mask1_data_r, 300)
                    #mask1_data_l = psylab.signal.atten(mask1_data_l, 300)
                elif ild.endswith("n"):
                    # Natural ILDs; apply the magnitude spectrum of the specified ir
                    # MITs KEMAR HRTF set is mono and thus symmetry is assumed
                    # Positive az's are to the right
                    # Use the same ear for both sides but flip the az. 
                    # So for a source az = +60, left = 60, right = 360-60 = 300
                    m_az = int(ild.rstrip("n"))
                    hrtf_r_ind = int((m_az % 360)/5)
                    hrtf_l_ind = int((360-(m_az % 360))/5)
                    
                    hrtf_l_data = dat[:, hrtf_l_ind]
                    hrtf_r_data = dat[:, hrtf_r_ind]

                    cue_data_l = ir.apply_mag_spec(cue_data,hrtf_l_data)
                    cue_data_r = ir.apply_mag_spec(cue_data,hrtf_r_data)
                    targ_data_l = ir.apply_mag_spec(targ_data,hrtf_l_data)
                    targ_data_r = ir.apply_mag_spec(targ_data,hrtf_r_data)
                    mask1_data_l = ir.apply_mag_spec(mask2_data,hrtf_r_data)
                    mask1_data_r = ir.apply_mag_spec(mask2_data,hrtf_l_data)
                    mask2_data_l = ir.apply_mag_spec(mask2_data,hrtf_l_data)
                    mask2_data_r = ir.apply_mag_spec(mask2_data,hrtf_r_data)

                    #mask1_data_l = np.convolve(mask1_data,hrtf_l_data)
                    #mask1_data_r = np.convolve(mask1_data,hrtf_r_data)
                else:
                    #mask1_data_l = mask1_data.copy()
                    this_ild = float(ild) / 2
                    cue_data_l = psylab.signal.atten(cue_data_l,-this_ild)
                    cue_data_r = psylab.signal.atten(cue_data_r,this_ild)
                    targ_data_l = psylab.signal.atten(targ_data_l,-this_ild)
                    targ_data_r = psylab.signal.atten(targ_data_r,this_ild)
                    mask1_data_l = psylab.signal.atten(mask1_data_l,this_ild)
                    mask1_data_r = psylab.signal.atten(mask1_data_r,-this_ild)
                    mask2_data_l = psylab.signal.atten(mask2_data_l,-this_ild)
                    mask2_data_r = psylab.signal.atten(mask2_data_r,this_ild)
                    #mask2_data_r = mask2_data.copy()


                targ_data_l, targ_data_r, mask1_data_l, mask1_data_r, mask2_data_l, mask2_data_r = psylab.signal.zeropad(
                targ_data_l, targ_data_r, mask1_data_l, mask1_data_r, mask2_data_l, mask2_data_r)

                # inf ild maskers good here 

                #plt.plot(targ_data_l)
                #plt.plot(targ_data_r)

                # # Can't combine and do normal enhancement, because targets
                # # and maskers must be vocoded separately, as per Ihlefeld
                #t_m_data_l = targ_data + mask1_data_l + mask2_data_l
                #t_m_data_r = targ_data + mask1_data_r + mask2_data_r
                #
                #
                #if mloc == '30e':
                #    sig_l_bp = psylab.signal.filter_bank(t_m_data_l, fs, 6, ild_bands)
                #    sig_r_bp = psylab.signal.filter_bank(t_m_data_r, fs, 6, ild_bands)
                #
                #    sig_out_l = np.zeros(t_m_data_l.shape)
                #    sig_out_r = np.zeros(t_m_data_r.shape)
                #
                #    for i in np.arange(sig_l_bp.shape[1]):
                #        band_l = sig_l_bp[:,i]
                #        band_r = sig_r_bp[:,i]
                #
                #        itd_track, cc = ild.gen_itd_track(band_l, band_r, fs=fs, 
                #                            max_itd_us=int(ild_max_itd_us),
                #                            wsize=ild_wsize, 
                #                            stepsize=ild_stepsize, 
                #                            thresh=0, 
                #                            return_cc=True,
                #                            )
                #
                #        gain_l, gain_r = ild.gen_gain_tracks(band_l.size, itd_track, ild_table_l, ild_table_r,
                #                            max_ild_db=ild_max_ild_db, 
                #                            max_itd_us=int(ild_max_itd_us), 
                #                            wsize=ild_wsize, 
                #                            stepsize=ild_stepsize, 
                #                            attack=ild_attack,
                #                            )
                #                            #nr_cc=cc, nr_thresh=nr_thresh[ii], nr_gain=nr_gain)
                #
                #        # TODO: Why are the gain tracks positive? Check gen_gain_tracks
                #        # AHA!: Look at line 385 in ild.py; mystery solved. FIXED (removed the negative signs in ild.py)
                #        sig_out_l += ild.apply_gain_track(band_l, gain_l)
                #        sig_out_r += ild.apply_gain_track(band_r, gain_r)
                #
                #    t_m_data_l = sig_out_l
                #    t_m_data_r = sig_out_r

            # Add location cue; combine target and masker
            sig_l = np.concatenate((cue_data_l, cue_isi_data, targ_data_l + mask1_data_l))
            sig_r = np.concatenate((cue_data_r, cue_isi_data, targ_data_r + mask1_data_r))
            if side == 'r':
                sig_out = np.column_stack((sig_r, sig_l))
            else:
                sig_out = np.column_stack((sig_l, sig_r))

            filename = f"{k}_{this_sf_path_out_sub}.wav"
            this_sf_filepath_out = os.path.join(this_sf_path_out, filename)
            this_sf_filepath_out_sub = os.path.join(this_sf_path_out_sub, filename)
            print(f"  {filename}")
            sfwrite(this_sf_filepath_out, sig_out, fs, format='wav')
            #targ_str += f"{','.join(targ_words)},{','.join(map(str, targ_times))}\n"
            #mask_str += f"{','.join(mask2_words)},{','.join(map(str, mask2_times))}\n"
            targ_str += ','.join([','.join(map(str, i)) for i in zip(targ_words, targ_times)]) + "\n"
            mask_str += ','.join([','.join(map(str, i)) for i in zip(mask2_words, mask2_times)]) + "\n"
            fid_all_times.write(targ_str)
            fid_all_times.write(mask_str)

            # Don't need these anymore, now using all_times
            #color_line_t = f"{this_sf_filepath_out_sub},{','.join(color_times_t)}"
            #fid_color_times_t.write(color_line_t+"\n")
            #
            #if len(color_times_m) > 0:
            #    color_line_m = f"{this_sf_filepath_out_sub},{','.join(color_times_m)}"
            #    fid_color_times_m.write(color_line_m+"\n")
            #
            #word_line = f"{this_sf_filepath_out_sub},{','.join(targ_words)}"
            #fid_target_words.write(word_line+"\n")


if __name__ == '__main__':

    subjectID = None
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hs:", ["help", "subjectID="])
    except (getopt.error, msg):
        print("for help use --help")
        sys.exit(2)
    for var, val in opts:
        if var in ("--help", "-h"):
            print(__doc__)
            sys.exit(0)
        elif var in ("--subjectID", "-s"):
            subjectID = val

    if not subjectID:
        subjectID = input('Enter subject ID: ')
    if subjectID == '':
        print("No subject ID found, exiting...")
        sys.exit(0)

    main(subjectID)
