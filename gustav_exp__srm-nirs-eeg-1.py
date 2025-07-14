# -*- coding: utf-8 -*-

# Copyright (c) 2010-2012 Christopher Brown
#
# This file is part of Psylab.
#
# Psylab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Psylab is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Psylab.  If not, see <http://www.gnu.org/licenses/>.
#
# Bug reports, bug fixes, suggestions, enhancements, or other 
# contributions are welcome. Go to http://code.google.com/p/psylab/ 
# for more information and to contribute. Or send an e-mail to: 
# cbrown1@pitt.edu.
# 
# Psylab is a collection of Python modules for handling various aspects 
# of psychophysical experimentation. Python is a powerful programming  
# language that is free, open-source, easy-to-learn, and cross-platform, 
# thus making it extremely well-suited to scientific applications. 
# There are countless modules written by other scientists that are  
# freely available, making Python a good choice regardless of your  
# particular field. Consider using Python as your scientific platform.
# 

# A Gustav experiment file!

import os, sys
import time
import numpy as np
import matplotlib.pyplot as plt
import psylab
import gustav
#from gustav.forms.curs import rt as theForm
from gustav.forms import rt as theForm
import medussa as m
#import spyral
import timeit

try:
    import triggers
except Exception as e:
    ret = input("Triggers module or Cedrus cpod hardware not found. Continue? (y/N) ")
    if ret == 'y':
        triggers = None
    else: 
        raise e

def setup(exp):

    
    # Machine-specific settings
    machine = psylab.config.local_settings(conf_file='config/psylab.conf')
#    exp.user.machine_name = machine.get_str('name')
    #exp.user.pa_id = machine.get_list_int('audio_id')
    workdir = machine.get_path('workdir')
#    exp.stim.stimdir = machine.get_path('stimdir')
#    exp.user.pa_id = 7,7,2
    dev_name = machine.get_str('audiodev_name')
    dev_ch = machine.get_int('audiodev_ch')
    devs = m.get_available_devices()
    dev_id = None
    for i,di in enumerate(devs):
        name = psylab.string.as_str(di.name)
        if name.startswith(dev_name):
            dev_id = i,i,dev_ch
    if dev_id:
        exp.stim.audiodev = m.open_device(*dev_id)
    else:
        raise Exception(f"The audio device {dev_name} was not found")


#    exp.user.pa_id = 2,2,4
#    workdir = '/home/cbrown/work/Python'

    # General Experimental Variables
    exp.name = 'srm-nirs-eeg-1'
    exp.method = 'constant' # 'constant' for constant stimuli, or 'adaptive' for a staircase procedure (SRT, etc)
    # TODO: move logstring and datastring vars out of exp and into either method or experiment, so they can be properly enumerated at startup

    exp.logFile = os.path.join(workdir,'logs','$name_$date.log') # Name and date vars only on logfile name
    exp.logConsoleDelay = True
    exp.dataFile = os.path.join(workdir,'data','$name.csv')
    exp.recordData = True # Data is saved manually in present_trial
    exp.dataString_header = "# A datafile created by Gustav!\n# \n# Experiment: $name\n# \n# $note\n# \n# $comments\n# \n\nS,Trial,Date,Block,Condition,@currentvars[],Soundfile,Ear,Times\n"
    exp.dataString_post_trial = "$subj,$trial,$date,$block,$condition,$currentvars[],$stim[file],$user[side],$user[response]\n"
    exp.logString_pre_exp = "\nExperiment $name running subject $subj started at $time on $date\n"
    exp.logString_post_exp = "\nExperiment $name running subject $subj ended at $time on $date\n"
    exp.logString_pre_block = "\n  Block $block of $blocks started at $time; Condition: $condition ; $currentvarsvals[' ; ']\n"
    exp.logString_post_trial = "    Trial $trial, target stimulus: $user[trial_stimbase], KWs correct: $response / possible: $user[trial_kwp] ($user[block_kwc] / $user[block_kwp]: $user[block_pc] %)\n"
    exp.logString_post_block = "  Block $block of $blocks ended at $time; Condition: $condition ; $currentvarsvals[' ; ']\n"
    exp.frontend = 'tk'
    exp.debug = False
    # CUSTOM: A comma-delimited list of valid single-char responses. This experiment is designed to have 
    # the experimenter do the scoring, and enter the score on each trial.
    exp.validKeys = '0,1,2,3,4,5,6,7,8,9'.split(',')
    exp.quitKey = '/'
    exp.note = 'Vocoded speech in spatially-separated speech or noise maskers.'
    exp.comments = '''
    Intended for nirs data collection, to extend data from Zhang and Ihlefeld 2021.
    Replicates speech v noise conditions, and infinite ILDs (speech-oppo). 
    Adds 10, 20 & 30 dB ILDs. All conditions are symmetrical maskers.
    The other major change is stim dur = 24s, as opposed to 15, to account for 
    the slow change found in nirs-im-6, apparently from symmetrical maskers.
    The ask is to hit a key when a target color (not object) word is heard.
    '''

    exp.stim.breath_blocks = 5
    exp.stim.breath_block_breaths = 3
    exp.stim.hale_dur = 5
    exp.stim.hold_dur = 15
    exp.stim.atten = 10

    if not exp.subjID:
        exp.subjID = exp.term.get_input(parent=None, title = "Gustav!", prompt = 'Enter a Subject ID:')

    exp.stim.fs = 44100.
    exp.stim.basedir = os.path.join("stim","bobmike",exp.name,f"s_{exp.subjID}")
    exp.stim.stimfiles = {}
    exp.var.factorial['masker'] = []
    exp.stim.practfiles = psylab.folder.consecutive_files(
                    path=os.path.join(exp.stim.basedir,"practice"),
                    repeat=True,
                    file_ext=".WAV;.wav",
            )
    if triggers:
        print("Generating triggers:")
        exp.stim.trigger_dict = {'Inhale': 1, 'Exhale': 2, 'Hold_Breath': 3}
        ntriggers = 3
        exp.stim.trigfile = os.path.join(workdir,'data',f"{exp.name}_triggers.csv")
        if os.path.exists(exp.stim.trigfile):
            tf = open(exp.stim.trigfile, 'a+')
        else:
            tf = open(exp.stim.trigfile, 'a+')
            tf.write("S,Trig,Condition\n")

        exp.stim.trigfile_cond = os.path.join(workdir,'data',f"{exp.name}_{exp.subjID}_condition_order.csv")
        exp.stim.tf_c = open(exp.stim.trigfile_cond, 'a+')
        #exp.stim.tf_c.write("Condition\n")


        for f in os.scandir(exp.stim.basedir):
            if f.is_dir() and f.name != "practice":
                exp.stim.stimfiles[f.name] = psylab.folder.consecutive_files(
                        path=f.path,
                        file_ext=".WAV;.wav",
                )
                exp.var.factorial['masker'].append(f.name)
                ntriggers += 1
                exp.stim.trigger_dict[f.name] = ntriggers
                tf.write(f"{exp.subjID},{exp.var.factorial['masker'].index(f.name) + 4},{f.name}\n")

        for cond,n in exp.stim.trigger_dict.items():
            print(f"Trigger {n}: {cond}")
            tf.write(f"{exp.subjID},{n},{cond}\n")

    """EXPERIMENT VARIABLES
        There are 2 kinds of variables: factorial and covariable

        Levels added as 'factorial' variables will be factorialized with each
        other. So, if you have 2 fact variables A & B, each with 3 levels, you
        will end up with 9 conditions: A1B1, A1B2, A1B3, A2B1 etc..

        Levels added as 'covariable' variables will simply be listed (in parallel
        with the corresponding levels from the other variables) in the order
        specified. So, if you have 2 'covariable' variables A & B, each with 3
        levels, you will end up with 3 conditions: A1B1, A2B2, and A3B3. All
        'covariable' variables must have either the same number of levels, or
        exactly one level. When only one level is specified, that level will
        be used in all 'covariable' conditions. Eg., A1B1, A2B1, A3B1, etc.

        You can use both types of variables in the same experiment, but both
        factorial and covariable must contain exactly the same set of variable
        names. factorial levels are processed first, covariable levels are added
        at the end.

        Both factorial and covariable are Python ordered dicts, where the keys 
        are variable names, and the values are lists of levels. During the 
        experiment, you have access to the current level of each variable. For 
        example, if you have the following variable:
        
        exp.var.factorial['target'] = ['Male', 'Female']
        
        Then, you can find out what the level is at any point in the experiment 
        with exp.var.current['target'], which would return either 'Male' or 
        'Female' depending on what the condition happened to be. This is 
        probably most useful to generate your stimuli, eg., in the pre_trial 
        function. 
    """

    #exp.var.factorial['masker'] = [
    #                                  'm_speech__ild_0__ftype_bb',
    #                                  'm_speech__ild_10__ftype_bb',
    #                                  'm_speech__ild_20__ftype_bb',
    #                                  'm_speech__ild_30__ftype_bb',
    #                                  'm_speech__ild_inf__ftype_bb',
    #                                  'm_speech__ild_70__ftype_bb',
    #                                  'm_speech__ild_70__ftype_nat',
    #                                  'm_noise__ild_itd__ftype_bb',
    #                                 ]



    """CONSTANT METHOD VARIABLES
        The method of constant limits requires three variables to be set.
            trialsperblock
            startblock [crash recovery]
            starttrial [crash recovery]
    """
    exp.var.constant = {
        'trialsperblock' : 1,
        'startblock' : 1,
        'starttrial' : 1,
        }
    
    """CONDITION PRESENTATION ORDER
        Use 'prompt' to prompt for condition on each block, 'random' to randomize
        condition order, 'menu' to be able to choose from a list of conditions at
        the start of a run, 'natural' to use natural order (1:end), or a
        print-range style string to specify the order ('1-10, 12, 15'). You can
        make the first item in the print range 'random' to randomize the specified
        range.
    """
    #if np.random.randint(2) == 1:
    #    exp.var.order = ",".join([str(i) for i in np.tile(np.arange(len(exp.var.factorial['masker']))+1, 6)])
    #else:
    #    exp.var.order = ",".join([str(i) for i in np.flipud(np.tile(np.arange(len(exp.var.factorial['masker']))+1, 6))])

    order = np.arange(len(exp.var.factorial['masker']))+1
    np.random.shuffle(order)
    for i in range(4-1):
        this = np.arange(len(exp.var.factorial['masker']))+1
        looking = True
        while looking:
            np.random.shuffle(this)
            if this[0] != order[-1]:
                order = np.append(order, this)
                looking = False
    exp.var.order = ",".join(str(item) for item in list(order))


    """IGNORE CONDITIONS
        A list of condition numbers to ignore. These conditions will not be
        reflected in the total number of conditions to be run, etc. They will
        simply be skipped as they are encountered during a session.
    """
    exp.var.ignore = []

    """USER VARIABLES
        Add any additional variables you need here
    """
    # Placeholders. Values will be set in pre or post trial, then used for logging / saving data.
    exp.user.trial_kwp = 0
    exp.user.trial_stimbase = ''
    exp.user.block_kwp = 0.
    exp.user.block_kwc = 0.
    exp.user.block_pc = 0.
    if triggers:
        exp.user.triggers = triggers.xid()


def pre_exp(exp):
    try:
        exp.interface = theForm.Interface()
        exp.interface.update_Title_Center(exp.name)
        exp.interface.update_Title_Right(f"S {exp.subjID}", redraw=False)
        exp.interface.update_Prompt("Hit a key to begin", show=True, redraw=True)
        ret = exp.interface.get_resp()

        practice = True
        while practice:
            exp.interface.update_Prompt("Practice? (y/N)", show=True, redraw=True)
            ret = exp.interface.get_resp()
            if ret != 'y':
                exp.interface.update_Prompt("End practice", show=True, redraw=True)
                practice = False
            else:
                exp.interface.update_Prompt("Hit a key when you hear [red, green, blue, white]", show=True, redraw=True)
                exp.stim.out, exp.stim.fs = m.read_file(exp.stim.practfiles.get_filename(fmt='full'))
                s = exp.stim.audiodev.open_file(exp.stim.practfiles.get_filename(fmt='full'))
                s.mix_mat = psylab.signal.atten(s.mix_mat, exp.stim.atten)
                #exp.stim.out = psylab.signal.atten(exp.stim.out, exp.stim.atten)
                dur_ms = len(exp.stim.out) / exp.stim.fs * 1000
                resp_percent = []
                if not exp.debug:
                    #s = exp.stim.audiodev.open_array(exp.stim.out,exp.stim.fs)
                    this_wait_ms = 500
                    this_elapsed_ms = 0
                    s.play()
                    start_ms = exp.interface.timestamp_ms()
                    while s.is_playing:
                        ret = exp.interface.get_resp(timeout=this_wait_ms/1000)
                        this_current_ms = exp.interface.timestamp_ms()
                        this_elapsed_ms = this_current_ms - start_ms
                        this_elapsed_percent = this_elapsed_ms / dur_ms * 100
                        if ret:
                            resp_percent.append(this_elapsed_percent)
                            #responses.append(str(this_elapsed_ms/1000))

                        progress = psylab.string.prog(this_elapsed_percent, width=50, char_done="=", spec_locs=resp_percent, spec_char="X")
                        exp.interface.update_Prompt(progress, show=True, redraw=True)
                    #exp.interface.show_Notify_Left(show=False, redraw=True)
            #    m.play_array(stim.out,exp.stim.fs) #,output_device_id=exp.user.audio_id)
    #            exp.interface.showPlaying(False)

        if not ret == exp.quitKey:
            exp.interface.update_Prompt("Breathing Exercise? (y/N)", show=True, redraw=True)
            ret = exp.interface.get_resp()
            if ret != 'y':
                exp.interface.update_Prompt("No Breathing Exercise", show=True, redraw=True)
            else:
                for i in range(exp.stim.breath_blocks):
                    exp.interface.update_Prompt(f"{exp.stim.hale_dur} sec inhale, {exp.stim.hale_dur} sec exhale ({exp.stim.breath_block_breaths} times), then {exp.stim.hold_dur} sec breath hold\n({i+1} of {exp.stim.breath_blocks}; hit a key to start)", show=True, redraw=True)
                    ret = exp.interface.get_resp()
                    if ret == exp.quitKey:
                        exp.run.gustav_is_go = False
                    else:
                        for j in range(exp.stim.breath_block_breaths):
                            hale_cur = 0
                            prompt = f"Inhale ({j+1}/{exp.stim.breath_block_breaths})..."
                            time_init = exp.interface.timestamp_ms()
                            if triggers:
                                # exp.stim.trigger_dict = {'Inhale': 1, 'Exhale': 2, 'Hold_Breath': 3}
                                exp.interface.update_Status_Right(f"trigger {exp.stim.trigger_dict['Inhale']}")
                                exp.user.triggers.trigger(exp.stim.trigger_dict['Inhale'])
                            while hale_cur < exp.stim.hale_dur:
                                hale_cur = np.minimum((exp.interface.timestamp_ms() - time_init) / 1000, exp.stim.hale_dur)
                                hale_cur = np.maximum(hale_cur, 0)
                                progress = psylab.string.prog(hale_cur, width=50, char_done="=", maximum=exp.stim.hale_dur)
                                this_prompt = f"{prompt}\n{progress} {np.int32(hale_cur)} / {exp.stim.hale_dur}"
                                exp.interface.update_Prompt(this_prompt, show=True, redraw=True)
                                time.sleep(.2)
                            hale_cur = exp.stim.hale_dur
                            prompt = f"Exhale ({j+1}/{exp.stim.breath_block_breaths})..."
                            time_init = exp.interface.timestamp_ms()
                            if triggers:
                                # exp.stim.trigger_dict = {'Inhale': 1, 'Exhale': 2, 'Hold_Breath': 3}
                                exp.interface.update_Status_Right(f"trigger {exp.stim.trigger_dict['Exhale']}")
                                exp.user.triggers.trigger(exp.stim.trigger_dict['Exhale'])
                            while hale_cur > 0:
                                hale_cur = np.minimum(exp.stim.hale_dur - ((exp.interface.timestamp_ms() - time_init) / 1000), exp.stim.hale_dur)
                                hale_cur = np.maximum(hale_cur, 0)
                                progress = psylab.string.prog(hale_cur, width=50, char_done="=", maximum=exp.stim.hale_dur)
                                this_prompt = f"{prompt}\n{progress} {exp.stim.hale_dur - np.int32(hale_cur)} / {exp.stim.hale_dur}"
                                exp.interface.update_Prompt(this_prompt, show=True, redraw=True)
                                time.sleep(.2)
                        hale_cur = 0
                        prompt = f"Inhale (then hold)..."
                        time_init = exp.interface.timestamp_ms()
                        if triggers:
                            # exp.stim.trigger_dict = {'Inhale': 1, 'Exhale': 2, 'Hold_Breath': 3}
                            exp.interface.update_Status_Right(f"trigger {exp.stim.trigger_dict['Inhale']}")
                            exp.user.triggers.trigger(exp.stim.trigger_dict['Inhale'])
                        while hale_cur < exp.stim.hale_dur:
                            hale_cur = np.minimum((exp.interface.timestamp_ms() - time_init) / 1000, exp.stim.hale_dur)
                            hale_cur = np.maximum(hale_cur, 0)
                            progress = psylab.string.prog(hale_cur, width=50, char_done="=", maximum=exp.stim.hale_dur)
                            this_prompt = f"{prompt}\n{progress} {np.int32(hale_cur)} / {exp.stim.hale_dur}"
                            exp.interface.update_Prompt(this_prompt, show=True, redraw=True)
                            time.sleep(.2)
                        hold_cur = exp.stim.hold_dur
                        prompt = f"Hold ({exp.stim.hold_dur} sec)..."
                        time_init = exp.interface.timestamp_ms()
                        if triggers:
                            # exp.stim.trigger_dict = {'Inhale': 1, 'Exhale': 2, 'Hold_Breath': 3}
                            exp.interface.update_Status_Right(f"trigger {exp.stim.trigger_dict['Hold_Breath']}")
                            exp.user.triggers.trigger(exp.stim.trigger_dict['Hold_Breath'])
                        while hold_cur > 0:
                            hold_cur = np.minimum(exp.stim.hold_dur - ((exp.interface.timestamp_ms() - time_init) / 1000), exp.stim.hold_dur)
                            hold_cur = np.maximum(hold_cur, 0)
                            progress = psylab.string.prog(hold_cur, width=50, char_done="=", maximum=exp.stim.hold_dur)
                            this_prompt = f"{prompt}\n{progress} {exp.stim.hold_dur - np.int32(hold_cur)} / {exp.stim.hold_dur}"
                            exp.interface.update_Prompt(this_prompt, show=True, redraw=True)
                            time.sleep(.2)
                        hale_cur = exp.stim.hale_dur
                        prompt = f"Exhale ({j+1}/{exp.stim.breath_block_breaths})..."
                        time_init = exp.interface.timestamp_ms()
                        if triggers:
                            # exp.stim.trigger_dict = {'Inhale': 1, 'Exhale': 2, 'Hold_Breath': 3}
                            exp.interface.update_Status_Right(f"trigger {exp.stim.trigger_dict['Exhale']}")
                            exp.user.triggers.trigger(exp.stim.trigger_dict['Exhale'])
                        while hale_cur > 0:
                            hale_cur = np.minimum(exp.stim.hale_dur - ((exp.interface.timestamp_ms() - time_init) / 1000), exp.stim.hale_dur)
                            hale_cur = np.maximum(hale_cur, 0)
                            progress = psylab.string.prog(hale_cur, width=50, char_done="=", maximum=exp.stim.hale_dur)
                            this_prompt = f"{prompt}\n{progress} {exp.stim.hale_dur - np.int32(hale_cur)} / {exp.stim.hale_dur}"
                            exp.interface.update_Prompt(this_prompt, show=True, redraw=True)
                            time.sleep(.2)

        if ret == exp.quitKey:
            #exp.interface.update_Prompt("Hit a key when you hear [red, green, blue, white]", show=False, redraw=True)
            exp.gustav_is_go = False

    except Exception as e:
        exp.interface.destroy()
        raise e


def pre_block(exp):
    try:
        exp.user.block_kwp = 0
        exp.user.block_kwc = 0
        exp.user.block_pc = 0.
        exp.user.pract = 1
        exp.interface.update_Status_Left(f"Block {exp.run.block+1} of {exp.run.nblocks}")
    except Exception as e:
        exp.interface.destroy()
        raise e

"""PRE_TRIAL
    This function gets called on every trial to generate the stimulus, and do
    any other processing you need. All settings and variables are available. 
    For the current level of a variable, use exp.var.current['varname'].
"""
def pre_trial(exp):
    try:
        exp.stim.file = exp.stim.stimfiles[exp.var.current['masker']].get_filename(fmt='full')
        exp.stim.trigger = exp.stim.trigger_dict[exp.var.current['masker']]
        exp.interface.update_Status_Center(exp.var.current['masker'], redraw=True) # Use condition # (1,2) as trigger #

        exp.stim.tf_c.write(f"{exp.var.current['masker']}\n")

        out,fs = m.read_file(exp.stim.file)
        out = psylab.signal.atten(out, exp.stim.atten)
        trigs = audio_triggers(out.shape[0], fs, {0: 1}, nbits = 1, click_len = 2, zeropad_cols = 0)
        exp.stim.out = np.hstack((out, trigs))

    except Exception as e:
        exp.interface.destroy()
        raise e


def present_trial(exp):
    # This is a custom present_trial that records keypress times during playback

    exp.interface.update_Status_Right(f"Trigger {exp.stim.trigger}", redraw=True) # Use condition # (1,2) as trigger #
    exp.interface.update_Prompt("Hit [L/R/B] to start", show=True, redraw=True)
    wait = True
    while wait:
        ret = exp.interface.get_resp()
        if ret in ['l', 'r', 'b', exp.quitKey]:
            wait = False

    if ret == exp.quitKey:
        exp.run.gustav_is_go = False
    else:
        try:
            exp.user.side = ret
            exp.interface.update_Prompt("Hit a key when you hear [blue, red, green, white]", show=True, redraw=True)
            #exp.interface.update_Notify_Left('Playing', show=False, redraw=True)
            #exp.interface.update_Status_Right("", redraw=True)
            responses = []
#            exp.interface.showPlaying(True)
            if not exp.debug:
                s = exp.stim.audiodev.open_array(exp.stim.out,exp.stim.fs)
                #exp.interface.show_Notify_Left(show=True, redraw=True)
                dur_ms = len(exp.stim.out) / exp.stim.fs * 1000
                this_wait_ms = 500
                this_elapsed_ms = 0
                resp_percent = []
                s.play()
                if triggers:
                    exp.user.triggers.trigger(exp.stim.trigger)
                start_ms = exp.interface.timestamp_ms()
                while s.is_playing:
                    ret = exp.interface.get_resp(timeout=this_wait_ms/1000)
                    this_current_ms = exp.interface.timestamp_ms()
                    this_elapsed_ms = this_current_ms - start_ms
                    this_elapsed_percent = this_elapsed_ms / dur_ms * 100
                    if ret:
                        responses.append(str(np.round(this_elapsed_ms/1000, 3)))
                        resp_percent.append(this_elapsed_ms / dur_ms * 100)

                    progress = psylab.string.prog(this_elapsed_percent, width=50, char_done="=", spec_locs=resp_percent, spec_char="X")
                    exp.interface.update_Prompt(progress, show=True, redraw=True)
                exp.user.response = ",".join(responses)
#                fid = open(exp.dataFile, 'a')
#                word_line = f"{exp.stim.file},{','.join(responses)}"
#                fid.write(word_line+"\n")
#            exp.interface.showPlaying(False)
        except Exception as e:
            exp.interface.destroy()
            raise e


#    m.play_array(stim.out,exp.stim.fs) #,output_device_id=exp.user.audio_id)
#    exp.interface.showPlaying(False)

"""CUSTOM PROMPT
    If you want a custom response prompt, define a function for it
    here. exp.run.response should receive the response as a string, and
    if you want to cancel the experiment, set both exp.run.block_on and
    exp.run.pylab_is_go to False
"""

def prompt_response(exp):
    pass
    # while True:
    #     ret = exp.interface.get_resp()
    #     if ret in exp.validKeys:
    #         exp.run.response = ret
    #         break
    #     elif ret in ['/','q']:
    #         exp.run.gustav_is_go = False
    #         break

def post_trial(exp):
    #if not exp.gustav_is_go:
    st = np.random.randint(13, 18)
    exp.interface.update_Prompt(f"Waiting {st} sec...", show=True, redraw=True)
    time.sleep(st)
    exp.interface.update_Prompt("", show=True, redraw=True)
#    try:
#        if exp.run.gustav_is_go:
#            exp.user.block_kwp += int(exp.user.trial_kwp)
#            exp.user.block_kwc += int(exp.run.response)
#            exp.user.block_pc = round(float(exp.user.block_kwc) / float(exp.user.block_kwp) * 100, 1)
#            trial_pc = round(float(exp.run.response) / float(exp.user.trial_kwp) * 100, 1)
#            #t = 'Trial: %s / %s (%s %%)\n' % (exp.run.response, exp.user.trial_kwp, trial_pc)
#            t = f"Trial: {exp.run.response} / {exp.user.trial_kwp}\n"
#        else:
#            t = 'Trial unfinished (Exp cancelled)\n'
#        #blockPercent = round(exp.user.block_kwc / exp.user.block_kwp * 100, 1)
#        #exp.interface.updateInfo_TrialScore(t + 'Total: %0.0f / %0.0f (%0.1f %%)' % 
#        #(exp.user.block_kwc, exp.user.block_kwp, exp.user.block_pc))
#    #    exp.interface.updateInfo_TrialScore(f"{t}Total: {exp.user.block_kwc:.1f} / {exp.user.block_kwp:.1f}")
#    except Exception as e:
#        exp.interface.destroy()
#        raise e

def post_block(exp):
    pass
#    exp.interface.updateInfo_BlockScore(f"Prev Condition # {exp.run.condition+1}\nScore: {exp.user.block_kwc:.1f} / {exp.user.block_kwp:.1f} ({exp.user.block_pc} %%)")

def post_exp(exp):
#    pass
#    exp.interface.dialog.isPlaying.setText("Finished")
#    exp.interface.showPlaying(True)
    exp.interface.destroy()

def byte_to_bits(byte, nbits=4):
    # Returns a list of bools (1's or 0's) indicating which bits are set in the specified byte
    # Returned list will be of length bits 
    return [int(i) for i in f"{{0:0{nbits}b}}".format(byte)]

def audio_triggers(n, fs, triggers, nbits = 4, click_len = 2, zeropad_cols = 0):
    """Generates an array of ones and zeros for use as triggers delivered as audio clicks
        (ie., suitable for Triggy from Cortech). 

        The idea is to generate your audio tracks, then generate triggers and append to the audio

        Note: To reiterate, this doesn't actually send triggers, it only generates an array
        of audio clicks to be sent to Triggy which will interpret them to generate the triggers

        Parameters
        ----------
        dur_s : float
            The duration of the audio, in seconds
        fs : int
            The sampling frequency
        triggers : dict
            Specifies which triggers to fire and when. Keys are times (floats) and vals as 
            trigger numbers (ints).
        nbits : int
            The number of bits available. If you are using Triggy with the audio adapter,
            you have 4 channels and thus 4 bits. This means you have 15 triggers available.
        click_len : int
            The length of the clicks, in samples. 
        zeropad_cols : int
            If the audio channels devoted to triggers are the last available (eg., channels
            5,6,7,8) and you are only using the first few for stimulus delivery (1,2), then
            you can zeropad the trigger array so that it can be easily appended to the audio
            array (in this example, zeropad_cols == 2). 

        Returns
        -------
        out : array
            An array of length dur_s * fs and width nbits, that specifies which triggers to 
            fire and when.  

    """
    out = np.zeros((n, nbits))
    for time, trigger in triggers.items():
        if trigger > 2**nbits-1:
            raise ValueError(f'Trigger {trigger} cannot be represented; nbits=={nbits} so only 2**nbits (0-{2**nbits-1}) triggers are available')
        sample = np.int32(np.round(np.float32(time)*np.float32(fs)))
        bits = byte_to_bits(trigger, nbits=nbits)
        bits.reverse()
        for n in range(len(bits)):
            out[sample:sample+click_len, n] = -bits[n]
    if zeropad_cols > 0:
        pad = np.atleast_2d(np.zeros((n, zeropad_cols)))
        out = np.hstack((pad, out))
    return out


if __name__ == '__main__':
    argv = sys.argv[1:]
    argv.append(f"--experimentFile={os.path.realpath(__file__)}")
    gustav.gustav.main(argv)

