#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy3 Experiment Builder (v3.0.7),
    on Fri Aug  9 18:28:16 2019
If you publish work using this script please cite the PsychoPy publications:
    Peirce, JW (2007) PsychoPy - Psychophysics software in Python.
        Journal of Neuroscience Methods, 162(1-2), 8-13.
    Peirce, JW (2009) Generating stimuli for neuroscience using PsychoPy.
        Frontiers in Neuroinformatics, 2:10. doi: 10.3389/neuro.11.010.2008
"""

from __future__ import absolute_import, division
from psychopy import locale_setup, sound, gui, visual, core, data, event, logging, clock
from psychopy.constants import (NOT_STARTED, STARTED, PLAYING, PAUSED,
                                STOPPED, FINISHED, PRESSED, RELEASED, FOREVER)
import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy import (sin, cos, tan, log, log10, pi, average,
                   sqrt, std, deg2rad, rad2deg, linspace, asarray)
from numpy.random import random, randint, normal, shuffle
import os  # handy system and path functions
import sys  # to get file system encoding
import time
import collections
import zmq
import pupil_socket

# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))
os.chdir(_thisDir)

# Store info about the experiment session
psychopyVersion = '3.0.7'
expName = 'CRFstudy'  # from the Builder filename that created this script
expInfo = {'participant': '', 'session': '001'}
dlg = gui.DlgFromDict(dictionary=expInfo, title=expName)
if dlg.OK == False:
    core.quit()  # user pressed cancel
expInfo['date'] = data.getDateStr()  # add a simple timestamp
expInfo['expName'] = expName
expInfo['psychopyVersion'] = psychopyVersion

# Create zmq socket
socket = pupil_socket.ZMQsocket()
socket.connect()
socket.start_recording(expInfo['participant'])
time.sleep(5)

# Data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
filename = _thisDir + os.sep + u'data/%s_%s_%s' % (expInfo['participant'], expName, expInfo['date'])

# An ExperimentHandler isn't essential but helps with data saving
thisExp = data.ExperimentHandler(name=expName, version='',
    extraInfo=expInfo, runtimeInfo=None,
    originPath='/Users/danbaker/Google Drive/Current work/Research/Melatonin/CRFstudy.py',
    savePickle=True, saveWideText=True,
    dataFileName=filename)
# save a log file for detail verbose info
logFile = logging.LogFile(filename+'.log', level=logging.EXP)
logging.console.setLevel(logging.WARNING)  # this outputs to the screen, not a file

endExpNow = False  # flag for 'escape' or other condition => quit the exp

# Start Code - component code to be run before the window creation

# Setup the Window
win = visual.Window(
    size=(800, 600), fullscr=False, screen=1, 
    winType='pyglet', allowGUI=False, allowStencil=False,
    monitor='Iiyama60Hz800x600', color=[-1, -1, -1], colorSpace='rgb',
    blendMode='avg', useFBO=True, 
    units='height')

# store frame rate of monitor if we can measure it
expInfo['frameRate'] = win.getActualFrameRate()
if expInfo['frameRate'] != None:
    frameDur = 1.0 / round(expInfo['frameRate'])
else:
    frameDur = 1.0 / 60.0  # could not measure, so guess

conditionlistL = np.array([0.06, 0.12, 0.24, 0.48, 0.96, 0.0, 0.0, 0.0, 0.0, 0.0,
                           0.06, 0.12, 0.24, 0.48, 0.96, 0.06, 0.12, 0.24, 0.48, 0.96,
                           0.06, 0.12, 0.24, 0.48, 0.96, 0.48, 0.48, 0.48, 0.48, 0.48,
                           0.06, 0.12, 0.24, 0.48, 0.96, 0.0, 0.0, 0.0, 0.0, 0.0,
                           0.06, 0.12, 0.24, 0.48, 0.96, 0.06, 0.12, 0.24, 0.48, 0.96,
                           0.06, 0.12, 0.24, 0.48, 0.96, 0.48, 0.48, 0.48, 0.48, 0.48])       # target opacities (e.g. contrasts)
conditionlistR = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.06, 0.12, 0.24, 0.48, 0.96,
                           0.06, 0.12, 0.24, 0.48, 0.96, 0.06, 0.12, 0.24, 0.48, 0.96,
                           0.48, 0.48, 0.48, 0.48, 0.48, 0.06, 0.12, 0.24, 0.48, 0.96, 
                           0.0, 0.0, 0.0, 0.0, 0.0, 0.06, 0.12, 0.24, 0.48, 0.96,
                           0.06, 0.12, 0.24, 0.48, 0.96, 0.06, 0.12, 0.24, 0.48, 0.96,
                           0.48, 0.48, 0.48, 0.48, 0.48, 0.06, 0.12, 0.24, 0.48, 0.96])       # target opacities (e.g. contrasts)
# stimulus frequencies in Hz for each condition
freqlistL = np.array([2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 
                      2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 
                      2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 
                      1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 
                      2.0, 2.0, 2.0, 2.0, 2.0, 1.6, 1.6, 1.6, 1.6, 1.6,
                      2.0, 2.0, 2.0, 2.0, 2.0, 1.6, 1.6, 1.6, 1.6, 1.6])
freqlistR = np.array([2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 
                      2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 
                      2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 
                      1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 
                      1.6, 1.6, 1.6, 1.6, 1.6, 2.0, 2.0, 2.0, 2.0, 2.0, 
                      1.6, 1.6, 1.6, 1.6, 1.6, 2.0, 2.0, 2.0, 2.0, 2.0])
    
# Initialize components for Routine "trial"
trialClock = core.Clock()

sound_1 = sound.Sound('A', secs=0.01, octave=1, hamming=False, stereo=True)
sound_1.setVolume(1)

left = visual.GratingStim(
    win=win, name='left',
    tex='sin', mask='circle',
    ori=0, pos=(-0.25, 0), size=(0.25, 0.25), sf=0, phase=0.0,
    color=[1,1,1], colorSpace='rgb', opacity=0.5,blendmode='avg',
    texRes=512, interpolate=True, depth=0.0)
right = visual.GratingStim(
    win=win, name='right',
    tex='sin', mask='circle',
    ori=0, pos=(0.25, 0), size=(0.25, 0.25), sf=0, phase=0.0,
    color=[1,1,1], colorSpace='rgb', opacity=0.5,blendmode='avg',
    texRes=512, interpolate=True, depth=-1.0)
Lcross = visual.ShapeStim(
    win=win, name='Lcross', vertices='cross',
    size=(0.02, 0.02),
    ori=0, pos=(-0.25, 0),
    lineWidth=1, lineColor=[-1,-1,-1], lineColorSpace='rgb',
    fillColor=[-1,-1,-1], fillColorSpace='rgb',
    opacity=1, depth=-2.0, interpolate=True)
Rcross = visual.ShapeStim(
    win=win, name='Rcross', vertices='cross',
    size=(0.02, 0.02),
    ori=0, pos=(0.25, 0),
    lineWidth=1, lineColor=[-1,-1,-1], lineColorSpace='rgb',
    fillColor=[-1,-1,-1], fillColorSpace='rgb',
    opacity=1, depth=-3.0, interpolate=True)

#TF = 2
ifims = 1000*frameDur
nframes = 12*1/frameDur
frametimes = np.linspace(0,12000,int(nframes+1))
#targetwaveform = np.sin(2 * TF * frametimes * pi/1000)

# Create some handy timers
globalClock = core.Clock()  # to track the time since experiment started
routineTimer = core.CountdownTimer()  # to track time remaining of each (non-slip) routine 
trialcounter = -1

condlist = []
for x in range(1,61):
    condlist.append(collections.OrderedDict([('cond',x)]))

# set up handler to look after randomisation of conditions etc
trials = data.TrialHandler(nReps=1, method='random', 
    extraInfo=expInfo, originPath=-1,
    trialList=condlist,
    seed=None, name='trials')
thisExp.addLoop(trials)  # add the loop to the experiment
thisTrial = trials.trialList[0]  # so we can initialise stimuli with some values
# abbreviate parameter names if possible (e.g. rgb = thisTrial.rgb)
if thisTrial != None:
    for paramName in thisTrial:
        exec('{} = thisTrial[paramName]'.format(paramName))

for thisTrial in trials:
    currentLoop = trials
    
    
    # abbreviate parameter names if possible (e.g. rgb = thisTrial.rgb)
    if thisTrial != None:
        for paramName in thisTrial:
            exec('{} = thisTrial[paramName]'.format(paramName))
            
    trialcounter = trialcounter + 1
    # ------Prepare to start Routine "trial"-------
    t = 0
    trialClock.reset()  # clock
    frameN = -1
    continueRoutine = True
    routineTimer.add(15.000000)
    delayframes = (3*round(expInfo['frameRate']))   # 3 second delay
    # update component parameters for each repeat
    # keep track of which components have finished
    trialComponents = [left, right, Lcross, Rcross, sound_1]
    for thisComponent in trialComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    
    left.opacity = 0.5
    right.opacity = 0.5
    
    # -------Start Routine "trial"-------
    while continueRoutine and routineTimer.getTime() > 0:
        # get current time
        t = trialClock.getTime()
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        if frameN>delayframes:
            left.opacity = (conditionlistL[cond-1]*np.sin(freqlistL[cond-1]*(frameN/round(expInfo['frameRate']))*2*pi)+1)/2
            right.opacity = (conditionlistR[cond-1]*np.sin(freqlistR[cond-1]*(frameN/round(expInfo['frameRate']))*2*pi)+1)/2
        
        # *left* updates
        if t >= 0.0 and left.status == NOT_STARTED:
            # keep track of start time/frame for later
            left.tStart = t
            left.frameNStart = frameN  # exact frame index
            left.setAutoDraw(True)
        frameRemains = 0.0 + 15.0- win.monitorFramePeriod * 0.75  # most of one frame period left
        if left.status == STARTED and t >= frameRemains:
            left.setAutoDraw(False)
        
        # *right* updates
        if t >= 0.0 and right.status == NOT_STARTED:
            # keep track of start time/frame for later
            right.tStart = t
            right.frameNStart = frameN  # exact frame index
            right.setAutoDraw(True)
        frameRemains = 0.0 + 15.0- win.monitorFramePeriod * 0.75  # most of one frame period left
        if right.status == STARTED and t >= frameRemains:
            right.setAutoDraw(False)
        
        # *Lcross* updates
        if t >= 0.0 and Lcross.status == NOT_STARTED:
            # keep track of start time/frame for later
            Lcross.tStart = t
            Lcross.frameNStart = frameN  # exact frame index
            Lcross.setAutoDraw(True)
        frameRemains = 0.0 + 15- win.monitorFramePeriod * 0.75  # most of one frame period left
        if Lcross.status == STARTED and t >= frameRemains:
            Lcross.setAutoDraw(False)
        
        # *Rcross* updates
        if t >= 0.0 and Rcross.status == NOT_STARTED:
            # keep track of start time/frame for later
            Rcross.tStart = t
            Rcross.frameNStart = frameN  # exact frame index
            Rcross.setAutoDraw(True)
        frameRemains = 0.0 + 15- win.monitorFramePeriod * 0.75  # most of one frame period left
        if Rcross.status == STARTED and t >= frameRemains:
            Rcross.setAutoDraw(False)
        
        
        # start/stop sound_1
        if t >= 3.0 and sound_1.status == NOT_STARTED:
            # keep track of start time/frame for later
            sound_1.tStart = t  # not accounting for scr refresh
            sound_1.frameNStart = frameN  # exact frame index
            win.timeOnFlip(sound_1, 'tStartRefresh')  # time at next scr refresh
            win.callOnFlip(sound_1.play)  # screen flip
            frameRemains = 3.0 + 0.01- win.monitorFramePeriod * 0.75  # most of one frame period left
            if sound_1.status == STARTED and t >= frameRemains:
                # keep track of stop time/frame for later
                sound_1.tStop = t  # not accounting for scr refresh
                sound_1.frameNStop = frameN  # exact frame index
                win.timeOnFlip(sound_1, 'tStopRefresh')  # time at next scr refresh
                if 0.1 > 0.5:  # don't force-stop brief sounds
                    sound_1.stop()
            
            
        # check for quit (typically the Esc key)
        if endExpNow or event.getKeys(keyList=["escape"]):
            core.quit()
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in trialComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
            if frameN==delayframes:
                expInfo['trialonset'] = time.time()
                expInfo['condition'] = cond
                expInfo['leftcond'] = conditionlistL[cond-1]
                expInfo['rightcond'] = conditionlistR[cond-1] 
                expInfo['leftfreq'] = freqlistL[cond-1]
                expInfo['rightfreq'] = freqlistR[cond-1] 
                 
    # -------Ending Routine "trial"-------
    for thisComponent in trialComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    thisExp.nextEntry()
    sound_1.stop()
    
# completed 6 repeats of 'trials'
    
win.flip()

# these shouldn't be strictly necessary (should auto-save)
thisExp.saveAsWideText(filename+'.csv')
thisExp.saveAsPickle(filename)
logging.flush()

socket.stop_recording()

# make sure everything is closed down
thisExp.abort()  # or data files will save again on exit
win.close()

