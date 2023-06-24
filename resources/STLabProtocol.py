#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to run monocular and binocular targeted photoreceptor experiments.


@author: jtm
"""
import os
import os.path as op
import shutil
import glob
from time import sleep
import pickle

import pandas as pd
import parallel
from pyplr.pupil import PupilCore
from pyplr.stlabscene import SpectraTuneLabScene
from pyplr.protocol import timer

import make_stimuli


#-------------------------------#
#    Experiment configuration   #
#-------------------------------#
# Experiments to choose from
EXPERIMENTS_TO_RUN = dict(
    enumerate(['S-cone', 'M-cone', 'L-cone', 'Melanopsin', 'L-M', 'Luminance']))

# Experiment pauses at end of trial if this file exists
PAUSEFILE = op.join(op.expanduser('~'), 'PAUSEFILE')

# Timing parameters in s
ADAPTATION_PERIOD = 2*60
INCLUDE_ADAPTATION_PERIOD = True
STIMULUS_DURATION = 12

# Experiment parameters
DOUBLE_TRIALS = False

# Set to 1. to run all 60 trials
FRAC_TRIALS = 1.0

# inter-trial interval
ITI = 3

# Hardware flags
USE_STLAB = True
USE_PUPIL_CORE = True
USE_EEG = False

#-------------------------------#
# Functions


class OverwriteError(Exception):
    pass


def choose_experiment():
    print('Which experiment to run?\n')
    print(*EXPERIMENTS_TO_RUN.items(), sep='\n')
    experiment = ''
    while experiment not in EXPERIMENTS_TO_RUN.keys():
        experiment = int(input('\n> Enter the number: '))
        if experiment not in EXPERIMENTS_TO_RUN.keys():
            print('Enter a valid experiment...')
    return EXPERIMENTS_TO_RUN[experiment]

#-------------------------------#


def main():
    try:
        print('\n')
        print(f"{'='*80 : ^80}")
        print(f"{'MonBin Photoreceptors': ^80}")
        print(f"{'='*80 : ^80}")
        print('\n')
        sleep(.1)

        # Ask experimenter to choose the experiment
        experiment = choose_experiment()

        # Ask for subject ID and block number
        subject_id = input('> Enter subject_id: ')
        block = input('> Enter block: ')

        # Experiment results folder
        results_folder = f'./{experiment}/results/'
        if not op.exists(results_folder):
            os.makedirs(results_folder)
            print(f'> Created results folder: {results_folder}')

        # Experiment stimulus folder
        stim_folder = f'./{experiment}/stims/'
        if not op.exists(stim_folder):
            os.makedirs(stim_folder)
            print(f'Created stims folder: {stim_folder}')

        # Subject directory
        subject_dir = op.join(results_folder, f'{subject_id}')
        if not op.exists(subject_dir):
            os.mkdir(subject_dir)
            print(f'> Created subject folder: {subject_dir}')
        else:
            print(f'> Subject folder already exists: {subject_dir}')

        # Stimuli are adjusted to account for age-related changes in spectral
        # sensitivity of the photoreceptors. If the stimuli for a given
        # experiment and observer age do not exist, call the script and create
        # them. This doesn't happen for the luminance condition.
        subject_age = int(input("> Enter subject age: "))
        if experiment=='Luminance':
            pass
        
        elif experiment=='L-M':
            stim_out_dir = op.join(stim_folder, str(subject_id))
            if not op.exists(stim_out_dir):
                print(f"> No stimuli for age {subject_id}")
                make_stimuli.main(stim_out_dir, experiment, subject_age)
            else:
                print(f"Stimuli exist for age {subject_id}")
        else:
            stim_out_dir = op.join(stim_folder, str(subject_age))
            if not op.exists(stim_out_dir):
                print(f"> No stimuli for age {subject_age}")
                make_stimuli.main(stim_out_dir, experiment, subject_age)
            else:
                print(f"Stimuli exist for age {subject_age}")

        # Do not allow silent overwrite of existing data.
        if op.exists(op.join(subject_dir, f'block_{block}_trial_log.csv')):
            raise OverwriteError(
                f'\nWARNING: The log file "block_{block}_trial_log.csv" '
                'already exists for this subject. Delete manually if you '
                'really wish to overwrite.\n'
            )
        
        if experiment=='L-M':
            stim_folder = op.join(stim_folder, str(subject_id))
        elif not experiment == 'Luminance':
            stim_folder = op.join(stim_folder, str(subject_age))
        else:
            pass
        print(f"Stimulus folder: {stim_folder}")


        # Background spectra for each device
        with open('./calibration/STLAB_1_background.pkl', 'rb') as fp:
            STLAB_1_BACKGROUND = pickle.load(fp)
        with open('./calibration/STLAB_2_background.pkl', 'rb') as fp:
            STLAB_2_BACKGROUND = pickle.load(fp)

        # The experiment will pause at the end of a trial if ~/PAUSEFILE
        # exists, and can then be resumed by hitting ENTER in the ipython
        # console. To pause the experiment, type 'pause' in a terminal (aliased
        # to "touch ~/PAUSEFILE") or manually create the file. The file is
        # removed when the experiment resumes, but here we check just in case.
        if os.path.isfile(PAUSEFILE):
            print(f'> Removing {PAUSEFILE}.')
            os.remove(PAUSEFILE)

        # Open LPT1 or /dev/parport0
        if USE_EEG:
            POST_SCENE_LAUNCH_CODE = 99
            pport = parallel.Parallel()
            pport.setData(0)
            print('> Opened parallel port and set data pins to 0x00')

        # Connect to Pupil Core and make sure Annotation Capture plugin is
        # active
        if USE_PUPIL_CORE:
            pcore = PupilCore()
            pcore.annotation_capture_plugin(should='start')
            print('> Connected to Pupil Core and enabled Annotation Capture')
 
        # Connect to STLAB and authenticate
        if USE_STLAB:
            d = SpectraTuneLabScene.from_config()
            print('\n')

            # Multicast addresses used for playing video files. These should
            # already be configured, but if not, they can be set manually with
            # d.set_multicast_address(...). The broadcast address targets all
            # luminaires
            STLAB_1, STLAB_2 = 1021, 1022
            STLAB_BROADCAST = 1023

            # Half max background for adaptation
            d.set_spectrum_a(STLAB_1_BACKGROUND, STLAB_1)
            print('> Setting background spectrum on STLAB:')
            print(f'\t{STLAB_1_BACKGROUND}')

            d.set_spectrum_a(STLAB_2_BACKGROUND, STLAB_2)
            print('> Setting background spectrum on STLAB:')
            print(f'\t{STLAB_2_BACKGROUND}')

        # Start recording with Pupil Core
        if USE_PUPIL_CORE:
            pcore.command(f'R {op.abspath(subject_dir)}')

        # Get conditions, double up and randomise
        trials = pd.read_csv('datasource.csv')

        if not DOUBLE_TRIALS:
            trials = (trials.sample(frac=FRAC_TRIALS)
                      .reset_index(drop=True))
        else:
            trials = (pd.concat([trials, trials])
                      .sample(frac=FRAC_TRIALS)
                      .reset_index(drop=True))

        # Save the trial log
        trials.to_csv(op.join(subject_dir, f'block_{block}_trial_log.csv'))
        with open(op.join(subject_dir, 'age.txt'), 'w') as fh:
            fh.write(f'Subject: {subject_id}, Age: {subject_age}')

        # Adaptation peridod
        if USE_STLAB and INCLUDE_ADAPTATION_PERIOD:
            _ = input("Press Enter to start adaptation period: ")
            timer(1, ADAPTATION_PERIOD, '> Background spectrum adaptation...')

        # Go straight on to experiment
      #  _ = input("Press Enter to start experiment: ")

        # Main trial loop starts here
        for trial_num, trial in trials.iterrows():
            print(f"{'='*80 : ^80}")
            print(f"{'TRIAL ' + str(trial_num) : ^80}")
            print(f"{'='*80 : ^80}")
            print(f"> Ocular condition: {trial.ocular_condition}")
            print(f"> Condition code: {trial.condition_code}")
            print('\n')
            print(f"{'Left eye' : ^40}{'Right eye' : ^40}")
            print(f"{'--------' : ^40}{'---------' : ^40}")
            print(f"{trial.videoL : ^40}{trial.videoR: ^40}")
            print(f"{str(trial.contrastL) + '% contrast' : ^40}"
                  f"{str(trial.contrastR) + '% contrast' : ^40}")
            print(f"{str(trial.frequencyL) + ' Hz frequency' : ^40}"
                  f"{str(trial.frequencyR) + ' Hz frequency' : ^40}")
            print('\n')

            # Video files are kept in the stims folder with informative names.
            # When required, we copy them into the current working directory
            # and change the names to video1.json and video2.json, which is the
            # required format for multiple uploads.
            # Let's always follow the convention below:
            # video1.json plays on STLAB_1, which stimulates the left eye
            # video2.json plays on STLAB_2, which stimulates the right eye

            # First, get rid of any old video files in the current working dir.
            for f in glob.glob('./*.json'):
                os.remove(f)
                print(f'> Removed {f} from the current working directory.')

            # Now copy accross the required files and rename. Note that in
            # MonL/MonR conditions only one video file is played, so we catch
            # the error and set the value to None, which means that a video
            # will not be played for that eye when we issue the scene command.
            try:
                video_1 = shutil.copyfile(
                    op.join(stim_folder, trial.videoL), 'video1.json')
            except Exception:
                video_1 = None
            try:
                video_2 = shutil.copyfile(
                    op.join(stim_folder, trial.videoR), 'video2.json')
            except Exception:
                video_2 = None

            if USE_STLAB:
                # Upload and cache the video files.
                for video_file in [video_1, video_2]:
                    if video_file is not None:
                        d.upload_video(video_file)
                        _ = d.get_video_file_metadata(video_file)
                print('> Video file(s) uploaded and cached on LIGHT HUB.')

            else:
                print('> Video files would be uploaded here...')

            # Send trigger to Pupil Capture and EEG before launching scene
            # command
            if USE_PUPIL_CORE:
                custom_fields = {
                    'condition_code': trial.condition_code,
                    'ocular_condition': trial.ocular_condition,
                    'contrastL': trial.contrastL,
                    'contrastR': trial.contrastR,
                    'frequencyL': trial.frequencyL,
                    'frequencyR': trial.frequencyR,
                    'videoL': trial.videoL,
                    'videoR': trial.videoR
                }
                pre_scene_annotation = pcore.new_annotation(
                    label='pre_scene',
                    custom_fields=custom_fields
                )
                pcore.send_annotation(pre_scene_annotation)

            if USE_EEG:
                pport.setData(trial.condition_code)

            #-----------------------------------#
            # Launch the video file(s) on STLAB #
            #-----------------------------------#
            if USE_STLAB:
                # Launch the video files
                d.scene(STLAB_1, STLAB_2, video_1, video_2)
            else:
                print('> Scene command would be used to launch video'
                      + 'files here...')

            # Send trigger to Pupil Capture and EEG when scene command has
            # returned
            if USE_EEG:
                pport.setData(POST_SCENE_LAUNCH_CODE)

            if USE_PUPIL_CORE:
                post_scene_annotation = pcore.new_annotation(
                    label='post_scene',
                    custom_fields=custom_fields
                )
                pcore.send_annotation(post_scene_annotation)

            # Wait until stimulus has finished
            timer(1, STIMULUS_DURATION, '> Administering stimulus...')

            if USE_EEG:
                pport.setData(0)
                print('> Reset parallel port data pins to 0x00.')

            print(f"{'='*80 : ^80}")
            print(f"{'END OF TRIAL ' + str(trial_num) : ^80}")
            print(f"{'='*80 : ^80}")
            print('\n')
            timer(1, ITI, '> Intertrial interval...')
            print('\n')

            # Check for pause
            if os.path.isfile(PAUSEFILE):
                input('> Experiment paused. Hit ENTER to continue...')
                print(f'> Removing {PAUSEFILE} and resuming experiment.')
                os.remove(PAUSEFILE)

    except KeyboardInterrupt:
        print('> Experiment terminated by user.')

    except OverwriteError as e:
        print(e)

    except Exception as e:
        print(e)
        print('> Experiment not working.')
        print('\t> Is Pupil Capture running and the eye tracker connected?')
        print('\t> Is STLAB / LIGHT HUB receiving power and connected?')
        print('\t> Is the parallel port connected? ')

    finally:
        print('\n')
        print(f"{'='*80 : ^80}")
        print(f"{'END OF EXPERIMENT': ^80}")
        print(f"{'='*80 : ^80}")
        print('\n')
        print('Cleaning up...')

        for f in glob.glob('./*.json'):
            os.remove(f)
            print(f'> Removed {f} from the current working directory.')

        if USE_STLAB:
            d.stop_video(STLAB_1)
            d.stop_video(STLAB_2)
            d.turn_off(STLAB_BROADCAST)
            d.logout()
            print('> Turned off STLAB and logged out of the LIGHT HUB.')

        if USE_PUPIL_CORE:
            pcore.command('r')
            print('> Stopped recording in Pupil Capture.')

        if USE_EEG:
            pport.setData(0)
            print('> Reset parallel port data pins to 0x00.')


if __name__ == '__main__':
    main()
