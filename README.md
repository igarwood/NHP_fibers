# Data for: Multifunctional fibers enable modulation of cortical and deep brain activity during cognitive behavior in macaques
---

Simultaneous neuronal recordings and neuromodulation were enabled by multifunctional fiber neurotechnology. This dataset contains neurophysiological and behavioral data recorded from a macaque performing a working memory task. Recordings were performed in the premotor cortex or putamen. During the recordings, there was intracranial neuromodulation. 

There are 5 recording types and corresponding files
- pmc_gaba: Recordings in the premotor cortex with intracranial GABA infusion(s) (n = 4)
- pmc_sal: Recordings in the premotor cortex with intracranial saline infusion(s) (n = 3)
- pmc_acsf: Recordings in the premotor cortex with intracranial aCSF infusion(s) (n = 2)
- pmc_noic: Recordings in the premotor cortex with no intracranial infusion (n = 4)
- put_GABA: Recording in the putamen with intracranial GABA infusion
- pmc_noic: Recordings in the putamen with no intracranial infusion.

Device data is stored in device_data.zip (see README_code.md for more information)

Additionally, we provide preprocessed data in experimental_data.zip and processed_data.zip (see README_code.md for more information)

Further detail regarding experimental details can be found in the associated manuscript (Garwood, et al.)


## Description of experimental data and file structure

Each data file contains 2 variables and 4 structures:
Variables:
- fs_lfp: The sampling rate of local field potential (lfp) recordings
- fs_spike: The sampling rate of single unit activity recordings
Structures:
- infusion (when applicable; empty otherwise): Data related to the intracranial infusions
    - start: The start time(s) of intracranial infusion(s) in seconds, relative to the beginning of the recording
    - end: The end time(s) of intracranial infusion(s) in seconds, relative to the beginning of the recording
    - rate: The infusion rate(s) of intracranial infusion(s) in nL/min
    - drug: String corresponding to the drug infused
- lfp: Data related to local field potential oscillations
    - lfp: 4xN matrix of LFP data
    - lfp_elec: The LFP electrode used for subsequent analysis
- sua: Data related to single unit activity
    - unit_locs: Indices, relative to the beginning of the recording, corresponding to timepoints when each (sorted) unit fired an action potential
    - unit_spikes: The corresponding (unsorted) event index for each sorted event, relative to the first event index
    - spike_locs_all: The indices of all valid spike events, relative to the beginning of the recording
    - spikes_all: The (unsorted) waveforms of all spike events, in chronological order
    - spike_length: The duration of each spike waveform (in ms)
    - spike_threshold: The threshold used to define spike events, relative to the standard deviation of the raw spike data at baseline
    - noise_thresh: The amplitude beyond which threshold crossings were considered to be an artifact (in uV)
 - task_info: Data related to the behavioral task
     - fs_task: Sampling rate used to convert trial indices to time (in sec)
     - trials: Trial start and end indices, relative to the beginning of the recording
     - samples: Start and stop indices corresponding to when the sample was visible on the screen
     - matches: Start and stop indices corresponding to when the match was visible on the screen.
     - correct_trials: Binary variable, where 0 indicates incorrect trials, 1 indicates correct trials
     - complete_trials: Binary variable, where 0 indicates incomplete trials, 1 indicates complete trials
     - no_fix: Binary variable, where 0 indicates successful fixation, 1 indicates unsuccessful fixation
     - sample_id: The ID of the sample (1-3)
     - correct_loc: The location of the correct selection
     - reaction time: The time between when the match appeared and when a selection was made (in s)
  
Note that one aCSF file was included in behavioral analysis but not electrophysiology analysis, because the infusion rate + volume exceeded that of the GABA recordings.

## Code/Software
Code will be available in a permanently archived repository located at https://github.com/igarwood/NHP_fibers

See the main README_code.md file in https://github.com/igarwood/NHP_fibers

MIT license. 
