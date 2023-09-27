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

Device data is stored in device_data.zip 

Additionally, we provide preprocessed data in experimental_data.zip and processed_data.zip 

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

We additionally provide the output of fooof analysis in experimental_data.zip/data_for_fooof
  - The data files have the following naming convention: fooof_output_exptype_6.mat where 6 corresponds to the FOOOF parameter for the max number of spectral peaks in each spectra
  - Each data file contains fooof parameters estimated from 2N spectral samples, where samples 1:N occurred before intracranial infusions and samples N+1:2N occurred after intracranial infusions (Garwood 2023, Methods)
  - Each data file contains 6 variables:
      - aperiodic_signal: two columns corresponding to the aperiodic exponent and aperiodic offset across 2N samples
      - fooof_spectrum: 2Nxn_freq matrix corresponding to the estimated spectra of 2N samples
      - fooof_freq: n_freq vector of frequencies corresponding to the columns in fooof_spectrum (where frequencies between 59.5 and 60.5 Hz have been removed)
      - peak_cf: 2Nx6 matrix of center frequencies across up to 6 spectral peaks; if peak_cf(n,p+1:6) = 0, only p spectral peaks were identified in that sample
      - peak_pow: 2Nx6 matrix of power across up to 6 spectral peaks; if peak_pow(n,p+1:6) = 0, only p spectral peaks were identified in that sample
      - session_tags: 2Nx1 vector containing the session number corresponding to each spectral sample (also saved in the files session_tags_exptype.mat)
  - freqs.csv
      - contains the frequencies of the original spectral samples

## Description of device data and file structure

Device data is organized into three subfolders: DMA, fluidic, and impedance.
- DMA: Dynamic material analysis data
  - Contains DMA data for three fiber samples and one steel sample in .csv format.
- fluidic: fluidic rate characterization data
  - Contains 30 .mat files in the format of fluidic infusion data from 10-100 nl/min
    - Naming convention: rate_nlm_trial_x.mat where "nlm" refers to nanoliters/minute as the units of time
    - Each file contains two variables: time, and volume (cumulative volume infused over time).
  - Subfolder 'raw_videos' contains the video recordings corresponding to each trial in the parent directory.
- impedance: impedance spectroscopy data
  - Contains 24 .txt files containing impedance spectroscopy data (frequency vs. impedance)
    - Naming convention: devx_elecy.txt or devx_postAC_elecy.txt (postAC indicates spectroscopy performed after autoclave sterilization)

## Description of preprocessed data and file structure

Processed data file includes spike sorting results and example SS-GLM and AR models
- spike_data: folder containing spike sorting results
  - Naming convention: spikes_cluster_location_exptypex.mat where 'x' indicates that experiment number
  - Each file contains 3 variables:
    - s_time: experiment time sampled at 30 kHz
    - spikes_cluster: the spike cluster id of every threshold crossing
    - spikes_locs: n_units by 1 cell array containing the index of every spike of the corresponding unit. Single unit times are given by s_time(spike_locs{unit_num})
- SSGLM models follow the naming convention location_exptypex_unity_trialvariant_Rz.mat where 'x' indicates the experiment number, 'y' indicates the unit id, and 'z' indicates the model order.
    - Additional subscripts include '_nohist.mat' indicating the model has no history terms and '_stationary.mat', indicating that the parameters are stationary across trials.
    - Each file contains estimated SSGLM parameters (Garwood 2023, Methods, README_code.md)
- AR models follow the naming convention location_exptypex_lfp_trialvariant_Rz.mat

## Code/Software
Code will be available in a permanently archived repository located at https://github.com/igarwood/NHP_fibers

See the main README_code.md file in https://github.com/igarwood/NHP_fibers

MIT license. 
