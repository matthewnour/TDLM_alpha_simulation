# TDLM_alpha_simulation
Effect of background alpha on TDLM sequenceness periodicity and SNR (simulation implemented in MATLAB)

Simulations, illustrating effect of background MEG alpha oscillation on sequenceness periodicity and signal-noise ratio (SNR). 

Main wrapper script: 
run_simulate_replay.m
(Set code path in line 9, simPath)

Details of simulation: 
For each simulated participant we generated 1 minute of synthetic MEG data (272 sensors, sampled at 100 Hz, incorporating cross- and auto-correlation sensor relationships). 
Synthetic data contained (1) forward sequential reactivations of sensor patterns for 8 states at 50 ms lag, in a manner that replays a hypothetical task structure, as hypothesised in our task, and (2) an additive alpha-band amplitude modulation to the MEG time series (frequency ~ N(μ=10,σ=0.2)). We simulated 25 experiments (each of 20 participants) for 5 alpha ‘strength’ parameter levels, ranging from minimal (5 arbitrary units, a.u.) to extreme (25 a.u.). As in the main paper, for each experiment we run the full TDLM pipeline for each participant, controlling for 10 Hz oscillation, before averaging the sequenceness x lag effect over participants. For each experiment we then estimated (1) the power spectral density (PSD) of the group mean sequenceness x lag effect (using a discrete Fourier transform), and (2) the effect magnitude of the forward sequenceness effect at 50 ms lag (i.e. the ground truth effect).  

Results presented in multi-panel subplot, showing:
- Relationship between background alpha power and sequenceness periodicity (Mean and SEM over n=25 experiments).
- Relationship between background alpha power and signal-to-noise ratio (SNR) for ground truth forward sequenceness effect at 50 ms lag (mean and SEM over n=25 experiments). Effect normalized to experiment-specific significance threshold derived by permutation testing (Liu et al., 2019), marked by dashed horizontal line). 
- Mean power spectra of sequenceness x lag effect as a function of background alpha strength (each line is the mean of all subjects and all experiments). 

Matthew Nour, London, April 2021
