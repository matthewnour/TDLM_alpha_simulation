% Simulations, illustrating effect of background MEG alpha oscillation on sequenceness periodicity and signal-noise ratio (SNR). 
%
% Matthew Nour, London, April 2021

close all
clear all
clc

simPath = '';
addpath(simPath)

% alpha strengths per experiment
alphas = [5:5:25];

% number of experimetns each strength
nExp = 25; % each of 20 subjects

sqn_store = nan(nExp, length(alphas), 20, 60);  % [experiment, alphaStrength, subj, lag]
GroupThresh_store = nan(nExp, length(alphas));  % [experiment, alphaStrength]
psd_store = nan(nExp, length(alphas), 20, 31);  % [experiment, alphaStrength, subj, freq]

for nE = 1:nExp
    disp(['Experiment: ' num2str(nE)])
    tic
    
    for ind = 1:length(alphas)
        
        sf = []; sb = [];
        
        alphaStrength = alphas(ind);
        
        simulate_replay
        
        % store the subj-specific fwd-bwd at this experiment and alphaStrength
        dtp = squeeze(sf(:,1,2:end)-sb(:,1,2:end));
        sqn_store(nE, ind, :, :) = dtp; % [experiment, alphaStrength, ground truth, lag]
        npThresh = squeeze(max(abs(mean(sf(:,2:end,2:end)-sb(:,2:end,2:end),1)),[],3));
        GroupThresh_store(nE, ind, :, :) = max(npThresh);
        
        % subj-specific spectrogram
        L = size(dtp,2);
        fft_freq_domain = [0:samplerate/L:samplerate/2];
        pxx1 = abs(fft(dtp', L));
        pxx1 = (1/(samplerate*L)) * pxx1.^2; % 'Power'
        pxx1 = pxx1(1:L/2+1, :);
        pxx1(2:end-1, :) = 2*pxx1(2:end-1, :);
        
        psd_store(nE, ind, :, :) = pxx1';
        
    end
    
    toc
end %exp
disp('done')

%% Plot
simulation_settings = {};
simulation_settings.nSensors = nSensors;
simulation_settings.nTrainPerStim = nTrainPerStim;
simulation_settings.nNullExamples = nNullExamples;
simulation_settings.nSamples = nSamples;
simulation_settings.nSequences = nSequences;
simulation_settings.maxLag = maxLag;
simulation_settings.cTime = cTime;
simulation_settings.nSubj = nSubj;
simulation_settings.gamA = gamA;
simulation_settings.gamB = gamB;
simulation_settings.pInds = pInds;
simulation_settings.samplerate = samplerate;
simulation_settings.nstates = nstates;
simulation_settings.bins = bins;
simulation_settings.TF = TF;

summary_plots

%% ----------------------------------------------------------
function squashed = squash(x)
    squashed = x(:);
end
