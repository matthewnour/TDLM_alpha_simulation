%% Specs
TF = [0,1,0,0,0,0,0,0;0,0,1,0,0,0,0,0;0,0,0,1,0,0,0,0;0,0,0,0,0,0,0,0;0,0,0,0,0,1,0,0;0,0,0,0,0,0,1,0;0,0,0,0,0,0,0,1;0,0,0,0,0,0,0,0]; % transition matrix
TR = TF';

nSubj = 20;                         % number of subjects to simulate

nSensors = 272;
nTrainPerStim = 18;                 % how many training examples for each stimulus
nNullExamples = nTrainPerStim*8;    % how many null examples to use
nSamples = 6000;                    % 60 seconds of unlabelled data to predict
nSequences = 2000;                  % how many real sequences to put in the data
maxLag = 60;                        % evaluate time lag up to 600ms
gamA = 10;
gamB = 0.6;
[~, pInds] = uperms([1:8],29);
uniquePerms=pInds;
nShuf = size(uniquePerms,1);
samplerate=100;
nstates=8;
bins=10;                            % mean alpha wavelength

sf = cell(nSubj,1);  sb = cell(nSubj,1);

%% Core function
parfor iSj = 1:nSubj
    
    sf{iSj} = nan(1, nShuf, maxLag+1);
    sb{iSj} = nan(1, nShuf, maxLag+1);
      
    % generate dependence of the sensors
    A = randn(nSensors);
    [U,~] = eig((A+A')/2);
    covMat = U*diag(abs(randn(nSensors,1)))*U';
    
    % generate the true patterns
    commonPattern = randn(1,nSensors);
    patterns = repmat(commonPattern, [8 1]) + randn(8, nSensors);
    
    % make training data
    trainingData = 4*randn(nNullExamples+8*nTrainPerStim, nSensors) + [zeros(nNullExamples,nSensors); ...
        repmat(patterns, [nTrainPerStim 1])];
    trainingLabels = [zeros(nNullExamples,1); repmat((1:8)', [nTrainPerStim 1])];
    
    % train classifiers on training data
    betas = nan(nSensors, 8); intercepts = nan(1,8);
    
    for iC=1:8
        [betas(:,iC), fitInfo] = lassoglm(trainingData, trainingLabels==iC, 'binomial', 'Alpha', 1, 'Lambda', 0.006, 'Standardize', false);
        intercepts(iC) = fitInfo.Intercept;
    end
    
    % make long unlabelled data with or without real sequences in it
    X = nan(nSamples, nSensors);
    X(1,:) = randn([1 nSensors]);
    
    for iT=2:nSamples
        X(iT,:) = 0.95*(X(iT-1,:) + mvnrnd(zeros(1,nSensors), covMat));
    end
    
    % Add alpha amplitude modulation
    for isen=1:nSensors
        % jitter and offset on alpha frequency, freq ~ N(10, 0.2)
        freq = bins + randn(1)*0.2;
        X(:,isen) = X(:,isen) + alphaStrength*[cos(2*pi*freq/samplerate*(1:length(X))+2*pi*rand(1)-pi)]';
    end
    
    % Injecting Sequence
    for iRS = 1:nSequences
        seqTime = randi([40 nSamples-40]);
        state = false(8,1); state(randi(8)) = true;
        
        for iMv=1:2
            if sum(state)==0
                X(seqTime,:) = X(seqTime,:);
            else
                X(seqTime,:) = X(seqTime,:) + patterns(state,:);
                state = (state'*TF)'; state2 = false(8,1); state2(find(rand < cumsum(state), 1, 'first')) = true; state = state2; % advance states
                seqTime = seqTime + round(gamrnd(gamA,gamB));
            end
        end
    end
      
    % make predictions with trained models
    preds = 1./(1+exp(-(X*betas + repmat(intercepts, [nSamples 1]))));
    
    % calculate sequenceness (TDLM)
    for iShuf = 1:nShuf
        rp = uniquePerms(iShuf,:);
        T1 = TF(rp,rp); T2 = T1';
        X=preds;
        
        nbins=maxLag+1;
        
        warning off % toeplitz warning
        dm=[toeplitz(X(:,1),[zeros(nbins,1)])];
        dm=dm(:,2:end);
        
        for kk=2:nstates
            temp=toeplitz(X(:,kk),[zeros(nbins,1)]);
            temp=temp(:,2:end);
            dm=[dm temp];
        end
        
        warning on
        
        Y=X;
        betas = nan(nstates*maxLag, nstates);
        
        % 1st level GLM (controlling for coactivations and alpha correction)
        for ilag=1:bins
            temp_zinds = (1:bins:nstates*maxLag) + ilag - 1;
            temp = pinv([dm(:,temp_zinds) ones(length(dm(:,temp_zinds)),1)])*Y;
            betas(temp_zinds,:)=temp(1:end-1,:);
        end
        
        % 2nd level GLM (evidnce for transitions of interest)
        betasnbins64=reshape(betas,[maxLag nstates^2]);
        bbb=pinv([T1(:) T2(:) squash(eye(nstates)) squash(ones(nstates))])*(betasnbins64');
        
        sf{iSj}(1,iShuf,2:end) = bbb(1,:);
        sb{iSj}(1,iShuf,2:end) = bbb(2,:);
  
    end
end

sf = cell2mat(sf);
sb = cell2mat(sb);




