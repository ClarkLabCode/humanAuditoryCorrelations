%% script to perform simple analysis of the correlation structures in human speech
%% this script requires the signal processing, image processing, and statistics and machine learning toolboxes

% sample English phrase for example figures
P = loadFile('/quoteSnippets/cummings_trimmed.m4a');
% change path to wherever this file is...

% calculate motion signal

P = extractMotion(P);

P.verbose = 1;
P.HRCtauHP = 2; % 50 ms
P = extractMotionEstimators(P);

P.HRCtauHP = 0; % just derivatives
P = extractMotionEstimators(P);

plotSpectrogram(P);

%% *************************
%% *************************
%% *************************
%% *************************
%% *************************
%% *************************
%% *************************
%% BELOW IS FOR FIGURES!
%% loading in a bunch of English data...

dEng = chooseAudioFiles('/Volumes/2023 SSD/LibriSpeech');
% this path should point towards the Librispeech home directory
% originally accessed from this link: https://www.openslr.org/12/

PmultEng = loadMultipleChosenFiles(dEng);
%%
for ii=1:length(PmultEng)
    temp = extractMotion(PmultEng(ii));
    temp.HRCtauHP = 1;
    Pmult2(ii) = temp; % do this with simplest HRC
    Pmult3(ii) = extractMotionEstimators(Pmult2(ii));
    if ~mod(ii,50)
        disp(['done with file ' num2str(ii)]);
    end
end
PmultEng = Pmult3; clear temp Pmult2 Pmult3;

QEng = compilePmultVals(PmultEng);
%%
regressMotionEstimators(QEng);

%% %%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%
%% loading in a bunch of Mandarin data...

dMan = chooseAudioFiles('/Users/DAC77/Documents/MATLAB/testing/test_ToneCorrelations/soundScapeAnalysis/test');
% this path should point to the home directory of the Mandarin speech
% corpus, retrieved from here: https://www.openslr.org/123/

PmultMan = loadMultipleChosenFiles(dMan);

for ii=1:length(PmultMan)
    temp = extractMotion(PmultMan(ii));
    temp.HRCtauHP = 0;
    Pmult2(ii) = temp;
    Pmult3(ii) = extractMotionEstimators(Pmult2(ii));
    if ~mod(ii,50)
        disp(['done with file ' num2str(ii)]);
    end
end
PmultMan = Pmult3; clear temp Pmult2 Pmult3;


QMan = compilePmultVals(PmultMan);

%%
regressMotionEstimators(QMan);

%% save all the figures!
if 1
    saveOpenFigures('savedFigs');
end



function Q = compilePmultVals(Pmult)
% function concatenates all the values for velocity estimates, etc through
% all elements in Pmult, creating a new single value structure Q containing
% those estimates, which can be used to plot histograms and regressions
% over the complete dataset...

    f = fieldnames(Pmult);
    Q = Pmult(1); % to get things started
    
    for ii=2:length(Pmult) % concatonate all the fields into one...
        P = Pmult(ii); % the one add
        currT = size(Q.power,2); % where we are now
        addT = size(P.power,2); % what we'll add
        ind1 = currT + 1;
        ind2 = currT + addT; 
        Q.power(:,ind1:ind2) = P.power;
        Q.powerBin(:,ind1:ind2) = P.powerBin;
        Q.uField(:,ind1:ind2) = P.uField;
        Q.uEst(ind1:ind2) = P.uEst;
        Q.timepoints(ind1:ind2) = P.timepoints;
        f2 = fieldnames(P.estimators);
        for jj=1:length(f2)
            f3 = fieldnames(P.estimators.(f2{jj}));
            for kk = 1:length(f3)
                Q.estimators.(f2{jj}).(f3{kk})(:,ind1:ind2) = P.estimators.(f2{jj}).(f3{kk});
            end
        end
        if ~mod(ii,100)
            disp(['concatonated ' num2str(ii) '...']);
        end
    end

end

function P = loadFile(folderFileName)

    [s,Fs] = audioread(folderFileName);
    P = extractSpectrogram(s,Fs);
    P.fileName = folderFileName;

end

function d = chooseAudioFiles(baseFolder)

    s = rng; % save rng settings
    rng(0); % so that it's choosing things at random in a reproducible way

    d = dir([baseFolder '/**/*.wav']); % goes through all subfolders
    if ~length(d) % in case it's not wav and instead flac 
        d = dir([baseFolder '/**/*.flac']); 
    end
    if ~length(d)
        disp('can''t find audio files');
        return;
    end

    keepInds = [];
    for ii=1:length(d) % clean out these stupid files starting with ._
        if ~strcmp('._',d(ii).name(1:2))
            keepInds = [keepInds,ii];
        end
    end
    d = d(keepInds); 
    
    % now, choose some set of them
    inds = randperm(length(d));
    loadSize = 1e8; % in bytes
    locSize = 0;
    for ii=1:length(inds)
        locSize(ii) = d(inds(ii)).bytes;
    end
    cumSize = cumsum(locSize);
    chooseNum = max(find(cumSize<=loadSize)); % how many to load in
    d = d(inds(1:chooseNum));
    
    rng(s); % return back the original rng settings
end

function Pmult = loadMultipleChosenFiles(d)

    % loads from dir structure d
    for ii=1:length(d)
        [s,Fs]=audioread([d(ii).folder '/' d(ii).name]);
        Pmult(ii) = extractSpectrogram(s,Fs);
        if ~mod(ii,50)
            disp(['file number ' num2str(ii) '...']);
        end
    end
    for ii=1:length(d)
        Pmult(ii).fileName = [d(ii).folder '/' d(ii).name];
    end
    disp(['loaded data from ' num2str(length(d)) ' files']);

end

function Pmult = loadMultipleFiles(folderName)

    d = dir([folderName '/*.wav']);
    if ~length(d)
        d = dir([folderName '/*.flac']);
    end
    if ~length(d)
        disp('unable to find audio files');
        return;
    end
    for ii=1:length(d)
        [s,Fs]=audioread([folderName '/' d(ii).name]);
        Pmult(ii) = extractSpectrogram(s,Fs);
        if ~mod(ii,50)
            disp(['file number ' num2str(ii) '...']);
        end
    end
    for ii=1:length(d)
        Pmult(ii).fileName = [folderName '/' d(ii).name];
    end
    disp(['loaded data from ' num2str(length(d)) ' files']);

end

function regressMotionEstimators(P)
    % function to do some regression to estimate the velocity from
    % correlations of various sorts
    totalTime = size(P.power,2)*(P.timepoints(2)-P.timepoints(1));
    disp(['Total time = ' num2str(totalTime/60) ' minutes']);

    %% first, do only the binary pairwise regressors!
    mhh = mean(P.estimators.binary.hhnet,1);
    mll = mean(P.estimators.binary.llnet,1);
    mhl = mean(P.estimators.binary.hlnet,1);
    mlh = mean(P.estimators.binary.lhnet,1);
    mfull = mean(P.estimators.binary.fullnet,1);

    doRegressionPlotting(P.uEst',[mhh',mll',mhl',mlh'],...
        {'binhh ','binll ','binhl ','binlh '},1e2); % since these tend to sum to 0, they need regularization


    %% second, HRC full
    hrc = mean(P.estimators.HRC.HRC,1);


    %% third, do four quadrants of HRC

    hrcpp = mean(P.estimators.HRC.HRCpp,1);
    hrcmm = mean(P.estimators.HRC.HRCmm,1);
    hrcpm = mean(P.estimators.HRC.HRCpm,1);
    hrcmp = mean(P.estimators.HRC.HRCmp,1);

    doRegressionPlotting(P.uEst',[hrcpp',hrcmm',hrcpm',hrcmp'],...
        {'HRCpp ','HRCmm ','HRCpm ','HRCmp '},0);
    
    function doRegressionPlotting(y,X,xnames,meth)
        % meth = 0 for standard regression, equal to regularizer for ridge
        % regression

        gridSize = ceil(sqrt(size(X,2))); % for plotting stuff

        rng(1);
        Nmax = 5e3;
        if length(y)<=Nmax
            ch = [1:Nmax];
        else
            ch = randi(length(y),[Nmax,1]);
        end

        figure;
        for ii=1:size(X,2)
            % only plot a fraction of all points for each plot
            

            subplot(gridSize,gridSize,ii); hold on;
            plot(y(ch),X(ch,ii),'k.');
            plot(10*[-1 1],[0 0],'k:');
            plot([0 0],10*[-1 1],'k:');
            xlabel('measured vel (a.u.)');
            ylabel(xnames{ii});
            C = corr(y,X(:,ii));
            title(['\rho = ' num2str(C)]);
            set(gca,'xlim',max(abs(y(ch)))*[-1 1]*1.1,'ylim',max(max(abs(X(ch,:))))*[-1 1]*1.1);
            box off; 
        end
        
        if 1 % don't plot this up for now
            if meth == 0
                [m,mint,r,rint,stats] = regress(y,[ones(size(X,1),1),X]);
                yhat = [ones(size(X,1),1),X]*m;
    %             disp(num2str(size(m)));
            else
                m = ridge(y,X,meth,0);
                yhat = m(1) + X*m(2:end);
            end
            C = corr(y,yhat);
            
            figure;
            plot(y(ch),yhat(ch),'.');
            xlabel('measured vel (oct/s, prob a.u.)');
            ylabel('full prediction');
            title(['\rho = ' num2str(C) '; predictors = ' [xnames{:}] '; m = ' num2str(m')]);
        end


    end
    

end

function P = extractMotionEstimators(P)
    % creates a bunch of motion estimators at each pixel in the image
    
    % regularize power so volume differences never matter for snippets...
    P.power = P.power - min(P.power(:));
    P.power = P.power ./ max(P.power(:)); % in limits of 0 to 1

    % FIRST: binary ones
    bw = imbinarize(P.power); % uses Otsu's method to find threshold that minimizes within class variance
    P.powerBin = bw;
    
    deltaT = 1;
    deltaF = 1;

    hhnet = matchCorrNet(bw,deltaF,deltaT,1,1);
    llnet = matchCorrNet(bw,deltaF,deltaT,0,0);
    hlnet = matchCorrNet(bw,deltaF,deltaT,1,0);
    lhnet = matchCorrNet(bw,deltaF,deltaT,0,1);

    P.estimators.binary.hhnet = hhnet;
    P.estimators.binary.llnet = llnet;
    P.estimators.binary.hlnet = hlnet;
    P.estimators.binary.lhnet = lhnet;
    P.estimators.binary.fullnet = hhnet + llnet - hlnet - lhnet;

    % SECOND: a set of continuous ones, using bandpass filters on the
    % stimulus power
    if ~isfield(P,'HRCtauHP')
        % this is default behavior if P doesn't have this field
        tauHP = find(P.timepoints - P.timepoints(1) == 0.05); % method fails if this interval doesn't exist.
        tauHP = tauHP - 1; % it's off by one using this method.
        [Bhi,Ahi] = butter(1,1/pi/(tauHP),'high'); % 0.05 s high pass filter to get changes in power.
        [Blo,Alo] = butter(1,1/pi/tauHP/2,'low'); % 0.1 s low pass filter to do the delay

    elseif P.HRCtauHP > 0
        % if it's given and >0, use it as tauHP
        tauHP = P.HRCtauHP;
        [Bhi,Ahi] = butter(1,1/pi/(tauHP),'high'); % chosen high pass filter to get changes in power.
        [Blo,Alo] = butter(1,1/pi/tauHP/2,'low'); % twice as long low pass filter to do the delay

    else
        % if it's set to 0, just do simplest thing
        Bhi = [1 -1]; Ahi = 1;
        Blo = [0 1]; Alo = 1;
    end
    % combining low with high pass gives a band pass filter...
    powerHigh = filter(Bhi,Ahi,P.power')';
    powerDelay = filter(Blo,Alo,powerHigh')'; % need to have this acting on the HP filtered power or you won't get negative delayed values
%     disp(num2str(tauHP));
    
    if isfield(P,'verbose')
        if P.verbose == 1
            % make a figure to show impulse responses
            x = [1;zeros(19,1)];
            fHigh = filter(Bhi,Ahi,x);
            fLow = filter(Blo,Alo,fHigh);
            figure; hold on;
            plot(P.timepoints(1:length(x)),fHigh,'k.-');
            plot(P.timepoints(1:length(x)),fLow,'r.-');
            legend('high pass','low pass');
            xlabel('time (s)');
            ylabel('filter amplitude');
        end
    end
    
    % compute HRC output and the four quadrants pp, mm, pm, and mp
    HRC = powerHigh(1:end-1,:).*powerDelay(2:end,:) - powerDelay(1:end-1,:).*powerHigh(2:end,:);
    HRCpp = rect(powerHigh(1:end-1,:)).*rect(powerDelay(2:end,:)) - ...
        rect(powerDelay(1:end-1,:)).*rect(powerHigh(2:end,:));
    HRCmm = rect(-powerHigh(1:end-1,:)).*rect(-powerDelay(2:end,:)) - ...
        rect(-powerDelay(1:end-1,:)).*rect(-powerHigh(2:end,:));
    % define: pm means p happened earlier, delay filtered signal is positive
    HRCpm = rect(-powerHigh(1:end-1,:)).*rect(powerDelay(2:end,:)) - ... 
        rect(powerDelay(1:end-1,:)).*rect(-powerHigh(2:end,:));
    HRCmp = rect(powerHigh(1:end-1,:)).*rect(-powerDelay(2:end,:)) - ...
        rect(-powerDelay(1:end-1,:)).*rect(powerHigh(2:end,:));

    % definitions above have inverted sign; want positive to be rising, so
    % invert here.
    P.estimators.HRC.HRC = -HRC;
    P.estimators.HRC.HRCpp = -HRCpp;
    P.estimators.HRC.HRCmm = -HRCmm;
    P.estimators.HRC.HRCpm = -HRCpm;
    P.estimators.HRC.HRCmp = -HRCmp;
    P.estimators.HRC.powerHigh = powerHigh;
    P.estimators.HRC.powerDelay = powerDelay;    

    function C = matchCorrNet(I,dy,dt,val1,val2)
        v1 = double(I==val1); % 0 or 1
        v2 = double(I==val2);
        % C can take on values of 1, 0, or -1, but only when v1 and v1 are
        % both 1 can it be non-zero; .* is an and-gate here
        C = v1(1:end-dy,1:end-dt).*v2(1+dy:end,1+dt:end) - ...
            v1(1+dy:end,1:end-dt).*v2(1:end-dy,1+dt:end);
        C = [zeros(size(C,1),1),C]; % add zeros to front to make it the same size as input I in dim 2
    end
end

function P = extractMotion(P)

    % extract 1d motion using 2d algorithms from Matlab
    % I is intensity field in x (rows) and time (columns)

    % create flow object; several options here
    opticFlow = opticalFlowHS('VelocityDifference',0.00001); % this seems to work as well as any
    % opticFlow = opticalFlowLK;
    % opticFlow = opticalFlowFarneback;

    I = P.power;
    
    for ii=1:size(I,2)
    
        grayIm = I(:,ii)*ones(1,30); % make a 2d gray scale for this thing; crashes when this is too narrow for LK method
        flow = estimateFlow(opticFlow,grayIm);
        uField(:,ii) = flow.Vy(:,2); % take the middle one
    
    end;

    P.uField = uField;
    P.uEst = mean(uField,1); % take mean over frequencies; but mean indicates that units don't really make sense.
    P.uUnits = log(P.freq(2)/P.freq(1))/log(2) / (P.timepoints(2)-P.timepoints(1)); % units in octaves/s
    P.uField = P.uField*P.uUnits; % put everything in these units already!
    P.uEst = P.uEst*P.uUnits;

end

function plotSpectrogram(P,varargin)

    if length(varargin)
        lastT = varargin{1}; % time to end plot in seconds
    end
    
    % if the velocity field has been computed, do this...
    if isfield(P,'uField')
        figure;
        subplot(3,1,[1 2]);
        imagesc([0,max(P.timepoints)], [],-P.power);
        ylabel('frequency (Hz)');
        xlabel('time (s)');
        set(gca,'xlim',[0 min(10,max(P.timepoints))],'ydir','normal');
        [tickLoc,tickVal] = findTicks(P.freq);
        set(gca,'ytick',tickLoc,'yticklabel',tickVal);
        box off;
        colormap('gray');

        % align this with the axes above
        subplot(3,1,3); 
        plot(P.timepoints,P.uEst*P.uUnits);
        hold on; plot([0 10],[0 0],'k:');
        ylabel('vel (oct/s, prob a.u.)');
        xlabel('time (s)')
        box off;
        set(gca,'xlim',[0 min(10,max(P.timepoints))],'ylim',[-1.1 1.1]*max(abs(P.uEst*P.uUnits)));
        drawnow;

    end

    % if there are regressors present, plot up the binary plot and plot up
    % the HRC plot...
    if isfield(P,'estimators')
        % first, binary estimators
        figure;
        subplot(3,1,[1 2]);
        imagesc([0,max(P.timepoints)], [],-P.powerBin);
        [tickLoc,tickVal] = findTicks(P.freq);
        set(gca,'xlim',[0 min(10,max(P.timepoints))],'ydir','normal');
        set(gca,'ytick',tickLoc,'yticklabel',tickVal)
        ylabel('freq (Hz)');
        colormap('gray');
        box off;

        subplot(3,1,3);
        pcolor([0 0 0;0 0 0; 0 1 1; 1 1 1 ; 1 1 1; 1 1 1]');
        colormap('gray');
        set(gca,'dataa',[1 1 1]);
        
        % make each figure in turn like the original, just to make life
        % easier
        figure;
        cmRB = makeRBcolormap(1);
        subplot(3,1,[1 2]);
        imagesc([0,max(P.timepoints)], [],P.estimators.binary.hhnet);
        set(gca,'xlim',[0 min(10,max(P.timepoints))],'ydir','normal');
        set(gca,'ytick',tickLoc,'yticklabel',tickVal);
        ylabel('freq (Hz)');
        title('hhnet');
        colormap(cmRB);
        box off;

        subplot(3,1,3);
        plot(P.timepoints,mean(P.estimators.binary.hhnet,1));
        hold on; plot([0 10],[0 0],'k:');
        set(gca,'xlim',[0 min(10,max(P.timepoints))],'ylim',[-1.1 1.1]*max(abs(mean(P.estimators.binary.hhnet,1))));
        ylabel('net ++');
        box off;
        
        figure;
        subplot(3,1,[1 2]);
        imagesc([0,max(P.timepoints)], [],P.estimators.binary.llnet);
        set(gca,'xlim',[0 min(10,max(P.timepoints))],'ydir','normal');
        set(gca,'ytick',tickLoc,'yticklabel',tickVal)
        ylabel('freq (Hz)');
        title('llnet');
        colormap(cmRB);
        box off;

        subplot(3,1,3);
        plot(P.timepoints,mean(P.estimators.binary.llnet,1));
        hold on; plot([0 10],[0 0],'k:');
        set(gca,'xlim',[0 min(10,max(P.timepoints))],'ylim',[-1.1 1.1]*max(abs(mean(P.estimators.binary.hhnet,1))));
        ylabel('net --');
        box off;

        figure;
        subplot(3,1,[1 2]);
        imagesc([0,max(P.timepoints)], [],P.estimators.binary.hlnet);
        set(gca,'xlim',[0 min(10,max(P.timepoints))],'ydir','normal');
        set(gca,'ytick',tickLoc,'yticklabel',tickVal)
        ylabel('freq (Hz)');
        title('hlnet');
        colormap(cmRB);
        box off;

        subplot(3,1,3);
        plot(P.timepoints,mean(P.estimators.binary.hlnet,1));
        hold on; plot([0 10],[0 0],'k:');
        set(gca,'xlim',[0 min(10,max(P.timepoints))],'ylim',[-1.1 1.1]*max(abs(mean(P.estimators.binary.hhnet,1))));
        ylabel('net +-');
        box off;

        figure;
        subplot(3,1,[1 2]);
        imagesc([0,max(P.timepoints)], [],P.estimators.binary.lhnet);
        set(gca,'xlim',[0 min(10,max(P.timepoints))],'ydir','normal');
        set(gca,'ytick',tickLoc,'yticklabel',tickVal)
        ylabel('freq (Hz)');
        title('lhnet');
        colormap(cmRB);
        box off;

        subplot(3,1,3);
        plot(P.timepoints,mean(P.estimators.binary.lhnet,1));
        hold on; plot([0 10],[0 0],'k:');
        set(gca,'xlim',[0 min(10,max(P.timepoints))],'ylim',[-1.1 1.1]*max(abs(mean(P.estimators.binary.hhnet,1))));
        ylabel('net -+');
        xlabel('time (s)');
        box off;

        figure;
        subplot(3,1,[1 2]);
        imagesc([0,max(P.timepoints)], [],P.estimators.binary.fullnet);
        set(gca,'xlim',[0 min(10,max(P.timepoints))],'ydir','normal');
        set(gca,'ytick',tickLoc,'yticklabel',tickVal)
        ylabel('freq (Hz)');
        title('fullnet');
        colormap(cmRB);
        box off;

        subplot(3,1,3);
        plot(P.timepoints,mean(P.estimators.binary.fullnet,1));
        hold on; plot([0 10],[0 0],'k:');
        set(gca,'xlim',[0 min(10,max(P.timepoints))],'ylim',[-1.1 1.1]*max(abs(mean(P.estimators.binary.fullnet,1))));
        ylabel('full HRC');
        xlabel('time (s)');
        box off;
 
        % second, HRC estimators, but we're not really planning to show
        % these. if we do, could copy figure code from above and show the
        % false color HRC plots
        figure;
        subplot(6,1,[1 2]);
        imagesc([0,max(P.timepoints)], [],-P.power);
        [tickLoc,tickVal] = findTicks(P.freq);
        set(gca,'xlim',[0 min(10,max(P.timepoints))],'ydir','normal');
        set(gca,'ytick',tickLoc,'yticklabel',tickVal)
        ylabel('freq (Hz)');
        colormap('gray');
        box off;

        subplot(6,1,3);
        plot(P.timepoints,mean(P.estimators.HRC.HRCpp,1));
        hold on; plot([0 10],[0 0],'k:');
        set(gca,'xlim',[0 min(10,max(P.timepoints))]);
        ylabel('net ++');
        box off;

        subplot(6,1,4);
        plot(P.timepoints,mean(P.estimators.HRC.HRCmm,1));
        hold on; plot([0 10],[0 0],'k:');
        set(gca,'xlim',[0 min(10,max(P.timepoints))]);
        ylabel('net --');
        box off;

        subplot(6,1,5);
        plot(P.timepoints,mean(P.estimators.HRC.HRCpm,1));
        hold on; plot([0 10],[0 0],'k:');
        set(gca,'xlim',[0 min(10,max(P.timepoints))]);
        ylabel('net +-');
        box off;

        subplot(6,1,6);
        plot(P.timepoints,mean(P.estimators.HRC.HRCmp,1));
        hold on; plot([0 10],[0 0],'k:');
        set(gca,'xlim',[0 min(10,max(P.timepoints))]);
        ylabel('net -+');
        xlabel('time (s)');
        box off;

    end


    function [ind,f] = findTicks(freq)
        Ffind = [100,200,400,800,1600,3200];
        for ii=1:length(Ffind)
            ind(ii) = find(freq == Ffind(ii));
            f(ii) = Ffind(ii);
        end
    end

    function cm = makeRBcolormap(varargin)

        if length(varargin)==0
            lowfrac = 0.5;
        else
            lowfrac=1-varargin{1};
        end
        A=lowfrac*100;
        dA=(100-A)/63;
        
        cm = [[ones(1,64)*100,[100:-dA:A]]',[A:dA:100,100:-dA:A]',[A:dA:100,100*ones(1,64)]'];
        cm = cm/100;
        cm = cm(end:-1:1,:);

    end

end

function [P,varargout] = extractSpectrogram(wav,Fs)

    % parameters here
    n = 40; % samples per second
    overlap = 0; % overlap between subsequent samples; 0.5 or 0.75 or 0?
    nTones = 20; % number of tones per octave; do want this to be matching? 45 was first round that looked good.
    scale = 100*2.^[0:1/nTones:6-1/nTones]; % in Hz; roughly matches what we did in pscyhophysics
    
    [spec,freq,timepoints] = spectrogram(wav,Fs/n,Fs/n*overlap,scale,Fs); % window is Fs/n = 1/n second; 0 overlap; 
    
    % new struture to use on everything
    P.power = abs(spec); % power is amplitude squared, but we want amplitude here, slightly compressed dynamic range and easier to deal with...
    P.freq =  freq;
    P.timepoints = timepoints;

end

function a = rect(x)

    a = max(x,0);

end

function saveOpenFigures(foldName)

    figHandles = findall(0,'Type','figure');

    for ii=1:length(figHandles)
        saveas(figHandles(ii),[foldName '/fig' num2str(figHandles(ii).Number)],'fig');
        saveas(figHandles(ii),[foldName '/fig' num2str(figHandles(ii).Number)],'pdf');

    end

end
