function [outStruct] = pipDeltaNoteFunction(params)

%% Variables set by params

Fs = 2e4; % freq

% Stimulus duration
T = params.stimDur; % in seconds

VERBOSE = params.verbose; % if 1, then plot things up for dubugging

% Number of pure tones per octave (spaced equally in log-space)
nTones = params.nTonesPerOctave; 

% Set scale
scale = 100*2.^[1:1/nTones:6-1/nTones];

% Duration of the notes
noteT = params.noteDur;

deltaNote = params.deltaNote; % how many notes to move with correlation in frequency dimension
deltaT = params.deltaT; % how many notes to move with corelation in time dimension

% Which stimulus type is presented
corrType = params.corrType; % 'corrPips' for this experiment

% Time between pips
pipDeltaT = params.pipDeltaT;

% Duration of each pip
pipDur = params.pipDur;

% Correlation (direction is simply set by positive or negative pip delta note)
P = params.corrParity; % correlation sign = +1 or -1

% Cosine envelope
COSINE_ENV = 0; % 0 or 1 -- whether to add a cosine fall off on frequencies, to kill those near top and bottom of range

% Size of the whole thing
eta = randn(nTones,T/noteT+2);

% Stimulus type options
switch corrType
        
    case 'ternScint'
        eta = double(eta>0)*2 - 1;
        playNote = ( (eta + P*circshift(eta,[deltaNote deltaT]))+2 ) * 1/4; % creates positive and negative windows
  
    case 'glider3conv'
        % make a converging 3 point glider, with c.o.m. going to higher
        % tones
        eta = double(eta>0)*2 - 1;
        playNote = eta;
        for tt = (1 + deltaT):size(eta,2)
            for jj = 1:size(eta,1)
                if jj>1
                    playNote(jj,tt) = P*playNote(jj,tt-1)*playNote(jj-1,tt-1);
                else
                    playNote(jj,tt) = P*playNote(jj,tt-1)*playNote(end,tt-1);
                end
            end
        end
        playNote = (playNote+1)/2;
        
    case 'glider3div'
        % make a divverging 3 point glider, with c.o.m. going to higher
        % tones
        eta = double(eta>0)*2 - 1;
        playNote = eta;
        for tt = (1 + deltaT):size(eta,2)
            for jj = 1:size(eta,1)
                if jj>1
                    playNote(jj,tt) = P*playNote(jj,tt-1)*playNote(jj-1,tt-1);
                else
                    playNote(jj,tt) = P*playNote(jj,tt-1)*playNote(end,tt-1);
                end
            end
        end
        playNote = playNote(end:-1:1,end:-1:1);
        playNote = (playNote+1)/2;
        
    case 'corrPips'
        dens = params.dens; % Hz, density of pips per tone, per second
        noteT = 1/Fs; % to make it compatible with toneEnv coming later
        playNote = zeros(size(eta,1),length([0:T*Fs])); % total number of samples
        playNote = double(rand(size(playNote))<(dens*1/Fs)); % creates Poisson distribution of events over all notes
        ff = find(playNote); 
        playNote(ff) = (double(randn(size(ff))<0)*2 - 1); % sets those events at random to +/- 1 
        playNote = playNote + P*circshift(playNote,[deltaNote,round(Fs*pipDeltaT)]); % create correlated pips
        playNote = filter(ones(1,round(pipDur*Fs)),1,playNote')';
        playNote = min(playNote,1); % don't let it get of range of [-1 1]
        playNote = max(playNote,-1);
        playNote = 1/2*(playNote+1); % shift to [0 1];
        
    otherwise
        disp('error: none chosen');
        return;
end

% First figure
if VERBOSE
    figure; 
    imagesc(playNote);
    xlabel('time number');
    ylabel('tone number');
    set(gca,'ydir','normal');
    title('amplitude of pitches over time');
    colorbar;
end

% Make the envelope
switch corrType
    case 'corrPips'
        toneEnv = playNote;
    otherwise
        toneEnv = interp1([1:1:size(eta,2)]*noteT,playNote',[0:T*Fs]/Fs+noteT,'nearest')'; % interpolates to get blocky windows
end
toneEnv = filtfilt(exp(-[0:1000]/10),1,toneEnv')'; % time of filter in samples
toneEnv = repmat(toneEnv,[5 1]); % make it the same for all octaves

% Another figure
if VERBOSE
figure;
imagesc([1 T],[min(scale) max(scale)],toneEnv);
set(gca,'ydir','normal');
xlabel('time (s)');
ylabel('pitch (Hz, log scale)');
title('amplitude pitches over time');
colorbar;
end

% Make the tones
t = [0:T*Fs]/Fs; % time in seconds
clear rawTones;
for ii=1:length(scale)
    rawTones(ii,:)=sin(t*scale(ii)*2*pi);
end

% Another figure
if VERBOSE
figure;
imagesc([1 T],[min(scale) max(scale)],rawTones.*toneEnv);
set(gca,'ydir','normal');
xlabel('time (s)');
ylabel('pitch (Hz, log scale)');
title('amplitude pitches over time');
colorbar;
end


% Scale the tones to have equal perceptual loudness
[spl,freq]=iso226(60); % iso226 standard
spli = interp1(freq,spl,scale,'linear'); % interpolate the amplitudes
amplitude = 10.^(spli/20); % relative amplitudes
if COSINE_ENV
    env = cos([0:(length(scale)-1)]/(length(scale)-1)*pi);
    amplitude = amplitude.*env;
end

fullRec = amplitude(:)'*(rawTones.*toneEnv);

if VERBOSE
    soundsc(fullRec,Fs);
end

outStruct.params = params;
outStruct.waveForm = fullRec;
outStruct.Fs = Fs;
outStruct.playNote = playNote;
