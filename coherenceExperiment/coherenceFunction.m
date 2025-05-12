function [outStruct] = coherenceFunction(params)

%% Variables set by params

% Stimulus duration
T = params.stimDur; % in seconds

VERBOSE = params.verbose; % if 1, then plot things up for dubugging

% Number of pure tones per octave (spaced equally in log-space)
nTones = params.nTonesPerOctave;

% Duration of the notes
noteT = params.noteDur; 

deltaNote = params.deltaNote; % how many notes to move with correlation in frequency dimension
deltaT = params.deltaT; % how many notes to move with corelation in time dimension

coh = params.coh; % coherence [0,1], with 1 being strongest signal (same as no coherence parameter used) and 0 being same as random

% Which type of stimulus is played
corrType = params.corrType;

% Correlation and direction
P = params.corrParity; % correlation sign = +1 or -1
DIRECTION = params.displacement; % +/-1, +1 is up, -1 is down for correlation displacement in tone space

% Cosine envelope
COSINE_ENV = 0; % 0 or 1 -- whether to add a cosine fall off on frequencies, to kill those near top and bottom of range

%% The following don't change

Fs = 2e4; %20 kHz sampling

% Set scale
scale = 100*2.^[1:1/nTones:6-1/nTones];

eta = randn(nTones,T/noteT+2); % size of the whole thing

switch corrType

    case 'ternRandom'
        eta2 = randn(nTones,T/noteT+2);
        etaB = double(eta > 0);
        eta2B = double(eta2 > 0);
        playNote = etaB+eta2B;
        playNote = playNote/2;
        
    case 'ternScint'
        eta = double(eta>0)*2 - 1;
        playNote = ( (eta + P*circshift(eta,[deltaNote deltaT]))+2 ) * 1/4; % creates positive and negative windows

    case 'ternScintCoh'
        eta = double(eta>0)*2 - 1;
        eta2 = P*circshift(eta,[deltaNote deltaT]);
        randinds = find(rand(size(eta2)) > coh);
        if length(randinds)
            eta2(randinds) = (rand(size(randinds)) > 0.5)*2 - 1;
        end
        playNote = ( eta + eta2 +2 ) * 1/4; % creates positive and negative windows

    otherwise
        disp('error: none chosen');
        return;
end

% Choose direction
if DIRECTION == 1

    % Do nothing
elseif DIRECTION == -1
    playNote = playNote(end:-1:1,:);
else
    disp('you messed something up with the DIRECTION variable');
    return;
end

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
toneEnv = interp1([1:1:size(eta,2)]*noteT,playNote',[0:T*Fs]/Fs+noteT,'nearest')'; % interpolates to get blocky windows
[Bf,Af]=butter(1,1/pi/500,'low');
toneEnv = filtfilt(Bf,Af,toneEnv')';
toneEnv = repmat(toneEnv,[5 1]); % make it the same for all octaves

if VERBOSE
    figure;
    imagesc([1 T],[min(scale) max(scale)],toneEnv);
    set(gca,'ydir','normal');
    xlabel('time (s)');
    ylabel('pitch (Hz, log scale)');
    set(gca,'ytick',[min(scale) max(scale)]);
    title('amplitude pitches over time');
    colorbar;
end

% Make the tones
t = [0:T*Fs]/Fs; % time in seconds
clear rawTones;
for ii=1:length(scale)
    rawTones(ii,:)=sin(t*scale(ii)*2*pi);
end

if VERBOSE
    figure;
    imagesc([1 T],[min(scale) max(scale)],rawTones.*toneEnv);
    set(gca,'ydir','normal');
    xlabel('time (s)');
    ylabel('pitch (Hz, log scale)');
    title('amplitude pitches over time');
    set(gca,'ytick',[min(scale) max(scale)]);
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
