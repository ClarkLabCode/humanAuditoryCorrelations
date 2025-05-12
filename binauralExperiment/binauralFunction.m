function [outStruct] = binauralFunction(params)

%% Variables set by params

outStruct.params = params;

Fs = 2e4; % Hz sampling

T = params.stimDur; % in seconds

REPEATOCTAVE = 0; % new

% Number of pure tones per octave (spaced equally in log-space)
nTones = params.nTonesPerOctave;

% Set scale
scale = 100*2.^[1:1/nTones:6-1/nTones];

% Duration of the notes
noteT = params.noteDur;

if REPEATOCTAVE
    eta = randn(nTones,T/noteT+2); % size of this whole thing for one octave
else
    eta = randn(length(scale),T/noteT + 2); % use the size of the whole scale
end

deltaNote = params.deltaNote; % how many notes to move with correlation in frequency dimension
deltaT = params.deltaT; % how many notes to move with corelation in time dimension

% Which type of stimulus is played
corrType = params.corrType;

% Set correlation and direction
P = params.corrParity; % correlation sign = +1 or -1
DIRECTION = params.displacement; % +/-1, +1 is up, -1 is down for correlation displacement in tone space
COSINE_ENV = 0; % 0 or 1 -- whether to add a cosine fall off on frequencies, to kill those near top and bottom of range

% New
corrPipType = [1 -1]; % = 0 to get all pip types, [+1, +1] to get ++ types correlated only, etc. 1/2 of pips are in a correlated pair; rest are uncorrelated but balanced in louder/softer

switch corrType
        
    case 'ternScint'
        eta = double(eta>0)*2 - 1;
        playNote = ( (eta + P*circshift(eta,[deltaNote deltaT])) + 2 ) * 1/4; % creates positive and negative windows

    case 'binauralScint'
        eta = double(eta>0); % binary case
        eta2 = (((circshift(eta,[deltaNote deltaT])* 2) - 1) * P + 1 ) / 2;
        playNoteLeft = eta; % each ear
        playNoteRight = eta2;

    case 'binauralScintRand'
        eta = double(eta>0); % binary case
        eta2 = double(randn(size(eta))>0);
        playNoteLeft = eta;
        playNoteRight = eta2;
        
    case 'ternRand'
        eta = double(eta>0)*2 - 1;
        eta2 = double(randn(size(eta))>0)*2 - 1;
        playNote = ( eta + eta2 + 2 ) * 1/4; % results in random ternary matrix with values (0,0.5,1) with prop (0.25, 0.5, 0.25)
        
    case 'glider3conv'

        % Make a converging 3 point glider, with c.o.m. going to higher
        % tones
        eta = double(eta>0)*2 - 1;
        playNote = eta;
        for tt = (1 + deltaT):size(eta,2)
            for jj = 1:size(eta,1)
                if jj>1
                    playNote(jj,tt) = P*playNote(jj,tt-1)*playNote(jj-1,tt-1);
                else
                    playNote(jj,tt) = P*playNote(jj,tt-1)*playNote(end,tt-1); % put on other side
                end
            end
        end
        playNote = (playNote+1)/2;
        
    case 'glider3div'

        % Make a diverging 3 point glider, with c.o.m. going to higher
        % tones
        eta = double(eta>0)*2 - 1;
        playNote = eta;
        for tt = (1 + deltaT):size(eta,2)
            for jj = 1:size(eta,1)
                if jj>1
                    playNote(jj,tt) = P*playNote(jj,tt-1)*playNote(jj-1,tt-1);
                else
                    playNote(jj,tt) = P*playNote(jj,tt-1)*playNote(end,tt-1); % put on other side
                end
            end
        end
        playNote = playNote(end:-1:1,end:-1:1);
        playNote = (playNote+1)/2;
        
    case 'corrPips'
        dens = 4; % Hz, density of pips per tone, per second, for initializing
        noteT = 1/Fs; % to make this compatible with toneEnv coming later
        playNote = zeros(size(eta,1),length([0:T*Fs])); % total number of samples
        pipDur = .04; % duration of each pip in seconds
        pipDeltaT = .04; % interval between pip starts in seconds

        if length(corrPipType) == 1
            % delta functions in here
            playNote = double(rand(size(playNote))<(dens*1/Fs)); % creates Poisson distribution of events over all notes
            ff = find(playNote); 
            playNote(ff) = (double(randn(size(ff))<0)*2 - 1); % sets those events at random to +/- 1 
            playNote = playNote + P*circshift(playNote,[deltaNote,round(Fs*pipDeltaT)]); % create correlated pips
        elseif length(corrPipType) == 2
            Ploc = prod(corrPipType); % local parity
            playNote1 = double(rand(size(playNote))<(dens*1/Fs/2)); % creates Poisson distribution of events over all notes
            playNote2 = double(rand(size(playNote))<(dens*1/Fs/2)); % non-modulated type
            playNote1 = playNote1*corrPipType(1);
            playNote1 = playNote1 + Ploc*circshift(playNote1,[deltaNote,round(Fs*pipDeltaT)]); % set up only that correlation
            if Ploc == 1
                playNote2 = playNote2*corrPipType(1)*-1; % non-correlated type has opposing contrast if positive parity
            else
                ff2 = find(playNote2);
                playNote2(ff2) = (double(randn(size(ff2))<0)*2 - 1); % sets those events at random to +/- 1 
            end
            playNote = playNote1 + playNote2; % add together the correlated and uncorrelated sets here
        else
            disp('corrPiptype error')
            return;
        end
        playNote = filter(ones(1,round(pipDur*Fs)),1,playNote')';
        playNote = min(playNote,1); % don't let it get of range of [-1 1]
        playNote = max(playNote,-1);
        playNote = 1/2*(playNote+1); % shift to [0 1];

    case 'posPips'
        dens = 4; % Hz, density of pips per tone, per second, for initializing
        noteT = 1/Fs; % to make this compatible with toneEnv coming later
        playNote = zeros(size(eta,1),length([0:T*Fs])); % total number of samples
        pipDur = .04; % duration of each pip in seconds
        pipDeltaT = .04; % interval between pip starts in seconds
        Ploc = 1; % local parity = 1
        playNote1 = double(rand(size(playNote))<(dens*1/Fs/2)); % creates Poisson distribution of events over all notes
        playNote1 = playNote1 + Ploc*circshift(playNote1,[deltaNote,round(Fs*pipDeltaT)]); % set up only that correlation
        playNote1 = filter(ones(1,round(pipDur*Fs)),1,playNote1')'; % make the windows the right length
        playNote = playNote1 * 1/2;

    case 'gridMoving'
        playNote = eye(size(eta,1));
        playNote = repmat(playNote,[1,ceil(size(eta,2)/size(eta,1))]);
        playNote = playNote(:,[1:size(eta,2)]);

    case 'gridOpp'
        playNote = eye(size(eta,1));
        playNote = repmat(playNote,[1,ceil(size(eta,2)/size(eta,1))]);
        playNote = playNote(:,[1:size(eta,2)]);
        playNote = playNote + playNote(end:-1:1,:);

    otherwise
        disp('error: none chosen');
        return;
end

% Choose direction
if DIRECTION == 1

    % Do nothing
elseif DIRECTION == -1
        playNoteLeft = playNoteLeft(end:-1:1,:);
        playNoteRight = playNoteRight(end:-1:1,:);
else
    disp('you messed something up with the DIRECTION variable');
    return;
end

% Make the envelope
switch corrType
    case 'corrPips'
        toneEnv = playNote;
    case 'posPips'
        toneEnv = playNote;
    otherwise
        toneEnvLeft = interp1([1:1:size(eta,2)]*noteT,playNoteLeft',[0:T*Fs]/Fs+noteT,'previous')'; % interpolates to get blocky windows
        toneEnvRight = interp1([1:1:size(eta,2)]*noteT,playNoteRight',[0:T*Fs]/Fs+noteT,'previous')'; % interpolates to get blocky windows
end
toneEnvLeft = filtfilt(exp(-[0:1000]/10),1,toneEnvLeft')'; % time of filter in samples
toneEnvRight = filtfilt(exp(-[0:1000]/10),1,toneEnvRight')'; % time of filter in samples

if REPEATOCTAVE
    toneEnvLeft = repmat(toneEnvLeft,[5 1]); % make it the same for all octaves
    toneEnvRight = repmat(toneEnvRight,[5 1]); % make it the same for all octaves
end

% Make the tones
t = [0:T*Fs]/Fs; % time in seconds
clear rawTones;
for ii=1:length(scale)
    rawTones(ii,:)=sin(t*scale(ii)*2*pi);
end

% Scale the tones to have equal perceptual loudness
[spl,freq]=iso226(60); % iso226 standard
spli = interp1(freq,spl,scale,'linear'); % interpolate the amplitudes
amplitude = 10.^(spli/20); % relative amplitudes
if COSINE_ENV
    env = cos([0:(length(scale)-1)]/(length(scale)-1)*pi);
    amplitude = amplitude.*env;
end

fullRecLeft = amplitude(:)'*(rawTones.*toneEnvLeft);
fullRecRight = amplitude(:)'*(rawTones.*toneEnvRight);

outStruct.waveFormLeft = fullRecLeft;
outStruct.waveFormRight = fullRecRight;
outStruct.Fs = Fs;
outStruct.playNoteLeft = playNoteLeft;
outStruct.playNoteLeft = playNoteLeft;
