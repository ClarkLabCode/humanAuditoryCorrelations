%% Script to create various tone types from psychophysical experiments

PLOTENV = 1;

% Which stimulus type is being played
corrType = 'ternScint'; % 'ternScint', 'ternRandom', 'glider3conv', 'glider3div', 'corrPips'

% ternScint = ternary pattern, short-range correlations -- tones are on,
% off, or intermediate

% ternRandom = random ternary pattern, control

% glider3conv = converging gliders with center of mass moving to higher
% tones by default (as in glider, which is same with 2pt)

% glider3div = diverging glider with center of mass moving to higher tones
% by default (as in glider, which is same with 2pt)

% corrPips = Poisson timed pips that correlate with P to higher or lower
% tones at some delay; for testing, additional parameters are contained in
% this case set of code

% Correlation and direction
P = 1; % correlation sign = +1 or -1
DIRECTION = 1; % +/-1, +1 is up, -1 is down for correlation displacement in tone space

%%%%%% basic parameters below
Fs = 2e4; % freq = 44.1e3 ideally, but 20e3 should mostly be fine, since it's above the Nyquist frequency of our highest frequency of ~6 kHz
T = 5; % in seconds

% Number of pure tones per octave (spaced equally in log-space)
nTones = 15;

% Set scale
scale = 100*2.^[1:1/nTones:6-1/nTones];

% Duration of the notes
noteT = 1/6; 

deltaNote = 1; % how many notes to move with correlation in frequency dimension
deltaT = 1; % how many notes to move with corelation in time dimension

corrPipTypeNumber = 0; % in range 0-4, with effects shown below
switch corrPipTypeNumber
    case 0
        corrPipType = 0; % = 0 to get all pip types, [+1, +1] to get ++ types correlated only, etc. 1/2 of pips are in a correlated pair; rest are uncorrelated but balanced in louder/softer
    case 1
        corrPipType = [1 1];
    case 2
        corrPipType = [1 -1];
    case 3
        corrPipType = [-1 1];
    case 4
        corrPipType = [-1 -1];
end


%%% start creating the envelopes = playNote variables
eta = randn(length(scale),T/noteT + 2); % use the size of the whole scale

% Stimulus type options
switch corrType
        
    case 'ternScint'
        eta = double(eta>0)*2 - 1; % binarize this
        playNote = ( (eta + P*circshift(eta,[deltaNote deltaT]))+2 ) * 1/4; % in range 0,1
    
    case 'ternRandom'
        eta2 = randn(size(eta));
        etaB = double(eta > 0);
        eta2B = double(eta2 > 0);
        playNote = etaB+eta2B;
        playNote = playNote/2; % in range 0,1
        
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
                    playNote(jj,tt) = P*playNote(jj,tt-1)*playNote(end,tt-1); % put on other side
                end
            end
        end
        playNote = (playNote+1)/2;
        
    case 'glider3div'
        % make a diverging 3 point glider, with c.o.m. going to higher
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
        noteT = 1/Fs; % to make it compatible with toneEnv coming later
        playNote = zeros(size(eta,1),length([0:T*Fs])); % total number of samples
        pipDur = 1/20; % duration of each pip in seconds
        pipDeltaT = 0.04; % interval between pip starts in seconds

        if length(corrPipType) == 1
            % delta functions in here
            playNote = double(rand(size(playNote))<(dens*1/Fs)); % creates Poisson distribution of events over all notes
            ff = find(playNote); 
            playNote(ff) = (double(randn(size(ff))<0)*2 - 1); % sets those events at random to +/- 1 
            playNote = playNote + P*circshift(playNote,[deltaNote,round(Fs*pipDeltaT)]); % create correlated pips
        elseif length(corrPipType) == 2
            Ploc = prod(corrPipType); % local parity
            playNote1 = double(rand(size(playNote))<(dens*1/Fs/2)); % creates Poisson distribution of events over all notes
            playNote2 = double(rand(size(playNote))<(dens*1/Fs/2)); % non-modulated type; note there end up being half as many pips here as one would want to match the correlated pairs.
            playNote1 = playNote1*corrPipType(1); % set these up, equal to +/-1
            playNote1 = playNote1 + Ploc*circshift(playNote1,[deltaNote,round(Fs*pipDeltaT)]); % set up only that correlation
            if Ploc == 1 % same parity correlations, fill in random pips with opposite sign
                playNote2 = playNote2*corrPipType(1)*-1; % non-correlated type has opposing contrast if positive parity
            else % fill in random pips with random signs
                ff2 = find(playNote2);
                playNote2(ff2) = (double(randn(size(ff2))<0)*2 - 1); % sets those events at random to +/- 1 
            end
            playNote = playNote1 + playNote2; % add together the correlated and uncorrelated sets here.
        else
            disp('corrPiptype error')
            return;
        end
        playNote = filter(ones(1,round(pipDur*Fs)),1,playNote')';
        playNote = min(playNote,1); % don't let it get of range of [-1 1]
        playNote = max(playNote,-1);
        playNote = 1/2*(playNote+1); % shift to [0 1];
        
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

% Make the envelope
switch corrType
    case 'corrPips'
        toneEnv = playNote;
    otherwise
        toneEnv = interp1([1:1:size(eta,2)]*noteT,playNote',[0:T*Fs]/Fs+noteT,'nearest')'; % interpolates to get blocky windows
end
toneEnv = filtfilt(exp(-[0:1000]/10)/sum(exp(-[0:1000]/10)),1,toneEnv')'; % time of filter in samples; 
% this is 10 samples. at 44kHz, this corresponds to a timescale of 0.2 ms,
% it's about 0.5 ms for 20kHz sampling.

% Make the tones
t = [0:T*Fs]/Fs; % time in seconds
clear rawTones;
for ii=1:length(scale)
    rawTones(ii,:)=sin(t*scale(ii)*2*pi); % can also add in random phase for each frequency here, with no perceptual effect
end

if PLOTENV
    figure;
    imagesc(t,scale,toneEnv);
    set(gca,'ydir','normal');
    xlabel('time (s)');
    ylabel('log frequency (Hz)');
    set(gca,'ytick',[scale(1),scale(end)]);
    title('envelope image');
    colormap('gray');
    colorbar;
end


% Scale the tones to have equal perceptual loudness -- this isn't really
% necessary for any of the percepts, we believe
[spl,freq]=iso226(60); % iso226 standard
spli = interp1(freq,spl,scale,'linear'); % interpolate the amplitudes
amplitude = 10.^(spli/20); % relative amplitudes


fullRec = amplitude(:)'*(rawTones.*toneEnv);
soundsc(fullRec,Fs);
