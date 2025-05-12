%% Start of script to run auditory psychophysical experiments (sweep coherence)

clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%
%% BEGIN PREP CODE %%
%%%%%%%%%%%%%%%%%%%%%

rng('shuffle'); % set seed for rand

%% First, let's get the subject's ID and store it

% Addpath data
subInfo = input('Enter the subject ID #:', 's');
subInfo = strcat(subInfo); % concatenate string entered
while exist(['DataCoherence/',subInfo,'.mat'],'file') % make sure no overwriting is happening
    subInfo = input('Subject # already exists! Enter new #:', 's'); % re-do if already used
    subInfo = strcat(subInfo);
end
mkdir(['stimuli/',subInfo]);

%% Second, we create the stimlus sequence

main_trial_types = 4; % how many different "types" of trials of interest

mainN = 60; % number of trials for EACH "main" trial type in experiment

coherence_vector = [0 0.125 0.25 0.5 0.75 1]; % vector of coherence values

tmpSeq = []; % initialize empty vector

% Build sequence vector
tmpSeq = [ones(mainN,1); 2*ones(mainN,1); 3*ones(mainN,1); 4*ones(mainN,1)];
tmp = rectpulse(coherence_vector,mainN/length(coherence_vector))';
tmpCoherenceSeq = repmat(tmp',1,4)';
randSeq = randperm(length(tmpSeq)); % get a randomization seed

% Randomize both sequences in the same way to maintain trial proporitions
trialTypeSeq = tmpSeq(randSeq);
CoherenceSeq = tmpCoherenceSeq(randSeq);

nTrials = length(trialTypeSeq); % total task N

%% Now that we have our shuffled sequence, let's pre-save the actual stimuli as WAVs

% Set fixed stimulus parameters first (these are key psychophysical variables)
params.stimDur = 1; % duration of stimulus in seconds (prevent ceiling effects)
params.verbose = 0; % include plotting (always keep OFF when using Psychtoolbox)
params.nTonesPerOctave = 15; % number of notes in each octave
params.noteDur = 1/6; % note duration (in seconds)
params.deltaT = 1; % delta time
params.deltaNote = 1; % delta note
Fs = 2e4; % freq

%% Now we use a for-loop to generate the tones using our stimulus sequence

for n = 1:nTrials
    
    % Vary coherence
    params.coh = CoherenceSeq(n);
    
    if trialTypeSeq(n) == 1 % TRIAL TYPE 1
        params.corrType = 'ternScintCoh';
        params.corrParity = 1; % positive correlation
        params.displacement = 1; % upward direction
        
    elseif trialTypeSeq(n) == 2 % TRIAL TYPE 2
        params.corrType = 'ternScintCoh';
        params.corrParity = 1; % positive correlation
        params.displacement = -1; % downward direction

    elseif trialTypeSeq(n) == 3 % TRIAL TYPE 3
        params.corrType = 'ternScintCoh';
        params.corrParity = -1; % negative correlation
        params.displacement = 1; % upward direction

    elseif trialTypeSeq(n) == 4 % TRIAL TYPE 4
        params.corrType = 'ternScintCoh';
        params.corrParity = -1; % negative correlation
        params.displacement = -1; % downward direction
           
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Call tone generating function %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Generate trial n's stimulus
    toneStim(n) = coherenceFunction(params);
    
    % Set WAV filename
    soundFilename{n} = ['stimuli/',subInfo,'/',subInfo,'_tone',num2str(n),'.wav'];

    % Scale waveform
    Y = toneStim(n).waveForm;
    Y = rescale(Y,-1,1); % rescale to prevent distorted noises

    % Write out WAV
    audiowrite(soundFilename{n},Y,Fs);
    
end

%%%%%%%%%%%%%%%%%%%%
%% MAIN TASK LOOP %%
%%%%%%%%%%%%%%%%%%%%

try
    %% Step 1

    % Sync bug
    Screen('Preference', 'SkipSyncTests', 1);

    % For psychportaudio (PPA) functioning
    AssertOpenGL;

    % Perform basic initialization of the sound driver
    InitializePsychSound;
    
    %% Step 2
    % Prepare "Screen"

    % Call default settings for setting up Psychtoolbox
    PsychDefaultSetup(2);
    
    % Get the screen numbers
    screens = Screen('Screens');

    % Select the external screen if it is present, else revert to the native
    screenNumber = max(screens);

    % Define black, white and grey
    black = BlackIndex(screenNumber);
    white = WhiteIndex(screenNumber);
    grey = white / 2;
    
    % Open an on screen window and color it grey
    [window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey);
    
    % Set the blend function for the screen
    Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

    % Get the size of the on screen window in pixels
    [screenXpixels, screenYpixels] = Screen('WindowSize', window);

    % Query the frame duration
    ifi = Screen('GetFlipInterval', window);
    
    % Get the centre coordinate of the window in pixels
    [xCenter, yCenter] = RectCenter(windowRect);

    % Set the text size
    Screen('TextSize', window, 22);
    
    %% Initialize fixation cross details

    % Screen Y fraction for fixation cross
    crossFrac = 0.0167;

    % Set the size of the arms for fixation cross
    fixCrossDimPix = windowRect(4) * crossFrac;

    % Set the coordinates for fixation cross
    xCoords = [-fixCrossDimPix fixCrossDimPix 0 0];
    yCoords = [0 0 -fixCrossDimPix fixCrossDimPix];
    allCoords = [xCoords; yCoords];

    % Set the line width for fixation cross
    lineWidthPix = 4;
    
    %% Initialize data frame for storing results

    % Initialize simple data struct
    data = struct('rt',nan(nTrials,1),'choice',nan(nTrials,1),'rt_violation',nan(nTrials,1),...
        'tone_start_t',nan(nTrials,1),'trial_start_t',nan(nTrials,1),'experiment_time',nan(nTrials,1),...
        'corrparity',nan(nTrials,1),'displacement',nan(nTrials,1),'coherence',nan(nTrials,1));
    
    % Define available keys to press
    upKey = KbName('UpArrow');
    downKey = KbName('DownArrow');
    RestrictKeysForKbCheck([upKey downKey]);
    
    %% Initialize events and fixed timing params

    e_play_tone = 1; % hold and play tone event code
    e_response = 0; % wait for response event code
    e_log = 0; % log response event code
    rt_max = 3 + params.stimDur(1); % maximum reaction time in seconds (from the beginning of the tone, so tone length is added)
    rt_min = 0.150 + params.stimDur(1); % minimum reaction time in seconds (from the beginning of the tone, so tone length is added)

    % Inter-trial interval
    iti = 0.750; % ITI (in seconds)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Initialize trial loop %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    trialnum = 1; % set first trial
    hold_start_t = GetSecs(); % log start of trial 1
    
    %% Step 3

    % Start the main loop
    while trialnum <= nTrials
        
        %% Trial has started
        
        if e_play_tone == 1 % enter the "play_tone" phase

            % show a fixation cross
            Screen('DrawLines', window, allCoords,...
                lineWidthPix, white, [xCenter yCenter], 2);
            Screen('Flip',window);
            
            %% Load sound for this trial %%

            % Specify filename of sound we want to load
            wavfilename = soundFilename{trialnum}; % load sound associated with this trial, from subject's directory in "stimuli"

            % Read WAV file from filesystem
            [y, freq] = psychwavread(wavfilename);
            wavedata = y';
            nrchannels = size(wavedata,1); % number of rows == number of channels
            
            % Make sure we have always two channels stereo output (PPA functioning issue)
            if nrchannels < 2
                wavedata = [wavedata ; wavedata];
                nrchannels = 2;
            end
            
            %% Buffer and play the sound %%
            
            % Return a handle to the audio device
            pahandle = PsychPortAudio('Open', [], [], 0, freq, nrchannels);
            
            % Fill the audio playback buffer with the audio data 'wavedata'
            PsychPortAudio('FillBuffer', pahandle, wavedata);
            repetitions = 1;

            % Start audio playback for 'repetitions' repetitions of the sound data,
            % start it immediately (0) and wait for the playback to start, return onset timestamp
            tone_start_t = PsychPortAudio('Start', pahandle, repetitions, 0, 1);
            
            % Toggle events
            e_play_tone = 0;
            e_response = 1;
            
            % Hold for stimulus playback
            WaitSecs(params.stimDur(1));
            
            % Stop playback
            PsychPortAudio('Stop', pahandle);
            
            % Close the audio device
            PsychPortAudio('Close', pahandle);
            
        end
        
        %% Registering responses %%

        if e_response == 1
            
            % Check the keyboard to see if a button has been pressed
            [keyIsDown,keyPress_t,keyCode] = KbCheck;
            
            % Log subject's reaction time (RT)
            RT = keyPress_t - tone_start_t;
            
            % Record up or down key based on participant input
            if keyCode(upKey)

                % Record "up" in results matrix
                choice = 1;
                e_response = 0;
                e_log = 1;
            elseif keyCode(downKey)

                % Record "down" in results matrix
                choice = -1;
                e_response = 0;
                e_log = 1;
            end
        end
        
        %% Log data and give warning if needed %%

        if e_log == 1

            %% Check for RT timeout %%

            if RT > rt_max
                rt_violation = 1;
                Screen('DrawText', window, 'PLEASE RESPOND FASTER!', xCenter,yCenter,[1 0 0]);
                Screen('Flip', window);
                WaitSecs(2); % punish with a delay
            elseif RT < rt_min
                rt_violation = 2;
                Screen('DrawText', window, 'PLEASE WAIT FOR THE TONE TO FINISH!', xCenter,yCenter,[1 0 0]);
                Screen('Flip', window); % flip
                WaitSecs(2); % punish with a delay
            else
                rt_violation = 0;
            end
            
            %% Save data

            % RT
            data.rt(trialnum) = RT;

            % Choice
            data.choice(trialnum) = choice;

            % RT violation
            data.rt_violation(trialnum) = rt_violation;

            % Trial start time
            data.trial_start_t(trialnum) = hold_start_t;

            % Stimulus presentation
            data.tone_start_t(trialnum) = tone_start_t;

            % Current experiment time
            data.experiment_time(trialnum) = GetSecs;

            % Variables for analysis
            data.corrparity(trialnum) = toneStim(trialnum).params.corrParity;
            data.displacement(trialnum) = toneStim(trialnum).params.displacement;
            data.coh(trialnum) = toneStim(trialnum).params.coh;
            
            %% Save everything %%

            save(['DataCoherence/',subInfo],'data'); % task data
            save(['DataCoherence/',subInfo,'_toneData'],'toneStim'); % stimulus data
            
            % Re-initialize loop
            trialnum = trialnum + 1; % iterate trial
            e_play_tone = 1; % time to start again
            hold_start_t = GetSecs();
            e_log = 0;
            
        end
    end
    
    % End message %
    Screen('DrawText', window, 'Great Job, Thanks!', xCenter-150,yCenter-150,[0 1 0]);
    Screen('Flip',window, [], 1);
    WaitSecs(3);
    clear Screen;
    ShowCursor;
    
    %% Zip all code and save for this subject

    fstruct = dir('*.m');
    for ii = 1:length(fstruct)
        code_filenames{ii}=fstruct(ii).name;
    end
    zip(['DataCoherence/' subInfo '_mfiles.zip'],code_filenames);
    
catch err % if 'try' hits a problem
    disp(err); % figure out what happened
    clear Screen % clear PTB screen if error
    Screen('CloseAll');
    ShowCursor;
    
    %% Zip all code and save for this subject
    fstruct = dir('*.m');
    for ii = 1:length(fstruct)
        code_filenames{ii}=fstruct(ii).name;
    end
    zip(['DataCoherence/' subInfo '_mfiles.zip'],code_filenames);
    
end






