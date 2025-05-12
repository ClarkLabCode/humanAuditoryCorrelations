%% Start of script to run auditory psychophysical experiments (binaural presentation)

clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%
%% BEGIN PREP CODE %%
%%%%%%%%%%%%%%%%%%%%%

rng('shuffle'); % set seed for rand

%% First, let's get the subject's ID and store it

% Addpath data
subInfo = input('Enter the subject ID #:', 's');
subInfo = strcat(subInfo); % concatenate string entered
while exist(['DataBinaural/',subInfo,'.mat'],'file') % make sure no overwriting is happening
    subInfo = input('Subject # already exists! Enter new #:', 's'); % re-do if already used
    subInfo = strcat(subInfo);
end
mkdir(['stimuli/',subInfo]);

%% Second, we create the stimlus sequence

main_trial_types = 1:4; % how many different "types" of trials of interest (not including random probes)

mainN = 20; % number of trials for each "main" trial type in experiment

randomN = 20; % number number of trials of random "probe" stimuli

tmpSeq = []; % initialize empty vector

% Build sequence vector
for j = main_trial_types
    tmpSeq = [ones(mainN,1)*main_trial_types(j);tmpSeq]; % vector of N trials of trial type j, appended onto itself each loop
end

% Add random trials on to the end
tmpSeq = [ones(randomN,1)*main_trial_types(end)+1;tmpSeq];

% Shuffle the vector
trialSeq = Shuffle(tmpSeq); % randomization function ("shuffle")

nTrials = length(trialSeq); % total task N

%% Now that we have our shuffled sequence, lets pre-save the actual stimuli as WAVs

% Set the FIXED stimulus parameters first (these are key psychophysical variables)
params.stimDur = 2; % duration of stimulus (in seconds)
params.nTonesPerOctave = 15; % number notes in each octave
params.noteDur = 1/6; % note duration (in seconds)
params.deltaNote = 1; % delta note
params.deltaT = 1; % delta time
Fs = 2e4; % freq

%% Now we use a for-loop to generate the tones using our stimulus sequence

for n = 1:nTrials
    
    if trialSeq(n) == 1 % TRIAL TYPE 1
        params.corrType = 'binauralScint';
        params.corrParity = 1; % positive correlation
        params.displacement = 1; % upward direction
        
    elseif trialSeq(n) == 2 % TRIAL TYPE 2
        params.corrType = 'binauralScint';
        params.corrParity = 1; % positive correlation
        params.displacement = -1; % downward direction
        
    elseif trialSeq(n) == 3 % TRIAL TYPE 3
        params.corrType = 'binauralScint';
        params.corrParity = -1; % negative correlation
        params.displacement = 1; % upward direction
        
    elseif trialSeq(n) == 4 % TRIAL TYPE 4
        params.corrType = 'binauralScint';
        params.corrParity = -1; % negative correlation
        params.displacement = -1; % downward direction
        
    elseif trialSeq(n) == 5 % TRIAL TYPE 5
        params.corrType = 'binauralScintRand'; % random
        params.corrParity = 1;
        params.displacement = 1;
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Call tone generating function %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Generate trial n's stimulus
    toneStim(n) = binauralFunction(params);
    
    % Set wav filename
    soundFilename{n} = ['stimuli/',subInfo,'/',subInfo,'_tone',num2str(n),'.wav'];

    % Rescale to prevent distorted noises
    Y1 = toneStim(n).waveFormLeft;
    Y1 = rescale(Y1,-1,1);

    Y2 = toneStim(n).waveFormRight;
    Y2 = rescale(Y2,-1,1);

    % Write out WAV
    audiowrite(soundFilename{n},[Y1(:),Y2(:)],Fs);
    
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

    % Get the center coordinate of the window in pixels
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
        'corrparity',nan(nTrials,1),'displacement',nan(nTrials,1));
    
    % Define avaliable keys to press
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
        
        %% trial has started
        
        if e_play_tone == 1 % enter the "play_tone" phase

            % show a fixation cross
            Screen('DrawLines', window, allCoords,...
                lineWidthPix, white, [xCenter yCenter], 2);
            Screen('Flip',window);
            
            %% Load sound for this trial %%

            % Specify filename of sound to load
            wavfilename = soundFilename{trialnum}; % load sound associated with this trial, from subject's directory in "stimuli"

            % Read WAV file from filesystem
            [y, freq] = psychwavread(wavfilename);
            wavedata = y';
            nrchannels = size(wavedata,1); % number of rows == number of channels
            
            % Make sure we have always 2 channels stereo output (PPA functioning issue)
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
                Screen('Flip', window);
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
            data.displacement(trialnum) = toneStim(trialnum).params.displacement;
            data.corrparity(trialnum) = toneStim(trialnum).params.corrParity;
            
            %% Save everything %%

            save(['DataBinaural/',subInfo],'data'); % task data
            save(['DataBinaural/',subInfo,'_toneData'],'toneStim'); % stimulus data
            
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
    zip(['DataBinaural/' subInfo '_mfiles.zip'],code_filenames);
    
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
    zip(['DataBinaural/' subInfo '_mfiles.zip'],code_filenames);
    
end






