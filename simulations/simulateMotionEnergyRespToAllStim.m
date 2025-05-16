% Simulation at 1 ms resolution of response of hypothetical motion energy
% units
% ME model as in Figure 4c.
% Stimuli as in Fig 1, 2, and 3

T = 1e6; % number time steps = 1000 s


%% Show sample stimuli

sampleDelays = [20 40];
figure;
for ii=1:length(sampleDelays)
    [sp,sn]=createStimFig3(2000,sampleDelays(ii));
    subplot(length(sampleDelays),2,2*(ii-1)+1);
    imagesc(sp'); ylabel('pos'); title(['F3 delay = ' num2str(sampleDelays(ii))]);
    colormap('gray');
    set(gca,'ydir','normal');
    subplot(length(sampleDelays),2,2*(ii-1)+2);
    imagesc(sn'); ylabel('neg');
    colormap('gray')
    set(gca,'ydir','normal');
end

%%

figure;
[sp,sn]=createStimFig1(2000,20);
subplot(2,1,1);
imagesc(sp'); ylabel('pos'); title(['fig 1 sample']);
colormap('gray')
set(gca,'ydir','normal');
subplot(2,1,2);
imagesc(sn'); ylabel('neg');
colormap('gray')
set(gca,'ydir','normal');

%%

figure;
[sp,sn]=createStimFig2(2000,20);
subplot(2,1,1);
imagesc(sp'); ylabel('pos'); title(['fig 2 sample']);
colormap('gray')
set(gca,'ydir','normal');
subplot(2,1,2);
imagesc(sn'); ylabel('neg');
colormap('gray')
set(gca,'ydir','normal');


%%

delay = 40;
pixelDur = round(1000/6);
[stimPos1,stimNeg1] = createStimFig1(T,pixelDur);
[stimPos2,stimNeg2] = createStimFig2(T,delay);
[stimPos3,stimNeg3] = createStimFig3(T,delay);

% These simulated responses are plotted in figure S4
% qualitatively, will get this sort of pattern no matter the details of the
% linear filter
offset = 1;
plotMEResponses(offset+stimPos1,offset+stimNeg1,'Fig 1')
plotMEResponses(offset+stimPos2,offset+stimNeg2,'Fig 2')
plotMEResponses(offset+stimPos3,offset+stimNeg3,'Fig 3 ++ and +-')
plotMEResponses(offset-stimPos3,offset-stimNeg3,'Fig 3 -- and -+')



function plotMEResponses(stimPos,stimNeg,str)

    % set up for filtering
    T = size(stimPos,1);

    Tfilt = 500; % duration of filter in ms
    tau1 = 40; % timescale in ms
    delayFilt = 40; % delay between filters in ms
    
    t = [0:Tfilt-1];
    fshape = t.^2.*exp(-t/tau1*2); % get peak at tau1
    fshape = fshape/max(fshape);
    shapeChoose = 'monolobed'; % 'bilobed','monolobed','derivative'
    switch shapeChoose
        case 'bilobed'
            f1 = [0,diff(fshape)]; %
            f1 = f1/max(f1);
            f2 = [zeros(1,delayFilt),f1(1:end-delayFilt)];
        case 'monolobed'
            f1 = fshape;
            f2 = [zeros(1,delayFilt),f1(1:end-delayFilt)];
        case 'derivative'
            f1 = [0,diff(fshape)];
            f1 = f1/max(f1);
            f2 = fshape;
    end



    meanMEpos(1) = runMEUnit(f1,f2,'squared',stimPos);
    meanMEpos(2) = runMEUnit(f2,f1,'squared',stimPos);

    meanMEneg(1) = runMEUnit(f1,f2,'squared',stimNeg);
    meanMEneg(2) = runMEUnit(f2,f1,'squared',stimNeg);
    
    
    figure;
    subplot(1,2,1);
    bar([meanMEpos(2)-meanMEpos(1),meanMEneg(2)-meanMEneg(1)]);
    title(['stimulus ' str])
    ylabel('ME output (au)');
    xlabel('stim type')
    set(gca,'xtick',[1 2],'xticklabel',{'positive','negative'});
    niceAxesLarge;
    
end


function [stimPos,stimNeg] = createStimFig1(T, duration);
    % create ternary stimulus a la figure 1
    N = ceil(T/duration); % must be whole number please

    b12 = (rand(N,1)<0.5)*2-1; % one binary to share
    b2 = (rand(N,1)<0.5)*2-1;
    b1 = (rand(N,1)<0.5)*2-1;

    stimPos(:,1) = (b12+b1)/2;
    stimPos(:,2) = (circshift(b12,[1 0])+b2)/2;

    stimNeg(:,1) = (b12+b1)/2;
    stimNeg(:,2) = (b2-circshift(b12,[1 0]))/2;

    stimNeg = upsample(stimNeg,duration);
    stimPos = upsample(stimPos,duration);

    % make the pips have duration dur
    stimPos = filter(ones(1,duration),1,stimPos);
    stimNeg = filter(ones(1,duration),1,stimNeg);


end

function [stimPos,stimNeg] = createStimFig3(T, delay)
% creates a set of stimuli with different delays and parities

    duration = 20; % ms, pip duration
    pipRate = 4; % Hz, correlated pip-pairs per second and 1x uncorrelated pips
    
    stim = zeros(T,2);
    stimPos = stim;
    stimNeg = stim;
    
    initPipInds = find(rand(size(stim(:,1)))<pipRate/1000);

    randPipInds12 = find(rand(size(stim))<pipRate*1/1000); % matches way stimulus was set up for these ones -- we halved the uncorrelated pips compared to correlated ones, likely unintentionally
    randPipC12 = double(rand(size(randPipInds12))>0.5)*2 - 1;

    randPipInds1 = find(rand(size(stim(:,1)))<pipRate*1/1000); % these are pips showing up as the second in the correlated pair coming from row 0
    randPipC1 = double(rand(size(randPipInds1))>0.5)*2 - 1;
    
    randPipInds2 = find(rand(size(stim(:,2)))<pipRate*1/1000); % these are pips showing up as the first in the correlated pair going to row 3
    randPipC2 = double(rand(size(randPipInds2))>0.5)*2 - 1;
    
    
    stimPos(initPipInds,1) = 1; % set the correlated pips
    stimNeg(initPipInds,1) = 1; % these will be negative correlated
    
    D = delay;
    stimPos(:,2) = circshift(stimPos(:,1),[D,0]);
    stimNeg(:,2) = -circshift(stimNeg(:,1),[D,0]);
    
    % set the uncorrelated pips
    stimPos(randPipInds12) = -1; % set all random pips to low
    stimNeg(randPipInds12) = randPipC12; % keep these random high-low

    % these are pips corresponding to the correlated ones coming into row 1
    % from row 0, so they are +1 for HH and -1 for HL
    stimPos(randPipInds1,1) = 1;
    stimNeg(randPipInds1,1) = -1;

    % these are pip corresponding to initial corrlated pip pairs initiating
    % in row 2, going to row 3
    stimPos(randPipInds2,2) = 1;
    stimNeg(randPipInds2,2) = 1;
    

    % make the pips have duration dur
    stimPos = filter(ones(1,duration),1,stimPos);
    stimNeg = filter(ones(1,duration),1,stimNeg);

    % get rid of stuff above 1 or below -1
    stimPos(stimPos>1) = 1;
    stimNeg(stimNeg>1) = 1;
    stimPos(stimPos<-1) = -1;
    stimNeg(stimNeg<-1) = -1;

end

function [stimPos,stimNeg] = createStimFig2(T, delay)
% creates a set of stimuli with different delays and parities

    duration = 20; % ms, pip duration
    pipRate = 4; % Hz, correlated pip-pairs per second and 1x uncorrelated pips
    
    stim = zeros(T,2);
    stimPos = stim;
    stimNeg = stim;
    
    initPipInds = find(rand(size(stim(:,1)))<pipRate/1000);
    initPipC = double(rand(size(initPipInds))>0.5)*2 - 1;

    randPipInds1 = find(rand(size(stim(:,1)))<pipRate*1/1000); % these are pips showing up as the second in the correlated pair coming from row 0
    randPipC1 = double(rand(size(randPipInds1))>0.5)*2 - 1;
    
    randPipInds2 = find(rand(size(stim(:,2)))<pipRate*1/1000); % these are pips showing up as the first in the correlated pair going to row 3
    randPipC2 = double(rand(size(randPipInds2))>0.5)*2 - 1;
    
    
    stimPos(initPipInds,1) = initPipC; % set the correlated pips
    stimNeg(initPipInds,1) = initPipC; % these will be negative correlated
    
    D = delay;
    stimPos(:,2) = circshift(stimPos(:,1),[D,0]);
    stimNeg(:,2) = -circshift(stimNeg(:,1),[D,0]);
    
    % these are pips corresponding to the correlated ones coming into row 1
    % from row 0, so they are +1 for HH and -1 for HL
    stimPos(randPipInds1,1) = randPipC1;
    stimNeg(randPipInds1,1) = randPipC1;

    % these are pip corresponding to initial corrlated pip pairs initiating
    % in row 2, going to row 3
    stimPos(randPipInds2,2) = randPipC2;
    stimNeg(randPipInds2,2) = randPipC2;
    

    % make the pips have duration dur
    stimPos = filter(ones(1,duration),1,stimPos);
    stimNeg = filter(ones(1,duration),1,stimNeg);

    % get rid of stuff above 1 or below -1
    stimPos(stimPos>1) = 1;
    stimNeg(stimNeg>1) = 1;
    stimPos(stimPos<-1) = -1;
    stimNeg(stimNeg<-1) = -1;

end

function meanResp = runMEUnit(f1,f2,nonlin,stim)

    % computed mean response of a unit based on parameters and stim

    rlin = filter(f1,1,stim(:,1)) + filter(f2,1,stim(:,2));
    switch nonlin
        case 'squared'
            resp = rlin.^2;
        case 'relu'
            resp = rlin;
            resp(resp<0) = 0;
        case 'relu2'
            resp = rlin;
            resp(resp<0) = 0;
            resp = resp.^2;
        otherwise
            disp('bad nonlinearity supplied');
            return;
    end
    meanResp = mean(resp(1000:end)); % ignore first second, generically


end



