% Simulation at 1 ms resolution of response of hypothetical units
% that detect individual high and low pairings in the sound

T = 1e7; % number time steps in ms

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


%% Figures shown in S2 and S3
%% qualitative results don't depend strongly on parameters, but the delay in 
%% the pattern must line up with imposed correlations to some degree
delay = 40;
pixelDur = 166;
[stimPos1,stimNeg1] = createStimFig1(T,pixelDur);
[stimPos2,stimNeg2] = createStimFig2(T,delay);
[stimPos3,stimNeg3] = createStimFig3(T,delay);

plotCountResponses(stimPos1,stimNeg1,delay,'Fig 1')
plotCountResponses(stimPos2,stimNeg2,delay,'Fig 2')
plotCountResponses(stimPos3,stimNeg3,delay,'Fig 3')
plotCountResponses(-stimPos3,-stimNeg3,delay,'Fig 3')



function plotCountResponses(stimPos,stimNeg,delay,str)

    T = size(stimPos,1);

    [stimPosCountHH(1),stimPosCountHH(2)] = runUnit(stimPos,delay,[1 1]);
    [stimPosCountHL(1),stimPosCountHL(2)] = runUnit(stimPos,delay,[1 -1]);
    [stimPosCountLL(1),stimPosCountLL(2)] = runUnit(stimPos,delay,[-1 -1]);
    [stimPosCountLH(1),stimPosCountLH(2)] = runUnit(stimPos,delay,[-1 1]);
    
    
    [stimNegCountHH(1),stimNegCountHH(2)] = runUnit(stimNeg,delay,[1 1]);
    [stimNegCountHL(1),stimNegCountHL(2)] = runUnit(stimNeg,delay,[1 -1]);
    [stimNegCountLL(1),stimNegCountLL(2)] = runUnit(stimNeg,delay,[-1 -1]);
    [stimNegCountLH(1),stimNegCountLH(2)] = runUnit(stimNeg,delay,[-1 1]);
    
    figure;
    subplot(2,2,1);
    bar([stimPosCountHH;stimPosCountLL;stimPosCountHL;stimPosCountLH]/T);
    title(['pos stimulus ' str])
    ylabel('occurence frequency (fraction of time)');
    xlabel('pattern counted')
    set(gca,'xtick',[1 2 3 4],'xticklabel',{'HH','LL','HL','LH'});
    legend('with','against')
    niceAxesLarge;
    
    subplot(2,2,2);
    bar([stimNegCountHH;stimNegCountLL;stimNegCountHL;stimNegCountLH]/T);
    title(['neg stimulus ' str])
    xlabel('pattern counted')
    set(gca,'xtick',[1 2 3 4],'xticklabel',{'HH','LL','HL','LH'});
    goodxlim = get(gca,'xlim');
    niceAxesLarge;

    subplot(2,2,3);
    bar([stimPosCountHH(1)-stimPosCountHH(2);stimPosCountLL(1)-stimPosCountLL(2);...
        stimPosCountHL(1)-stimPosCountHL(2);stimPosCountLH(1)-stimPosCountLH(2)]/T);    ylabel('occurence frequency (fraction of time)');
    xlabel('pattern counted')
    set(gca,'xtick',[1 2 3 4],'xticklabel',{'HH','LL','HL','LH'});
    set(gca,'xlim',goodxlim)
    niceAxesLarge;

    subplot(2,2,4);
    bar([stimNegCountHH(1)-stimNegCountHH(2);stimNegCountLL(1)-stimNegCountLL(2);...
        stimNegCountHL(1)-stimNegCountHL(2);stimNegCountLH(1)-stimNegCountLH(2)]/T);
    xlabel('pattern counted')
    set(gca,'xtick',[1 2 3 4],'xticklabel',{'HH','LL','HL','LH'});
    set(gca,'xlim',goodxlim)
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


function [count1,count2] = runUnit(stim,offset,pattern)
% counts how often particular patterns appear at specific offsets
    % counts instances in one direction and other direction
    count1 = sum( (stim(1:end-offset,1)==pattern(1)) .* (stim(1+offset:end,2)==pattern(2)) );
    count2 = sum( (stim(1:end-offset,2)==pattern(1)) .* (stim(1+offset:end,1)==pattern(2)) );

end

