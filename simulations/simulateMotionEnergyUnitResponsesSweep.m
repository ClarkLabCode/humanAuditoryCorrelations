% Simulation at 1 ms resolution of mean response of hypothetical unit with
% ME-like properties

T = 1e6; % number time steps = 100 s

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

figure;
subplot(2,1,1);
plot(t,f1,'r',t,f2,'k');
xlabel('time (ms)');
ylabel('filter amplitude');
legend('f1','f2');
set(gca,'xlim',[0 300])

subplot(2,1,2);
imagesc([f1;f2]);
xlabel('time (ms)');
ylabel('frequency');
colormap(makeRBcolormap);
set(gca,'clim',[-1 1]*max(abs(f1)));
set(gca,'xlim',[0 300],'ydir','normal')

figure;
imagesc([0:128]);
colormap(makeRBcolormap);
set(gca,'clim',[-128 128]);

figure; hold on;
x = [-1:.01:1];
plot(x,x.^2,'b');
plot([-1 1],[0 0],'k-');
plot([0 0],[-1 1],'k-');
set(gca,'xlim',[-1 1],'ylim',[-0.5 1.1]);


%% Show sample stimuli

sampleDelays = [20 40];
figure;
for ii=1:length(sampleDelays)
    [sp,sn]=createStim(2000,sampleDelays(ii));
    subplot(length(sampleDelays),2,2*(ii-1)+1);
    imagesc(sp{1}'); ylabel('pos'); title(['delay = ' num2str(sampleDelays(ii))]);
    colormap('gray')
    set(gca,'ydir','normal');
    subplot(length(sampleDelays),2,2*(ii-1)+2);
    imagesc(sn{1}'); ylabel('neg');
    colormap('gray')
    set(gca,'ydir','normal');
end


%% Make stimuli

delayVals = [-150:5:150];
[stimPos,stimNeg] = createStim(T,delayVals);


for ii=1:length(stimPos)
    respMeanPos(ii) = runUnit(f2,f1,'squared',stimPos{ii});
    respMeanNeg(ii) = runUnit(f2,f1,'squared',stimNeg{ii});
end

figure; 
subplot(2,1,1);
plot(delayVals,respMeanPos,'g.-','markersize',20); hold on;
plot(delayVals,respMeanNeg,'m.-','markersize',20);
legend('pos','neg');
xlabel('delay (ms)');
ylabel('response amplitude');
set(gca,'ylim',[0 max(get(gca,'ylim'))]);

subplot(2,1,2);
plot(delayVals,respMeanPos-respMeanPos(end:-1:1),'g.-','markersize',20); hold on;
plot(delayVals,respMeanNeg-respMeanNeg(end:-1:1),'m.-','markersize',20);
legend('pos','neg');
xlabel('delay (ms)');
ylabel('opponent response amplitude');
set(gca,'ylim',1.1*get(gca,'ylim'));


function [stimPos,stimNeg] = createStim(T, delays)
% creates a set of stimuli with different delays and parities

    dur = 20; % 10 ms pips
    pipRate = 4; %  correlated pip-pairs per second and 2x uncorrelated pips
    
    stim = zeros(T,2);
    
    initPipInds = find(rand(size(stim(:,1)))<pipRate/1000); 
    initPipC = double(rand(size(initPipInds))>0.5)*2 - 1; % contrasts between -1 and 1
    randPipInds = find(rand(size(stim))<pipRate*2/1000);
    randPipC = double(rand(size(randPipInds))>0.5)*2 - 1;
    
    stim(initPipInds,1) = initPipC; % set the correlated pips
    
    for ii=1:length(delays)
        D = delays(ii);
        stimPos{ii} = stim;
        stimNeg{ii} = stim;
        stimPos{ii}(:,2) = circshift(stimPos{ii}(:,1),[D,0]);
        stimNeg{ii}(:,2) = -circshift(stimNeg{ii}(:,1),[D,0]);
    
        % set the uncorrelated pips
        stimPos{ii}(randPipInds) = randPipC;
        stimNeg{ii}(randPipInds) = randPipC;
    
        % make the pips have duration dur
        stimPos{ii} = filter(ones(1,dur),1,stimPos{ii});
        stimNeg{ii} = filter(ones(1,dur),1,stimNeg{ii});
    
        % get rid of stuff above 1 or below -1
        stimPos{ii}(stimPos{ii}>1) = 1;
        stimNeg{ii}(stimNeg{ii}>1) = 1;
        stimPos{ii}(stimPos{ii}<-1) = -1;
        stimNeg{ii}(stimNeg{ii}<-1) = -1;
    end

end

function resp = runUnit(f1,f2,nonlin,stim)
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
    resp = mean(resp(1000:end)); % ignore first second, generically

end

