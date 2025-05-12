% Analyze supplementary pip delta time data (20 ms pips)

%% Loop through subjects and pull out data 

numsubs = 8;
counter=1;
for subnum = 1:numsubs
    load(['s00', num2str(subnum),'.mat']);
    uT = unique(data.deltaT);
    disp(['unique times = ' num2str(uT')]);
    counter2 = 1;
    for C = [1,-1]
        for D = [1 -1]
            for ii=1:length(uT)
                ff = find( (data.deltaT == uT(ii)) .* (data.corrparity == C) .* (data.displacement == D));
                meanChoice(ii,counter,((C+1)+(D+1)/2)+1) = mean(data.choice(ff));
                semChoice(ii,counter,((C+1)+(D+1)/2)+1) = std(data.choice(ff))/sqrt(length(ff));
            end
            counter2 = counter2 + 1; 
        end
    end
       counter=counter+1;
end

%% Plot positive correlation results

% Set x-axis values
t=[0 10 20 40 80 160 320];

% Turn anything that has negative direction into something that has
% negative time
[sortedPlusMinusT, sortedPlusMinusTIdx] = sort([-t, t]);

% Set plotting for positive or negative correlations here: (1,2) = -C
% (negative correlation), (3,4) = +C (positive correlation)
meanChoicePlusMinus = [(meanChoice(:,:,3)'+1)/2,(meanChoice(:,:,4)'+1)/2];
sortedmeanChoicePlusMinus = meanChoicePlusMinus(:,sortedPlusMinusTIdx);

% Take some steps to combine positive and negative time symmetrically
indZero = find(sortedPlusMinusT == 0); % give us the indices where == 0
VecZero = (sortedmeanChoicePlusMinus(:,indZero(1)) + sortedmeanChoicePlusMinus(:,indZero(2)))/2;
sortedmeanChoicePlusMinus(:,indZero(2)) = [];
sortedmeanChoicePlusMinus(:,indZero(1)) = VecZero;
sortedPlusMinusT(:,indZero(1)) = [];

% Calculate SEM
mC2=sortedmeanChoicePlusMinus;
sem=std(mC2,[],1)./sqrt(size(sortedmeanChoicePlusMinus,1));

% Define x-values
x_values = [-6:6];

% Plot results
figure;

% Plot individual data in faint, gray lines
plot(x_values,sortedmeanChoicePlusMinus', 'color',[0.8 0.8 0.8]);
hold on;

% Create y-values that cover standard error above and below mean
y1 = mean(sortedmeanChoicePlusMinus,1) + sem; % above mean
y2 = mean(sortedmeanChoicePlusMinus,1) - sem; % below mean

% plot mean + SEM line
plot(x_values(:),y1, 'w')
hold on;

% plot mean - SEM line
plot(x_values(:),y2, 'w')

% plot mean line
plot(x_values,(mean(sortedmeanChoicePlusMinus,1)), 'o-', 'LineWidth', 2, 'color', [0.83 0.44 0.5]);

% Shade in error region
patch([x_values(:); flipud(x_values(:))],[y1(:); flipud(y2(:))], 'r', 'FaceAlpha',0.1, 'EdgeColor','none');
hold off;

% Set titles and labels
xlabel('time interval (ms)')
ylabel('probability perceived up')
hold on
yline(0.5,'k','linewidth',2)
set(gca,'linewidth',3,'fontSize',16)
xticks(x_values);
xlim([-6 6]);
ylim([0 1]);
xticklabels({'-320','-160','-80','-40','-20','-10' '0','10','20','40','80','160','320'});
title('positive correlations');


%% Plot negative correlation results

% Set x-axis values
t=[0 10 20 40 80 160 320];

% Turn anything that has negative direction into something that has
% negative time
[sortedPlusMinusT, sortedPlusMinusTIdx] = sort([-t, t]);

% Set plotting for positive or negative correlations here: (1,2) = -C
% (negative correlation), (3,4) = +C (positive correlation)
meanChoicePlusMinus = [(meanChoice(:,:,1)'+1)/2,(meanChoice(:,:,2)'+1)/2];
sortedmeanChoicePlusMinus = meanChoicePlusMinus(:,sortedPlusMinusTIdx);

% Take some steps to combine positive and negative time symmetrically
indZero = find(sortedPlusMinusT == 0); % give us the indices where == 0
VecZero = (sortedmeanChoicePlusMinus(:,indZero(1)) + sortedmeanChoicePlusMinus(:,indZero(2)))/2;
sortedmeanChoicePlusMinus(:,indZero(2)) = [];
sortedmeanChoicePlusMinus(:,indZero(1)) = VecZero;
sortedPlusMinusT(:,indZero(1)) = [];

% Calculate SEM
mC2=sortedmeanChoicePlusMinus;
sem=std(mC2,[],1)./sqrt(size(sortedmeanChoicePlusMinus,1));

% Define x-values
x_values = [-6:6];

% Plot results
figure;

% Plot individual data in faint, gray lines
plot(x_values,sortedmeanChoicePlusMinus', 'color',[0.8 0.8 0.8]);
hold on;

% Create y-values that cover standard error above and below mean
y1 = mean(sortedmeanChoicePlusMinus,1) + sem; % above mean
y2 = mean(sortedmeanChoicePlusMinus,1) - sem; % below mean

% plot mean + SEM line
plot(x_values(:),y1, 'w')
hold on;

% plot mean - SEM line
plot(x_values(:),y2, 'w')

% plot mean line
plot(x_values,(mean(sortedmeanChoicePlusMinus,1)), 'o-', 'LineWidth', 2, 'color', [0.83 0.44 0.5]);

% Shade in error region
patch([x_values(:); flipud(x_values(:))],[y1(:); flipud(y2(:))], 'r', 'FaceAlpha',0.1, 'EdgeColor','none');
hold off;

% Set titles and labels
xlabel('time interval (ms)')
ylabel('probability perceived up')
hold on
yline(0.5,'k','linewidth',2)
set(gca,'linewidth',3,'fontSize',16)
xticks(x_values);
xlim([-6 6]);
ylim([0 1]);
xticklabels({'-320','-160','-80','-40','-20','-10' '0','10','20','40','80','160','320'});
title('negative correlations');

%% One-way, repeated measures ANOVA for positive correlation

% Load responses
load('positive_correlation.mat') % called "positive_correlation" in the workspace

% Rename things
pos_responses = positive_correlation;
data = pos_responses;

% Make a table
table = array2table(data);
table.Properties.VariableNames = {'pos1', 'pos2', 'pos3', 'pos4', 'pos5', 'pos6', 'pos7', 'pos8', 'pos9', 'pos10', 'pos11', 'pos12', 'pos13'};

% Define delay times
ts = [1 2 3 4 5 6 7 8 9 10 11 12 13]';

other_structure = ts;

withinDesign = array2table(other_structure);

withinDesign.Properties.VariableNames = {'DeltaTs'};
withinDesign.DeltaTs = categorical(withinDesign.DeltaTs);

% Fit model
rm = fitrm(table, 'pos1-pos13 ~ 1', 'WithinDesign', withinDesign);

% Save results
result = ranova(rm);

%% One-way, repeated measures ANOVA for negative correlation

% Load responses
load('negative_correlation.mat') % called "negative_correlation" in the workspace

% Rename things
neg_responses = negative_correlation;
data = neg_responses;

% Make a table
table = array2table(data);
table.Properties.VariableNames = {'neg1', 'neg2', 'neg3', 'neg4', 'neg5', 'neg6', 'neg7', 'neg8', 'neg9', 'neg10', 'neg11', 'neg12', 'neg13'};

% Define delay times
ts = [1 2 3 4 5 6 7 8 9 10 11 12 13]';

other_structure = ts;

withinDesign = array2table(other_structure);

withinDesign.Properties.VariableNames = {'DeltaTs'};
withinDesign.DeltaTs = categorical(withinDesign.DeltaTs);

% Fit model
rm = fitrm(table, 'neg1-neg13 ~ 1', 'WithinDesign', withinDesign);

% Save results
result = ranova(rm);