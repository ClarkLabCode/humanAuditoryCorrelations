% Analyze correlated pip data

%% Loop through subjects and pull out data 

numsubs = 10;
counter=1;
for subnum = 1:numsubs
    load(['subj', num2str(subnum),'.mat']);
    uP = unique(data.corrPipTypeNumber);
    disp(['unique notes = ' num2str(uP')]);
counter2 = 1;
for D = [1,-1]
    for ii=1:length(uP)
        ff = find((data.corrPipTypeNumber == uP(ii)) .* (data.direction == D));
        meanChoice(ii,counter,((D+3)/2)) = mean(data.choice(ff));
        semChoice(ii,counter,((D+3)/2)) = std(data.choice(ff))/sqrt(length(ff));
    end
    counter2 = counter2 + 1;  
end
counter = counter + 1;
end

%% Plot upward-directed results

% Set direction (+1 for upward-directed, -1 for downward-directed)
D=1;

% Turn -1 into 1 and 1 into 2
dir=(D+3)/2;

% Calculate SEMs
mC2=(meanChoice(:,:,dir)+1)/2;
sem=std(mC2,[],2)./sqrt(numsubs);

% Define x-axis values
conditions = [1:4];

% Consolidate things
data = ((meanChoice(:,:,dir)+1)/2)';
sems = sem;

% Caulcate means
group_means = mean(data);

% Define number of conditions and subjects
num_conditions = size(data, 2);
num_subjects = size(data, 1);

% Plot results
figure;
bar_handle = bar(conditions, group_means, 'basevalue',0.5);
hold on;
errorbar(conditions, group_means, sems, 'LineStyle', 'none', 'color', [0 0 0])
hold on;

% Plot individual subjects with dots and connecting lines
for i = 1:num_subjects
    x = 1:num_conditions;
    y = data(i, :);
    scatter(x, y, 'k', 'filled')
    plot(x, y, '-o', 'color', [0 0 0])
end
hold off;

yline(0.5,'k','linewidth',2)
set(gca,'linewidth',3,'fontSize',16)
xlabel('stimulus type')
ylabel('probability perceived rising')
title('upward-directed')
xticklabels({'(+,+)','(+,-)', '(-,+)', '(-,-)'})

%% Plot downward-directed results

% Set direction (+1 for upward-directed, -1 for downward-directed)
D=-1;

% Turn -1 into 1 and 1 into 2
dir=(D+3)/2;

% Calculate SEMs
mC2=(meanChoice(:,:,dir)+1)/2;
sem=std(mC2,[],2)./sqrt(numsubs);

% Define x-axis values
conditions = [1:4];

% Consolidate things
data = ((meanChoice(:,:,dir)+1)/2)';
sems = sem;

% Caulcate means
group_means = mean(data);

% Define number of conditions and subjects
num_conditions = size(data, 2);
num_subjects = size(data, 1);

% Plot results
figure;
bar_handle = bar(conditions, group_means, 'basevalue',0.5);
hold on;
errorbar(conditions, group_means, sems, 'LineStyle', 'none', 'color', [0 0 0])
hold on;

% Plot individual subjects with dots and connecting lines
for i = 1:num_subjects
    x = 1:num_conditions;
    y = data(i, :);
    scatter(x, y, 'k', 'filled')
    plot(x, y, '-o', 'color', [0 0 0])
end
hold off;

yline(0.5,'k','linewidth',2)
set(gca,'linewidth',3,'fontSize',16)
xlabel('stimulus type')
ylabel('probability perceived rising')
title('downward-directed')
xticklabels({'(+,+)','(+,-)', '(-,+)', '(-,-)'})

%% Perform stats: two sample t-tests for upward- vs. downward-directed stimuli

% Load and organize data
load('up_data.mat') % called "up_data" in the workspace
load('down_data.mat') % called "down_data" in the workspace

% Rename things
up = up_data;
down = down_data;

up_pos_pos = up(:,1);
up_pos_neg = up(:,2);
up_neg_pos = up(:,3);
up_neg_neg = up(:,4);

down_pos_pos = down(:,1);
down_pos_neg = down(:,2);
down_neg_pos = down(:,3);
down_neg_neg = down(:,4);

% Perform the two-sample t-tests
[h1_twosample,p1_twosample,ci1_twosample,stats1_twosample] = ttest2(up_pos_pos, down_pos_pos);
[h2_twosample,p2_twosample,ci2_twosample,stats2_twosample] = ttest2(up_pos_neg, down_pos_neg);
[h3_twosample,p3_twosample,ci3_twosample,stats3_twosample] = ttest2(up_neg_pos, down_neg_pos);
[h4_twosample,p4_twosample,ci4_twosample,stats4_twosample] = ttest2(up_neg_neg, down_neg_neg);

% Save results
results_twosample = [p1_twosample, p2_twosample, p3_twosample, p4_twosample];

%% Analyze net results (P(rising | upward directed) – P(rising | downward directed))

% Set up direction
D=1;

% Turn -1 into 1 and 1 into 2
dir=(D+3)/2;

% Calculate SEMs
mC2=(meanChoice(:,:,dir)+1)/2;
sem_rising=std(mC2,[],2)./sqrt(10);

% Define rising data and compute mean
data_rising = ((meanChoice(:,:,dir)+1)/2)';
group_means_rising = mean(data_rising);

% Set down direction
D=-1;

% Turn -1 into 1 and 1 into 2
dir=(D+3)/2;

% Calculate SEMs
mC2=(meanChoice(:,:,dir)+1)/2;
sem_falling=std(mC2,[],2)./sqrt(10);

% Define falling data and compute mean
data_falling = ((meanChoice(:,:,dir)+1)/2)';
group_means_falling = mean(data_falling);

% Compute net result, mean, and SEM
data_net = data_rising - data_falling;
group_mean_net = mean(data_net);
sem_net = std(data_net) ./ sqrt(numsubs);

% Define number of conditions
num_conds = size(data_net, 2);

% Plot results
figure;
bar_handle = bar(conditions, group_mean_net, 'basevalue',0);
hold on;
errorbar(conditions, group_mean_net, sem_net, 'LineStyle', 'none', 'color', [0 0 0])
hold on;

% Plot individual subjects with dots and connected lines
for i = 1:numsubs
    x = 1:num_conds;
    y = data_net(i, :);
    scatter(x, y, 'k', 'filled')
    plot(x, y, '-o', 'color', [0 0 0])
end
hold off;

yline(0,'k','linewidth',2)
set(gca,'linewidth',3,'fontSize',16)
xlabel('net condition')
ylabel('probability perceived rising')
title('net')
xticklabels({'(+,+)','(+,-)', '(-,+)', '(-,-)'})