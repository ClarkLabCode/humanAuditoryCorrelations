% Analyze ternary stimuli experiment

%% Organize data by each trial type

% Define sample size
subs = [1:10];
nsubs = length(subs);

% Loop through subjects
for j = 1:nsubs
    load(['s00',num2str(j),'_toneData.mat']); % load 'rich' file
    load(['s00',num2str(j),'.mat']); % load 'diet' file

    DF = data; % rename data frames

    % Pull out "random" probe trials
    for k = 1:320
        randTrial_boolean(k) = toneStim(k).params.corrType=="ternRandom";
    end

   
    DF.choice(DF.choice==-1) = 0; % Convert "falling" pitch into zero (instead of -1)

    POSavgs(j,:) = [mean(DF.choice(DF.corrparity==1 & DF.displacement==1 & randTrial_boolean'==0)) ... % +1 corr, +1 disp

        mean(DF.choice(DF.corrparity==1 & DF.displacement==-1 & randTrial_boolean'==0))]; % +1 corr, -1 disp

    NEGavgs(j,:) = [mean(DF.choice(DF.corrparity==-1 & DF.displacement==1 & randTrial_boolean'==0)) ... % -1 corr, +1 disp

        mean(DF.choice(DF.corrparity==-1 & DF.displacement==-1 & randTrial_boolean'==0))]; % -1 corr, -1 disp

    RANDavg(j) = mean(DF.choice(randTrial_boolean'==1)); % random trials!

end

%% Calculate means, standard deviations, and SEMs

% Compute means
PosGroupMean = mean(POSavgs);
NegGroupMean = mean(NEGavgs);
RandGroupMean = mean(RANDavg);

% Compute SEMs
sem1 = std(POSavgs(:,1))./sqrt(nsubs);
sem2 = std(POSavgs(:,2))./sqrt(nsubs);
sem3 = std(NEGavgs(:,1))./sqrt(nsubs);
sem4 = std(NEGavgs(:,2))./sqrt(nsubs);
sem5 = std(RANDavg(1,:))./sqrt(nsubs);

%% Plot group averages and individual subjects with dots/lines

% Define x-ticks
conditions = [1:5];

% Group data
data = horzcat(POSavgs, NEGavgs, RANDavg');
group_means = mean(data);
sems = [sem1, sem2, sem3, sem4, sem5];

% Define number of trial types and subjects
num_conditions = size(data, 2);
num_subjects = size(data, 1);

% Plot bar graph
figure;
bar_handle = bar(conditions, group_means, 'basevalue',0.5);
xticklabels({'up pos', 'down pos', 'up neg', 'down neg', 'rand'});
hold on;

% Errorbars
errorbar(conditions, group_means, sems, 'LineStyle', 'none', 'color', [0 0 0])
hold on;

% Individual subject dots and connecting lines
for i = 1:num_subjects
    x = 1:num_conditions;
    y = data(i, :);
    scatter(x, y, 'k', 'filled')
    plot(x, y, '-o', 'color', [0 0 0])
end
hold off;

yline(0.5,'k','linewidth',2) % line at random chance probability
set(gca,'linewidth',3,'fontSize',16)
xlabel('stimulus type')
ylabel('probability perceived rising')

%% Perform stats: one sample t-test see if each result is different from 0.5

% Re-organize data
up_pos = data(:,1);
down_pos = data(:,2);
up_neg = data(:,3);
down_neg = data(:,4);
random = data(:,5);

% Mean to test against
mu = 0.5;

% Perform the one-sample t-tests
[h1, p1, ci1, stats1] = ttest(up_pos, mu);
[h2, p2, ci2, stats2] = ttest(down_pos, mu);
[h3, p3, ci3, stats3] = ttest(up_neg, mu);
[h4, p4, ci4, stats4] = ttest(down_neg, mu);
[h5, p5, ci5, stats5] = ttest(random, mu);

% Consolidate results
results = [p1, p2, p3, p4, p5];

%% Make a scatterplot to see how well the probability of percept is in the direction of displacement (Figure S1)

% Define positive correlations in data
up_pos = POSavgs(:,1);
down_pos = POSavgs(:,2);

% Initialize empty array
down_pos_flipped = [];

% Flip the downwards directed positively correlated averages
for i = 1:height(down_pos)
    down_pos_flipped(i) = 1 - down_pos(i);
end

% Transpose to be column vector
down_pos_flipped = down_pos_flipped';

% Define negative correlations in data
up_neg = NEGavgs(:,1);
down_neg = NEGavgs(:,2);

% Initialize empty array
up_neg_flipped = [];

% Flip the upwards directed negatively correlated averages
for i = 1:height(up_neg)
    up_neg_flipped(i) = 1 - up_neg(i);
end

% Transpose to be column vector
up_neg_flipped = up_neg_flipped';

% Get neighbors in matrices next to each other for each correlation,
% respectively
pos_corr_flips = horzcat(up_pos, down_pos_flipped);
neg_corr_flips = horzcat(up_neg_flipped, down_neg);

% Find neighbor averages for each correlation, respectively
pos_corr_flips_avg = mean(pos_corr_flips,2);
neg_corr_flips_avg = mean(neg_corr_flips,2);

% Plot results
figure
scatter(pos_corr_flips_avg, neg_corr_flips_avg, "filled")
xlim([0.5,1])
ylim([0.5,1])
xlabel('positive correlations')
ylabel('negative correlations')
set(gca,'linewidth',3,'fontSize',16)

% Compute stats for above analysis (p-value and confidence intervals)
[R,P,RL,RU] = corrcoef(pos_corr_flips_avg, neg_corr_flips_avg);


