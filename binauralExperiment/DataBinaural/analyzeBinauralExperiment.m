% Analyze binaural experiment

%% Loop through subjects and reorganize data

% Define number of subjects
nsubs = 10;

% Loop through subjects
for j = 1:nsubs
    load(['s00',num2str(j),'_toneData.mat']); % load 'rich' file
    load(['s00',num2str(j),'.mat']); % load 'diet' file

    DF = data; % rename data frames

    % Pull out "random" probe trials
    for k = 1:100
        randTrial_boolean(k) = toneStim(k).params.corrType=="binauralScintRand";
    end
   
    DF.choice(DF.choice==-1) = 0; % convert "falling" pitch response to 0 (instead of -1)

    POSavgs(j,:) = [mean(DF.choice(DF.corrparity==1 & DF.displacement==1 & randTrial_boolean'==0)) ... % +1 corr, +1 disp

        mean(DF.choice(DF.corrparity==1 & DF.displacement==-1 & randTrial_boolean'==0))]; % +1 corr, -1 disp

    NEGavgs(j,:) = [mean(DF.choice(DF.corrparity==-1 & DF.displacement==1 & randTrial_boolean'==0)) ... % -1 corr, +1 disp

        mean(DF.choice(DF.corrparity==-1 & DF.displacement==-1 & randTrial_boolean'==0))]; % -1 corr, -1 disp

    RANDavg(j) = mean(DF.choice(randTrial_boolean'==1)); % random trials

end

%% Calculate means and SEMs

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

%% Redo graphing

% Define number of conditions
conditions = [1:5];

% Concatenate data
data = horzcat(POSavgs, NEGavgs, RANDavg');

% Compute means
group_means = mean(data);

% Consolidate SEMs
sems = [sem1, sem2, sem3, sem4, sem5];

% Define number of conditions and subjects
num_conditions = size(data, 2);
num_subjects = size(data, 1);

% Plot results
figure;
bar_handle = bar(conditions, group_means, 'basevalue',0.5);
hold on;

errorbar(conditions, group_means, sems, 'LineStyle', 'none', 'color', [0 0 0])
hold on;

% Plot each subject with a dot and respective connecting lines
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

%% Perform stats: one sample t-tests to compare against chance (0.5)

% Reoganize data
up_pos = data(:,1);
down_pos = data(:,2);
up_neg = data(:,3);
down_neg = data(:,4);
random = data(:,5);

% Mean to test against
mu = 0.5;

% Perform the one-sample t-test
[h1, p1, ci1, stats1] = ttest(up_pos, mu);
[h2, p2, ci2, stats2] = ttest(down_pos, mu);
[h3, p3, ci3, stats3] = ttest(up_neg, mu);
[h4, p4, ci4, stats4] = ttest(down_neg, mu);
[h5, p5, ci5, stats5] = ttest(random, mu);

% Consolidate p-value results
results = [p1, p2, p3, p4, p5];