% Analyze glider data

%% Loop through subjects and pull out data

numsubs = 10;
counter=1; % count is a random made up thing that helps exclude a subject without their data column getting turned into all zeros
for subnum = (1:numsubs)
    load(['subj', num2str(subnum),'.mat']);
    uG = unique(data.corrTypeNumber);
    disp(['unique types = ' num2str(uG')]);
    counter2 = 1;
    for C = [1 -1]
        for D = [1 -1]
            for ii=1:length(uG)
                ff = find( (data.corrTypeNumber == uG(ii)) .* (data.corrParity == C) .* (data.direction == D));
                meanChoice(ii,counter,((C+1)+(D+1)/2)+1) = mean(data.choice(ff));
                semChoice(ii,counter,((C+1)+(D+1)/2)+1) = std(data.choice(ff))/sqrt(length(ff));
            end
            counter2 = counter2 + 1; 
        end
    end
       counter=counter+1;
end

%% Plot results

% Load responses to each condition
% First column is converging, second column is diverging
load('up_pos.mat')
load('up_neg.mat')
load('down_pos.mat')
load('down_neg.mat')

% Make data into organized variables within the workspace
up_div_pos = up_pos(:,2);
up_div_neg = up_neg(:,2);
down_div_pos = down_pos(:,2);
down_div_neg = down_neg(:,2);
up_conv_pos = up_pos(:,1);
up_conv_neg = up_neg(:,1);
down_conv_pos = down_pos(:,1);
down_conv_neg = down_neg(:,1);

% Load SEMs of each condition
load('sem_up_pos.mat')
sem_up_pos = sem';
clear sem;

load('sem_up_neg.mat')
sem_up_neg = sem'; 
clear sem;

load('sem_down_pos.mat')
sem_down_pos = sem';
clear sem;

load('sem_down_neg.mat')
sem_down_neg = sem';
clear sem;

% Concatenate
data = horzcat(up_div_pos, up_div_neg, down_div_pos, down_div_neg, up_conv_pos, up_conv_neg, down_conv_pos, down_conv_neg);

% Make SEMs into organized variables within the workspace
sem_up_div_pos = sem_up_pos(:,2);
sem_up_div_neg = sem_up_neg(:,2);
sem_down_div_pos = sem_down_pos(:,2);
sem_down_div_neg = sem_down_neg(:,2);
sem_up_conv_pos = sem_up_pos(:,1);
sem_up_conv_neg = sem_up_neg(:,1);
sem_down_conv_pos = sem_down_pos(:,1);
sem_down_conv_neg = sem_down_neg(:,1);

% Concatenate SEMs
sems = horzcat(sem_up_div_pos, sem_up_div_neg, sem_down_div_pos, sem_down_div_neg, sem_up_conv_pos, sem_up_conv_neg, sem_down_conv_pos, sem_down_conv_neg);

% Organize everything for graphing
conditions = [1:8];
group_means = mean(data);
num_conditions = size(data, 2);
num_subjects = size(data, 1);

% Plot results
figure;
bar_handle = bar(conditions, group_means, 'basevalue',0.5);
hold on;

errorbar(conditions, group_means, sems, 'LineStyle', 'none', 'color', [0 0 0])
hold on;

% Plot individual subjects as dots with connected lines
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
ylim([0 1])
xticklabels({'up div pos', 'up div neg', 'down div pos', 'down div neg', 'up conv pos', 'up conv neg', 'down conv pos', 'down conv neg'})

%% Replot net graph --> CLEAR BEFORE RUNNING!

% Load responses to each condition
% First column is converging, second column is diverging
load('up_pos.mat')
load('up_neg.mat')
load('down_pos.mat')
load('down_neg.mat')

% Make them into organized variables within the workspace
up_div_pos = up_pos(:,2);
up_div_neg = up_neg(:,2);
down_div_pos = down_pos(:,2);
down_div_neg = down_neg(:,2);
up_conv_pos = up_pos(:,1);
up_conv_neg = up_neg(:,1);
down_conv_pos = down_pos(:,1);
down_conv_neg = down_neg(:,1);

% Calculate net values
pos_div = up_div_pos - down_div_pos;
neg_div = up_div_neg - down_div_neg;
pos_conv = up_conv_pos - down_conv_pos;
neg_conv = up_conv_neg - down_conv_neg;

% Concatenate into a dataframe
data = horzcat(pos_div, neg_div, pos_conv, neg_conv);

% Compute means
group_means = mean(data);

% Define conditions and sizes
conditions = 1:4;
num_conditions = size(data, 2);
num_subjects = size(data, 1);

% Calculate SEMs
sem1 = std(pos_div)./sqrt(num_subjects);
sem2 = std(neg_div)./sqrt(num_subjects);
sem3 = std(pos_conv)./sqrt(num_subjects);
sem4 = std(neg_conv)./sqrt(num_subjects);
sems = [sem1, sem2, sem3, sem4];

% Plot results
figure;
bar_handle = bar(conditions, group_means, 'basevalue',0);
hold on;

errorbar(conditions, group_means, sems, 'LineStyle', 'none', 'color', [0 0 0])
hold on;

% Plot individual subjects as dots with connected lines
for i = 1:num_subjects
    x = 1:num_conditions;
    y = data(i, :);
    scatter(x, y, 'k', 'filled')
    plot(x, y, '-o', 'color', [0 0 0])
end
hold off;

set(gca,'linewidth',3,'fontSize',16)
xlabel('stimulus type')
ylabel('net probability perceived rising')
xticklabels({'pos div', 'neg div', 'pos conv', 'neg conv'})

%% Perform stats: two-sample t-tests

% Compute stats
[h1,p1,ci1,stats1] = ttest2(pos_div,neg_div);
[h2,p2,ci2,stats2] = ttest2(pos_conv,neg_conv);

% Save results
results = [p1, p2];