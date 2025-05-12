% Analyze coherence experiment

%% Loop through subjects and pull out data

numsubs = 10;
counter=1;
for subnum = 1:numsubs
    load(['s00', num2str(subnum),'.mat']);
    cohTranspose = data.coh';
    uC = unique(cohTranspose);
    disp(['unique coherences = ' num2str(uC')]);
    counter2 = 1;
    for C = [1,-1]
        for D = [1 -1]
            for ii=1:length(uC)
                ff = find( (cohTranspose == uC(ii)) .* (data.corrparity == C) .* (data.displacement == D));
                meanChoice(ii,counter,((C+1)+(D+1)/2)+1) = mean(data.choice(ff));
                semChoice(ii,counter,((C+1)+(D+1)/2)+1) = std(data.choice(ff))/sqrt(length(ff));
            end
            counter2 = counter2 + 1; 
        end
    end
       counter=counter+1;
end

%% Plot all coherence results on same graph

% Define coherences
cohs = [0 0.125 0.25 0.5 0.75 1]; % set x-axis values

% Setting correlation and direction key
% 1 = -C,-D (down neg) ; 2 = -C,+D (up neg); 3 = +C,-D (down pos); 4 =
% +C,+D (up pos)

% It uses this formula:
% corr=((C+1)+(D+1)/2)+1;

% Plot stim 4 (+C,+D)
corr = 4; 
mC2=(meanChoice(:,:,corr)+1)/2;
sem=std(mC2,[],2)./sqrt(numsubs);
y1 = mean(((meanChoice(:,:,corr)+1)/2),2) + sem; % above mean
y2 = mean(((meanChoice(:,:,corr)+1)/2),2) - sem; % below mean
figure;
plot(cohs(:),y1, 'w')
hold on;
plot(cohs(:),y2, 'w')
patch([cohs(:); flipud(cohs(:))],[y1(:); flipud(y2(:))], 'r', 'FaceAlpha',0.1, 'EdgeColor','none');
plot(cohs,mean(((meanChoice(:,:,corr)+1)/2),2), 'LineWidth', 2, 'color', 'r')
xticks(cohs);
hold on

% Plot stim 3 (+C,-D)
corr = 3; 
mC2=(meanChoice(:,:,corr)+1)/2;
sem=std(mC2,[],2)./sqrt(numsubs);
y1 = mean(((meanChoice(:,:,corr)+1)/2),2) + sem;
y2 = mean(((meanChoice(:,:,corr)+1)/2),2) - sem;
plot(cohs(:),y1, 'w')
hold on;
plot(cohs(:),y2, 'w')
patch([cohs(:); flipud(cohs(:))],[y1(:); flipud(y2(:))], 'b', 'FaceAlpha',0.1, 'EdgeColor','none');
plot(cohs,mean(((meanChoice(:,:,corr)+1)/2),2),'k--','LineWidth', 3, 'color', 'b')
xticks(cohs);
hold on

% Plot stim 2 (-C,+D)
corr = 2; 
mC2=(meanChoice(:,:,corr)+1)/2;
sem=std(mC2,[],2)./sqrt(numsubs);
y1 = mean(((meanChoice(:,:,corr)+1)/2),2) + sem; % above mean
y2 = mean(((meanChoice(:,:,corr)+1)/2),2) - sem; % below mean
plot(cohs(:),y1, 'w')
hold on;
plot(cohs(:),y2, 'w')
patch([cohs(:); flipud(cohs(:))],[y1(:); flipud(y2(:))], 'g', 'FaceAlpha',0.1, 'EdgeColor','none');
plot(cohs,mean(((meanChoice(:,:,corr)+1)/2),2), 'LineWidth', 3, 'color', 'g')
xticks(cohs);
hold on

% Plot stim 2 (-C,-D)
corr = 1; 
mC2=(meanChoice(:,:,corr)+1)/2;
sem=std(mC2,[],2)./sqrt(numsubs);
y1 = mean(((meanChoice(:,:,corr)+1)/2),2) + sem; % above mean
y2 = mean(((meanChoice(:,:,corr)+1)/2),2) - sem; % below mean
plot(cohs(:),y1, 'w')
hold on;
plot(cohs(:),y2, 'w')
patch([cohs(:); flipud(cohs(:))],[y1(:); flipud(y2(:))], 'm', 'FaceAlpha',0.1, 'EdgeColor','none');
plot(cohs,mean(((meanChoice(:,:,corr)+1)/2),2), 'k--','LineWidth', 3, 'color', 'm')
xticks(cohs);
hold off
ylim([0,1])
yline(0.5,'k--','linewidth',2);
xlabel('coherence');
ylabel('probability perceived rising');
set(gca,'linewidth',3,'fontSize',16)

%% Perform a two-way, repeated measures ANOVA for top 2 curves (up pos and down neg)

% First, we need to get the data in an appropriate structure

% There'll be 4 columns of the table: subject, condition, coherence, P(up)

% First column -- subjects:

% Define the sequence of subjects from 1 to 10
subs = 1:10;

% Repeat each subject entry
repeats = repelem(subs, 12);

% Transpose to column vector
subject_list = repeats';

% Second column -- conditions:

% There are 2 distinct conditions (up_pos, down_neg)
conditions = 1:2;

% Repeat each entry 6 times
repeats = repelem(conditions, 6);

% Transpose to column vector
condition_list_1 = repeats';

% Repeat that for each subject to make the 240x1 be 240x10
condition_list_2 = repmat(condition_list_1, length(subs), 1);

% Turn numbers into categorical variables
condition_list = cell(size(condition_list_2));

% Loop over each element of the column vector
for i = 1:length(condition_list_2)
    if condition_list_2(i) == 1
        condition_list{i} = 'up_pos';
    elseif condition_list_2(i) == 2
        condition_list{i} = 'down_neg';
    end
end

condition_list = categorical(condition_list);

% Third column -- coherences:

% Define the initial sequence
coherences = [0, 1/8, 1/4, 1/2, 3/4, 1]';

% Repeat the sequence 20 times to fill a 120x1 column vector
coherence_list = repmat(coherences, 20, 1);

% Fourth column -- P(up):

% Pull out and reshape rising positive and falling negative data
corr = 4; % rising positive (up pos)
up_pos_subs = ((meanChoice(:,:,corr)+1)/2);

corr = 1; % falling negative (down neg)
down_neg_subs = ((meanChoice(:,:,corr)+1)/2);

% Get each subject's rising positive data
sub1_up_pos = up_pos_subs(:,1);
sub2_up_pos = up_pos_subs(:,2);
sub3_up_pos = up_pos_subs(:,3);
sub4_up_pos = up_pos_subs(:,4);
sub5_up_pos = up_pos_subs(:,5);
sub6_up_pos = up_pos_subs(:,6);
sub7_up_pos = up_pos_subs(:,7);
sub8_up_pos = up_pos_subs(:,8);
sub9_up_pos = up_pos_subs(:,9);
sub10_up_pos = up_pos_subs(:,10);

% Get each subject's falling negative data
sub1_down_neg = down_neg_subs(:,1);
sub2_down_neg = down_neg_subs(:,2);
sub3_down_neg = down_neg_subs(:,3);
sub4_down_neg = down_neg_subs(:,4);
sub5_down_neg = down_neg_subs(:,5);
sub6_down_neg = down_neg_subs(:,6);
sub7_down_neg = down_neg_subs(:,7);
sub8_down_neg = down_neg_subs(:,8);
sub9_down_neg = down_neg_subs(:,9);
sub10_down_neg = down_neg_subs(:,10);

% Consolidate
response_list = vertcat(sub1_up_pos, sub1_down_neg, sub2_up_pos, sub2_down_neg, ...
    sub3_up_pos, sub3_down_neg, sub4_up_pos, sub4_down_neg, sub5_up_pos, sub5_down_neg, ...
    sub6_up_pos, sub6_down_neg, sub7_up_pos, sub7_down_neg, sub8_up_pos, sub8_down_neg, ...
    sub9_up_pos, sub9_down_neg, sub10_up_pos, sub10_down_neg);

% Make a table
data = table(subject_list, condition_list, coherence_list, response_list);

% Restructure everything 
sub1_responses = horzcat(sub1_up_pos', sub1_down_neg');
sub2_responses = horzcat(sub2_up_pos', sub2_down_neg');
sub3_responses = horzcat(sub3_up_pos', sub3_down_neg');
sub4_responses = horzcat(sub4_up_pos', sub4_down_neg');
sub5_responses = horzcat(sub5_up_pos', sub5_down_neg');
sub6_responses = horzcat(sub6_up_pos', sub6_down_neg');
sub7_responses = horzcat(sub7_up_pos', sub7_down_neg');
sub8_responses = horzcat(sub8_up_pos', sub8_down_neg');
sub9_responses = horzcat(sub9_up_pos', sub9_down_neg');
sub10_responses = horzcat(sub10_up_pos', sub10_down_neg');

new_structure = vertcat(sub1_responses, sub2_responses, sub3_responses, sub4_responses, ...
    sub5_responses, sub6_responses, sub7_responses, sub8_responses, sub9_responses, sub10_responses);

new_table = array2table(new_structure);
new_table.Properties.VariableNames = {'up_pos1', 'up_pos2', 'up_pos3', 'up_pos4', 'up_pos5', 'up_pos6', 'down_neg1', 'down_neg2', 'down_neg3','down_neg4', 'down_neg5', 'down_neg6'};

% Create the within-subjects design
withinDesign = table([1 1 1 1 1 1 2 2 2 2 2 2]',[1 2 3 4 5 6 1 2 3 4 5 6]','VariableNames',{'Conditions','Coherences'});
withinDesign.Conditions = categorical(withinDesign.Conditions);
withinDesign.Coherences = categorical(withinDesign.Coherences);

% Create repeated measures model
rm = fitrm(new_table, 'up_pos1-down_neg6 ~ 1', 'WithinDesign', withinDesign);

% Perform anova (remove semicolon to view table generated by ranova)
result = ranova(rm, 'WithinModel', 'Conditions*Coherences');

%% Perform a two-way, repeated measures ANOVA for bottom 2 curves (down pos, up neg)

% First column -- subjects

% Define the sequence of subjects from 1 to 10
subs = 1:10;

% Repeat each subject entry
repeats = repelem(subs, 12);

% Transpose to column vector
subject_list = repeats';

% Second column -- conditions:

% There are 2 distinct conditions (down_pos, up_neg)
conditions = 1:2;

% Repeat each entry 6 times
repeats = repelem(conditions, 6);

% Transpose to column vector
condition_list_1 = repeats';

% Repeat that for each subject to make the 240x1 be 240x10
condition_list_2 = repmat(condition_list_1, length(subs), 1);

% Turn numbers into categorical variables
condition_list = cell(size(condition_list_2));

% Loop over each element of the column vector
for i = 1:length(condition_list_2)
    if condition_list_2(i) == 1
        condition_list{i} = 'down_pos';
    elseif condition_list_2(i) == 2
        condition_list{i} = 'up_neg';
    end
end

condition_list = categorical(condition_list);

% Third column -- coherences:

% Define the initial sequence
coherences = [0, 1/8, 1/4, 1/2, 3/4, 1]';

% Repeat the sequence 20 times to fill a 120x1 column vector
coherence_list = repmat(coherences, 20, 1);

% Fourth column -- P(up):

% Pull out and reshape falling positive and rising negative data
corr = 3; % falling positive
down_pos_subs = ((meanChoice(:,:,corr)+1)/2);

corr = 2; % rising negative
up_neg_subs = ((meanChoice(:,:,corr)+1)/2);

% Get each subject's down positive data
sub1_down_pos = down_pos_subs(:,1);
sub2_down_pos = down_pos_subs(:,2);
sub3_down_pos = down_pos_subs(:,3);
sub4_down_pos = down_pos_subs(:,4);
sub5_down_pos = down_pos_subs(:,5);
sub6_down_pos = down_pos_subs(:,6);
sub7_down_pos = down_pos_subs(:,7);
sub8_down_pos = down_pos_subs(:,8);
sub9_down_pos = down_pos_subs(:,9);
sub10_down_pos = down_pos_subs(:,10);

% Get each subject's up negative data
sub1_up_neg = up_neg_subs(:,1);
sub2_up_neg = up_neg_subs(:,2);
sub3_up_neg = up_neg_subs(:,3);
sub4_up_neg = up_neg_subs(:,4);
sub5_up_neg = up_neg_subs(:,5);
sub6_up_neg = up_neg_subs(:,6);
sub7_up_neg = up_neg_subs(:,7);
sub8_up_neg = up_neg_subs(:,8);
sub9_up_neg = up_neg_subs(:,9);
sub10_up_neg = up_neg_subs(:,10);

response_list = vertcat(sub1_down_pos, sub1_up_neg, sub2_down_pos, sub2_up_neg, ...
    sub3_down_pos, sub3_up_neg, sub4_down_pos, sub4_up_neg, sub5_down_pos, sub5_up_neg, ...
    sub6_down_pos, sub6_up_neg, sub7_down_pos, sub7_up_neg, sub8_down_pos, sub8_up_neg, ...
    sub9_down_pos, sub9_up_neg, sub10_down_pos, sub10_up_neg);

% Make a table
data = table(subject_list, condition_list, coherence_list, response_list);

% Restructure everything 
sub1_responses = horzcat(sub1_down_pos', sub1_up_neg');
sub2_responses = horzcat(sub2_down_pos', sub2_up_neg');
sub3_responses = horzcat(sub3_down_pos', sub3_up_neg');
sub4_responses = horzcat(sub4_down_pos', sub4_up_neg');
sub5_responses = horzcat(sub5_down_pos', sub5_up_neg');
sub6_responses = horzcat(sub6_down_pos', sub6_up_neg');
sub7_responses = horzcat(sub7_down_pos', sub7_up_neg');
sub8_responses = horzcat(sub8_down_pos', sub8_up_neg');
sub9_responses = horzcat(sub9_down_pos', sub9_up_neg');
sub10_responses = horzcat(sub10_down_pos', sub10_up_neg');

new_structure = vertcat(sub1_responses, sub2_responses, sub3_responses, sub4_responses, ...
    sub5_responses, sub6_responses, sub7_responses, sub8_responses, sub9_responses, sub10_responses);

new_table = array2table(new_structure);
new_table.Properties.VariableNames = {'down_pos1', 'down_pos2', 'down_pos3', 'down_pos4', 'down_pos5', 'down_pos6', 'up_neg1', 'up_neg2', 'up_neg3','up_neg4', 'up_neg5', 'up_neg6'};

% Create the within-subjects design
withinDesign = table([1 1 1 1 1 1 2 2 2 2 2 2]',[1 2 3 4 5 6 1 2 3 4 5 6]','VariableNames',{'Conditions','Coherences'});
withinDesign.Conditions = categorical(withinDesign.Conditions);
withinDesign.Coherences = categorical(withinDesign.Coherences);

% Create repeated measures model
rm = fitrm(new_table, 'down_pos1-up_neg6 ~ 1', 'WithinDesign', withinDesign);

% Perform anova (remove semicolon to view table generated by ranova)
result = ranova(rm, 'WithinModel', 'Conditions*Coherences');