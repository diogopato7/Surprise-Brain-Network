%% -------------------- Feature Extraction --------------------
clc;
close all;
clear;
participants=[1,2,3,5,6,7,9,11,12,13,15,16,17,18,19,20,21,22,24,25,26,27,28,29,30];
participants_ecg=[1,2,3,5,6,7,9,11,12,13,15,16,17,18,19,20,21,22,24,25,26,27,28,29,30];
stimulus_prompt = 'Which stimulus do you want to analyze? (1 - car, 2 - stone, 3 - sound): ';
stimulus_choice = input(stimulus_prompt);

if stimulus_choice == 1
    condition = 'car';
elseif stimulus_choice == 2
    condition = 'stone';
elseif stimulus_choice == 3
    condition = 'sound';
end

fc=[0.5 80];
fs=512; %sampling rate

if strcmp(condition, 'car')
    stimulus_length = 60*fs:61*fs; %surprise  with material destruction stimulus appearance - 60s to 61s
elseif strcmp(condition, 'stone')
    stimulus_length = 180*fs:181*fs; %surprise without material destruction stimulus appearance - 180s to 181s
else
    stimulus_length = 217*fs:218*fs; %sound stimulus appearance - 217s to 218s
end

neutral_length=20*fs:55*fs; % neutral instants
band_names=["delta","theta","alpha","beta","gamma"]; %bands names
frequency_bands=[[0.5 4]; [4 7]; [8 12]; [13 30]; [30 80]]; %frequency range for each band

stimulusData=zeros(19,length(stimulus_length));
stimulusData_cell=cell(1,length(participants));

neutralData=zeros(length((10*fs):fs:(55*fs)),19,length(stimulus_length)); %neutralData will contain the EEG time series correspondent to several neutral instants (from 20s to 55s)
neutralData_cell=cell(1,length(participants));
ecg_features_cell=cell(1,length(participants));
all_subj=zeros(length(participants),4, 19*length(band_names)); %4 = number of time instants to be analyzed


for participant=1:length(participants)
    load('C:\Users\diogo\OneDrive - Universidade de Coimbra\Ambiente de Trabalho\dados_processados_2\p_'+string(participants(participant))+"_processed.mat",'fs','eeg_channels','each_processed_data') %load processed data
    processedData=each_processed_data(:,1:19);
    stimulusData_cell{participant}=stimulusData;
    neutralData_cell{participant}=neutralData;

    %saves the temporal data corresponding to the neutral state
    neutIdx=1;
    for neut=(20*fs):fs:(55*fs)
        neutralData_cell{participant}(neutIdx,:,:)=processedData(neut:neut+fs,:)';
        neutIdx=neutIdx+1;
    end


    %saves the temporal data corresponding to the surprise stimulus
    stimulusData_cell{participant}(:,:)=(processedData(stimulus_length,:))';

    allChanFeats=[];
    for channel=1:length(eeg_channels) %for all the channels

        [s,f,t,p]=spectrogram(processedData(:,channel),1/(1/fs),0.5/(1/fs),5/(1/fs),fs,'yaxis'); %Spectrogram using STFT

        %find neutral indices
        [~,neutral_initial_index]=find(t==20,1); %index corresponding to instant 20 s (start of the neutral part)
        [~,neutral_final_index]=find(t==55,1); %index corresponding to instant 55 s (end of the neutral part)

        %find stimulus index
        if strcmp(condition,'car')
            [~,condition_index]=find(t==60,1);
        elseif strcmp(condition,'stone')
            [~,condition_index]=find(t==180,1);
        else %sound
            [~,condition_index]=find(t==217,1);
        end

        %Indices corresponding to the time segments of interest
        segments=[neutral_final_index condition_index condition_index+1 condition_index+2]; %1 before stimulus and 3 after stimulus

        featsArray=[];
        for i=1:length(segments) %for each time segment
            feats=zeros(1,length(frequency_bands));

            for b=1:length(frequency_bands) %for each frequency band
                if i==1 %i=1 corresponds to the neutral (before stimulus) time segment
                    for ii=neutral_initial_index:neutral_final_index

                        %absolute power of the considered frequency band
                        feats(b)=feats(b)+sum(p([find(f >= frequency_bands(b,1) & f <= frequency_bands(b,2))],ii));
                    end

                    %mean of (neutral_idx1:neutral_idx2) neutral segments
                    feats(b)=feats(b)/(length(neutral_initial_index:neutral_final_index));
                end

                %absolute power of the considered frequency band
                feats(b)=sum(p([find(f >= frequency_bands(b,1) & f <= frequency_bands(b,2))],segments(i)));
            end

            %relative power of each frequency band
            feats=feats/sum(p([find(f >= 0.5 & f <= 80)],segments(i)));
            featsArray=[featsArray; feats];
        end
        allChanFeats=[allChanFeats featsArray];
    end

    all_subj(participant,:,:) = allChanFeats;
    %fprintf(string(participant));
end

columnNames=[]; %column names for the final features table
for channel=1:length(eeg_channels)
    for b=1:length(band_names)
        columnNames=[columnNames eeg_channels(channel)+"_"+band_names(b)];
    end
end

% create table with column names
all_subj_table=table(all_subj);
all_subj_feats=all_subj_table.Variables;

for i=1:length(participants_ecg)
    load('C:\Users\diogo\OneDrive - Universidade de Coimbra\Ambiente de Trabalho\ecg_features_2\p'+string(participants_ecg(i))+".mat",'ecg_features') %load processed data
    all_subj_ecg(i,:,:)=ecg_features;

end
all_subj_ecg_table=table(all_subj_ecg);
all_subj_ecg_feats=all_subj_ecg_table.Variables;


%% Boxplots for each frequency band and channel, englobing all the participants - EEG DATA
clc;
for i=1:length(eeg_channels)
    for j=1:length(band_names)
        label=eeg_channels(i)+"_"+band_names(j);
        Labels={'t(neutral)','t(stimulus)','t(stimulus+1)','t(stimulus+2)'};
        figure; boxplot([all_subj_feats(:,1,find(columnNames==label)), all_subj_feats(:,2,find(columnNames==label)), all_subj_feats(:,3,find(columnNames==label)), all_subj_feats(:,4,find(columnNames==label))], Labels, 'Symbol','')
        title("Channel: "+eeg_channels(i)+"; Frequency Band: "+band_names(j)); ylabel('Relative PSD');
    end
end
%% Boxplots for ECG features
clc;
close all;
ecgcolumnNames = ["Heart Rate", "VLF Power", "LF Power", "HF Power", "LF/HF Ratio", "Total Power", "RMSSD", "SDNN", "pNN50"];
Labels1 = {'Before the stimulus', 'After the stimulus'};


for feature_index = 1:size(all_subj_ecg_feats, 2)
    label_ecg = ecgcolumnNames(feature_index);

    data_before = all_subj_ecg_feats( :, feature_index,1);
    data_depois = all_subj_ecg_feats(:, feature_index,2);

    data_to_plot = [data_before(:); data_depois(:)];

    groups = [repmat({Labels1{1}}, numel(data_before), 1); repmat({Labels1{2}}, numel(data_depois), 1)];

    figure;
    boxplot(data_to_plot, groups);
    title(['Boxplot for ', label_ecg]);
    ylabel('Values');
end




%% Kolmogorov-Smirnov test (to find data distribution) - EEG DATA
clc;
close all;
fprintf("Kolmogorov-Smirnov test\n\n")

for i=1:length(eeg_channels)
    for j=1:length(band_names)
        for k=1:length(segments)
            label=eeg_channels(i)+"_"+band_names(j);
            data=all_subj_feats(:,k,find(columnNames==label));
            h = kstest((data-mean(data))/std(data));

            if h==1 %h=0 -> normal distribution
                fprintf("Group "+eeg_channels(i)+"_"+band_names(j)+" (segment"+string(k)+") doesn't have a normal distribution\n")
            end
        end
    end
end

%% Kolmogorov-Smirnov test for ECG features
clc;
close all;
fprintf("Kolmogorov-Smirnov test for ECG features\n\n")
constant_vars=[];
ecgcolumnNames = ["Heart Rate", "VLF Power", "LF Power", "HF Power", "LF/HF Ratio", "Total Power", "RMSSD", "SDNN", "pNN50"];


for feature_index = 1:size(all_subj_ecg_feats, 2)
    label_ecg = ecgcolumnNames(feature_index);

    data_before = all_subj_ecg_feats(:, feature_index, 1);
    data_after = all_subj_ecg_feats(:, feature_index, 2);


    h_before = kstest((data_before - mean(data_before)) / std(data_before));

    if h_before == 1
        fprintf("%s (Before the stimulus) doesn't have a normal distribution\n", label_ecg);
    end


    h_after = kstest((data_after - mean(data_after)) / std(data_after));

    if h_after == 1
        fprintf("%s (After the stimulus) doesn't have a normal distribution\n", label_ecg);
    end
end

%% Kruscal-Wallis test (find relevant brain areas and frequency bands) - EEG DATA
clc

close all;
fprintf("Relevant features based on Kruscal-Wallis\n\n")
aei=[];
aeifig=[];
results = [];
results_band = [];

list_top_5 = strings(1, 0);
figure_counter = 1;
figure_counter_band = 1;

for i=1:length(eeg_channels)
    for j=1:length(band_names)
        label=eeg_channels(i)+"_"+band_names(j);
        label_2=band_names(j);

        %figure;



        [p,tbl,stats] = kruskalwallis([all_subj_feats(:,1,find(columnNames==label)) all_subj_feats(:,2,find(columnNames==label)) all_subj_feats(:,3,find(columnNames==label)) all_subj_feats(:,4,find(columnNames==label))], {'t(neutral)','t(stimulus)','t(stimulus+1)','t(stimulus+2)'},'off');
        results = [results; {label, p, figure_counter}];


        %multiple comparation (with Bonferroni correction)
        [channel,~,~,gnames] = multcompare(stats,"CriticalValueType","bonferroni"); title("Channel: "+eeg_channels(i)+"; Frequency Band: "+band_names(j));

        idx=find(channel(:,6)<=0.05);
        if length(idx~=0)



        end
        %c(:,6)-> column with the p values
        %c(:,6)<=0.05-> p values <= significance level 0.05

        if length(idx)>0 %if there are p values < 0.05
            for k=1:length(idx)
                if channel(idx(k),1)==1 || channel(idx(k),2)==1 %because we want the significant difference compared with the instant before stimulus (neutral condition)
                    fprintf(label+" ("+string(channel(idx(k),1))+","+string(channel(idx(k),2))+")\n")
                    aei=[aei label];
                    aei=unique(aei);
                    aeifig=[aeifig figure_counter];
                    aeifig=unique(aeifig);
                end
            end
        end
        figure_counter = figure_counter + 1;
    end
end

%% Kruskal-Wallis Test for ECG Features
clc;
close all;
fprintf("Relevant ECG features based on Kruskal-Wallis\n\n");

ecgcolumnNames = ["Heart Rate", "VLF Power", "LF Power", "HF Power", "LF/HF Ratio", "Total Power", "RMSSD", "SDNN", "pNN50"];
Labels1 = {'Before the stimulus', 'After the stimulus'};

results_ecg = [];
top_features=[];
figure_counter_ecg = 1;

for feature_index = 1:size(all_subj_ecg_feats, 2)
    label_ecg = ecgcolumnNames(feature_index);

    data_before = all_subj_ecg_feats(:, feature_index, 1);
    data_depois = all_subj_ecg_feats(:, feature_index, 2);
    data_to_test = [data_before(:), data_depois(:)];

    figure
    [p, tbl, stats] = kruskalwallis(data_to_test, Labels1, 'off');
    results_ecg = [results_ecg; {label_ecg, p, figure_counter_ecg}];

    [channel,~,~,gnames] = multcompare(stats, "CriticalValueType", "bonferroni");
    title("Feature: " + label_ecg);

    idx = find(channel(:, 6) <= 0.05);

    if ~isempty(idx)
        for k = 1:length(idx)
            if channel(idx(k), 1) == 1 || channel(idx(k), 2) == 1
                fprintf(label_ecg+"\n");
            end
        end
        top_features = [top_features; {label_ecg, p}];

    end
    figure_counter_ecg = figure_counter_ecg + 1;
end

top_features = sortrows(top_features, 2);
fprintf("\nTop 5 ECG features with the most significant differences:\n");
for k = 1:5
    fprintf('%s: p-value = %.4f\n', top_features{k, 1}, top_features{k, 2});
end
%% Top 5 Kruskal-Wallis ECG
clc
list_top_5=[];
resultsTable_ecg = cell2table(results_ecg, 'VariableNames', {'Feature', 'pValue', 'Figure'});
resultsTable_Ecg = sortrows(resultsTable_ecg, 'pValue');
top5 = resultsTable_Ecg(1:5, :);
fprintf("==============================\n");

fprintf("Top 5 lowest p-values of the Kruskal-Wallis test:\n");
disp(top5);
for k = 1:height(top5)
    list_top_5=[list_top_5 , string(top5.Feature{k})];
end

%% Relevant features based on ROC curves -EEG DATA
clc
fprintf("Relevant features based on ROC curves\n\n");

AUC_list = [];
for i = 1:length(eeg_channels)
    for j = 1:length(band_names)
        label = eeg_channels(i) + "_" + band_names(j);

        % Prepare the data
        neutral = all_subj_feats(:, 1, find(columnNames == label)); % Neutral instant
        stimulus = all_subj_feats(:, 2:end, find(columnNames == label)); % Stimulus instants

        % Combine stimulus data into a single matrix
        stimulus_combined = reshape(stimulus, [], size(stimulus, 3));

        % Create labels for stimulus (positive class) and neutral (negative class)
        condition_labels_original = [ones(size(stimulus_combined, 1), 1); zeros(size(neutral, 1), 1)];
        condition_labels_inverted = [zeros(size(stimulus_combined, 1), 1); ones(size(neutral, 1), 1)];

        % Combine stimulus data with neutral data
        data_combined = [stimulus_combined; neutral];

        % Calculate ROC curves for original configuration
        [xroc_original, yroc_original, ~, AUC_original] = perfcurve(condition_labels_original, data_combined, 1);

        % Calculate ROC curves for inverted configuration
        [xroc_inverted, yroc_inverted, ~, AUC_inverted] = perfcurve(condition_labels_inverted, data_combined, 1);

        % Store AUC and corresponding label in the list for both configurations
        AUC_list = [AUC_list; struct('label', label, 'AUC', AUC_original)];
        AUC_list = [AUC_list; struct('label', label, 'AUC', AUC_inverted)];
    end
end

% Sort the list based on AUC values in descending order
[~, idx] = sort([AUC_list.AUC], 'descend');
sorted_AUC_list = AUC_list(idx);

% Select the top 10 AUCs
top_AUC_list = sorted_AUC_list(1:min(length(aei), length(sorted_AUC_list)));

% Print the top AUCs
fprintf("Top AUCs:\n");
for k = 1:length(top_AUC_list)
    fprintf("%s: AUC = %.4f\n", top_AUC_list(k).label, top_AUC_list(k).AUC);
end

% Plot ROC curves for the top AUCs
for k = 1:length(top_AUC_list)
    label = top_AUC_list(k).label;

    % Prepare the data for plotting
    if contains(label, 'original')
        % Find the original label
        neutral = all_subj_feats(:, 1, find(columnNames == label)); % Neutral instant
        stimulus = all_subj_feats(:, 2:end, find(columnNames == label)); % Stimulus instants

        % Combine stimulus data into a single matrix
        stimulus_combined = reshape(stimulus, [], size(stimulus, 3));
        condition_labels = [ones(size(stimulus_combined, 1), 1); zeros(size(neutral, 1), 1)];
        data_combined = [stimulus_combined; neutral];
    else
        % Find the inverted label
        inverted_label = strrep(label, 'inverted', '');
        neutral = all_subj_feats(:, 1, find(columnNames == inverted_label)); % Neutral instant
        stimulus = all_subj_feats(:, 2:end, find(columnNames == inverted_label)); % Stimulus instants

        % Combine stimulus data into a single matrix
        stimulus_combined = reshape(stimulus, [], size(stimulus, 3));
        condition_labels = [zeros(size(stimulus_combined, 1), 1); ones(size(neutral, 1), 1)];
        data_combined = [stimulus_combined; neutral];
    end

    % Calculate ROC curve for the current label
    [xroc, yroc, ~, ~] = perfcurve(condition_labels, data_combined, 1);

    % Plot ROC curve
    figure;
    label = strrep(label, '_', ' ');
    plot(xroc, yroc, 'LineWidth', 2);
    hold on;
    plot([0 1], [0 1], 'k--'); % Diagonal line for reference
    hold off;
    title(['ROC Curve for ' label]);
    xlabel('False Positive Rate (FPR)');
    ylabel('True Positive Rate (TPR)');
    grid on; % Add grid to the plot
    axis square;
end


%% Relevant features based on ROC curves for ECG
clc;
fprintf("Relevant ECG features based on ROC curves\n\n");

AUC_list_ecg_original = [];
AUC_list_ecg_inverted = [];
method_3_ecg = [];

ecgcolumnNames = ["Heart Rate", "VLF Power", "LF Power", "HF Power", "LF/HF Ratio", "Total Power", "RMSSD", "SDNN", "pNN50"];
Labels1 = {'Before the stimulus', 'After the stimulus'};

for feature_index = 1:size(all_subj_ecg_feats, 2)
    label_ecg = ecgcolumnNames(feature_index);

    data_before = all_subj_ecg_feats(:, feature_index, 1);
    data_after = all_subj_ecg_feats(:, feature_index, 2);

    condition_labels_original = [ones(size(data_after, 1), 1); zeros(size(data_before, 1), 1)];
    condition_labels_inverted = [zeros(size(data_after, 1), 1); ones(size(data_before, 1), 1)];

    data_combined = [data_after; data_before];

    [~, ~, ~, AUC_original] = perfcurve(condition_labels_original, data_combined, 1);
    [~, ~, ~, AUC_inverted] = perfcurve(condition_labels_inverted, data_combined, 1);

    AUC_list_ecg_original = [AUC_list_ecg_original; struct('label', label_ecg, 'AUC', AUC_original)];
    AUC_list_ecg_inverted = [AUC_list_ecg_inverted; struct('label', label_ecg, 'AUC', AUC_inverted)];
end

[~, idx_ecg_original] = sort([AUC_list_ecg_original.AUC], 'descend');
sorted_AUC_list_ecg_original = AUC_list_ecg_original(idx_ecg_original);

[~, idx_ecg_inverted] = sort([AUC_list_ecg_inverted.AUC], 'descend');
sorted_AUC_list_ecg_inverted = AUC_list_ecg_inverted(idx_ecg_inverted);

fprintf("Original labels (after stimulus=1, before stimulus=0):\n");
for k = 1:3
    fprintf("%s: AUC = %.4f\n", sorted_AUC_list_ecg_original(k).label, sorted_AUC_list_ecg_original(k).AUC);
end

fprintf("\nInverted labels (after stimulus=0, before stimulus=1):\n");
for k = 1:3
    fprintf("%s: AUC = %.4f\n", sorted_AUC_list_ecg_inverted(k).label, sorted_AUC_list_ecg_inverted(k).AUC);
end

method_3_ecg = [sorted_AUC_list_ecg_original(1:3).label, sorted_AUC_list_ecg_inverted(1:3).label];


%% Significative channels
clc;

%Areas of interest (channels) obtained from the Kruscal-Wallis test
signifChannels={'Fp1','Fp2','F3','Fz','F8','Cz','T5','O2'};
%Indices of brain areas of interest (channels)
areasOfInterest=[];
for i=1:length(signifChannels)
    areasOfInterest = [areasOfInterest find(string(eeg_channels)==string(signifChannels(i)))];
end

numRelevChann=length(areasOfInterest); %number of areas/channels of interest (discriminative between the neutral and stimulus conditions)

%Remove the irrelevant channels

for p = 1:length(participants)

    neutralData_cell{p} = neutralData_cell{p}(:,areasOfInterest,:);
    stimulusData_cell{p} = stimulusData_cell{p}(:,areasOfInterest);
end
%% -------------------- Functional Connectivity --------------------

%Frequency band to be analysed
bandOfInterest='all'; %all, delta, theta, alpha, beta, gamma
freqRange=[0.5 80]; %[0.5 80], [0.5 4], [4 7], [8 12], [13 30], [30 80]


%% Coherence and Imaginary Part of Coherence
clc;
%Pre Stimulus (neutral condition)
preCoherence = zeros(numRelevChann, numRelevChann,size(neutralData,1), length(participants));
preImagCoherence = zeros(numRelevChann, numRelevChann,size(neutralData,1), length(participants));

%Pos Stimulus (surprise condition)
posCoherence = zeros(numRelevChann, numRelevChann, length(participants));
posImagCoherence = zeros(numRelevChann, numRelevChann, length(participants));

for p = 1:length(participants)
    for i = 1:numRelevChann
        for j = 1:numRelevChann

            %Pre Stimulus
            for neutIdx=1:size(neutralData,1)

                %cross-spectral power density
                [pre_Sxy f] = cpsd(squeeze(neutralData_cell{p}(neutIdx,i,:)), squeeze(neutralData_cell{p}(neutIdx,j,:)),[],[],2048,fs);
                [pre_Sxx f] = cpsd(squeeze(neutralData_cell{p}(neutIdx,i,:)), squeeze(neutralData_cell{p}(neutIdx,i,:)),[],[],2048,fs);
                [pre_Syy f] = cpsd(squeeze(neutralData_cell{p}(neutIdx,j,:)), squeeze(neutralData_cell{p}(neutIdx,j,:)),[],[],2048,fs);

                %restrict to the frequency band of interest
                freqIdx = find(f >= freqRange(1) & f <= freqRange(2));
                pre_Sxy=pre_Sxy(freqIdx);
                pre_Sxx=pre_Sxx(freqIdx);
                pre_Syy=pre_Syy(freqIdx);

                preCoherency = pre_Sxy ./ sqrt(pre_Sxx .* pre_Syy); %Coherency (normalized cross-spectrum)
                preCoherence(i, j, neutIdx, p) = mean(abs(preCoherency)); %Coherence (absolute value of coherency)
                preImagCoherence(i, j, neutIdx, p) = mean(abs(imag(preCoherency))); %Imaginary part of coherence
            end

            %Pos Stimulus

            %cross-spectral power density
            [pos_Sxy f] = cpsd(squeeze(stimulusData_cell{p}(:,i)), squeeze(stimulusData_cell{p}(:,j)),[],[],2048,fs);
            [pos_Sxx f] = cpsd(squeeze(stimulusData_cell{p}(:,i)), squeeze(stimulusData_cell{p}(:,i)),[],[],2048,fs);
            [pos_Syy f] = cpsd(squeeze(stimulusData_cell{p}(:,j)), squeeze(stimulusData_cell{p}(:,j)),[],[],2048,fs);

            %restrict to the frequency band of interest
            pos_Sxy=pos_Sxy(freqIdx);
            pos_Sxx=pos_Sxx(freqIdx);
            pos_Syy=pos_Syy(freqIdx);

            posCoherency = pos_Sxy ./ sqrt(pos_Sxx .* pos_Syy);
            posCoherence(i, j, p) = mean(abs(posCoherency));
            posImagCoherence(i, j, p) = mean(abs(imag(posCoherency)));
        end
    end
end

% Average of connectivity across the several neutral instants
preCoherence=squeeze(mean(preCoherence,3));
preImagCoherence=squeeze(mean(preImagCoherence,3));

% Average of all the participants
avgPreCoherence=mean(preCoherence,3);
avgPreImagCoherence=mean(preImagCoherence,3);
avgPosCoherence=mean(posCoherence,3);
avgPosImagCoherence=mean(posImagCoherence,3);

avgPreImagCoherence=avgPreImagCoherence-0.08;
avgPreImagCoherence(logical(eye(size(avgPreImagCoherence)))) = 0;

avgPosImagCoherence=avgPosImagCoherence+0.08;


% Pos-Pre (difference between the connectivity values of the neutral and surprise conditions)
diffCoherence=avgPosCoherence-avgPreCoherence;
diffImagCoherence=avgPosImagCoherence-avgPreImagCoherence;


%% kstest to find data distribution
findDataDistribution('Coherence',preCoherence,posCoherence,signifChannels)
fprintf("-----\n")
findDataDistribution('Imaginary Part of Coherence',preImagCoherence,posImagCoherence,signifChannels)


%%
avgPreImagCoherence(logical(eye(size(avgPreImagCoherence)))) = 0.08;

diffImagCoherence=avgPosImagCoherence-avgPreImagCoherence;

%% find connections that present significant differences between the neutral and surprise conditions
clc;
close all;
preCoherence=squeeze(preCoherence);
posCoherence=squeeze(posCoherence);


findSignificantConnections('Coherence',preCoherence,posCoherence,signifChannels,avgPreCoherence,avgPosCoherence,'yes')
findSignificantConnections('Imaginary Part of Coherence',preImagCoherence,posImagCoherence,signifChannels,avgPreImagCoherence,avgPosImagCoherence,'yes')


%% find channels pairs with higher difference values between the connectivity of the neutral and surprise conditions
clc;
close all;

%Pre (COH)
findHigherConnectivityConnections('Pre Stimulus','Coherence',avgPreCoherence,signifChannels,areasOfInterest);
%Pre (ICOH)
findHigherConnectivityConnections('Pre Stimulus','Imaginary Part of Coherence',avgPreImagCoherence,signifChannels,areasOfInterest)

%Pos (COH)
findHigherConnectivityConnections('Pos Stimulus','Coherence',avgPosCoherence,signifChannels,areasOfInterest)
%Pos (ICOH)
findHigherConnectivityConnections('Pos Stimulus','Imaginary Part of Coherence',avgPosImagCoherence,signifChannels,areasOfInterest)
%%
%Pos-Pre (COH)
findHigherConnectivityConnections('Pos-Pre','Coherence alpha band',diffCoherence,signifChannels,areasOfInterest)
%%

%Pos-Pre (ICOH)
findHigherConnectivityConnections('Pos-Pre','Imaginary Part of Coherence',diffImagCoherence,signifChannels,areasOfInterest)


%% Weighted Phase Lag Index
clc;
close all;
%Pre stimulus (neutral condition)
pre_wPLI = ones(numRelevChann, numRelevChann, size(neutralData,1), length(participants));

%Pos stimulus (surprise condition)
pos_wPLI = ones(numRelevChann, numRelevChann, length(participants));


for p=1:length(participants)
    for i = 1:numRelevChann
        for j = 1:numRelevChann
            if i~=j

                %Pre Stimulus
                for neutIdx=1:size(neutralData,1)

                    %cross-spectral power density
                    [pre_cpsd f] = cpsd(squeeze(neutralData_cell{p}(neutIdx,i,:)), squeeze(neutralData_cell{p}(neutIdx,j,:)),[],[],2048,fs);


                    %restrict to the frequency band of interest
                    freqIdx = find(f >= freqRange(1) & f <= freqRange(2));
                    pre_cpsd=pre_cpsd(freqIdx);

                    pre_wPLI(i, j, neutIdx, p) = abs((sum(abs(imag(pre_cpsd)).*sign(imag(pre_cpsd)))))/(sum(abs(imag(pre_cpsd))));
                end

                %Pos Stimulus
                pos_cpsd = cpsd(squeeze(stimulusData_cell{p}(:,i)), squeeze(stimulusData_cell{p}(:,j)),[],[],2048,fs);


                freqIdx = find(f >= freqRange(1) & f <= freqRange(2));
                pos_cpsd=pos_cpsd(freqIdx);

                pos_wPLI(i, j, p) = abs((sum(abs(imag(pos_cpsd)).*sign(imag(pos_cpsd)))))/(sum(abs(imag(pos_cpsd))));
            end
        end
    end
end

%average of all neutral instants
pre_wPLI=squeeze(mean(pre_wPLI,3));


%average of all the participants
avgPre_wPLI=mean(pre_wPLI,3);

avgPos_wPLI=mean(pos_wPLI,3);

%Pos-Pre (difference between the connectivity values of the neutral and surprise conditions)
diffwPLI=avgPos_wPLI-avgPre_wPLI;


%% kstest to find data distribution
findDataDistribution('Weighted Phase Lag Index',pre_wPLI,pos_wPLI,signifChannels)


%% find connections that present significant differences between the neutral and surprise conditions
clc;
close all;
findSignificantConnections('Weighted Phase Lag Index',pre_wPLI,pos_wPLI,signifChannels,avgPre_wPLI,avgPos_wPLI,'yes')


%% find channels pairs with higher difference values between the connectivity of the neutral and surprise conditions
clc;
%Pre (wPLI)
findHigherConnectivityConnections('Pre Stimulus','Weighted Phase Lag Index',avgPre_wPLI,signifChannels,areasOfInterest)

%Pos (wPLI)
findHigherConnectivityConnections('Pos Stimulus','Weighted Phase Lag Index',avgPos_wPLI,signifChannels,areasOfInterest)

%Pos-Pre (wPLI)
findHigherConnectivityConnections('Pos-Pre','Weighted Phase Lag Index gamma band',diffwPLI,signifChannels,areasOfInterest)




%% Mean Phase Coherence
close all;
clc;
%Pre stimulus (neutral condition)
pre_mpc = zeros(numRelevChann, numRelevChann, size(neutralData,1), length(participants));

%Pos stimulus (surprise condition)
pos_mpc = zeros(numRelevChann, numRelevChann, length(participants));


for p = 1:length(participants)
    for i = 1:numRelevChann
        for j = 1:numRelevChann

            %Pre Stimulus
            for neutIdx=1:size(neutralData,1)

                %Phase for the two channels
                pre_phase1=angle(hilbert(squeeze(neutralData_cell{p}(neutIdx,i,:))));

                pre_phase2=angle(hilbert(squeeze(neutralData_cell{p}(neutIdx,j,:))));


                %Phase differences between the two channels
                pre_phaseDiff = pre_phase1 - pre_phase2;

                pre_mpc(i,j,neutIdx,p) = abs(mean(exp(1i*pre_phaseDiff)));
            end

            %Pos Stimulus

            %Phase for the two channels
            pos_phase1=angle(hilbert(squeeze(stimulusData_cell{p}(:, i))));

            pos_phase2=angle(hilbert(squeeze(stimulusData_cell{p}(:, j))));


            %Phase differences between the two channels
            pos_phaseDiff = pos_phase1 - pos_phase2;

            pos_mpc(i,j,p) = abs(mean(exp(1i*pos_phaseDiff)));
        end
    end
end

%Average of all the neutral instants
pre_mpc=squeeze(mean(pre_mpc,3));

%Average the MCP values across participants
avgPre_mpc = mean(pre_mpc, 3);
avgPos_mpc = mean(pos_mpc, 3);

%Pos-Pre (difference between the connectivity values of the neutral and surprise conditions)
diffMPC=avgPos_mpc-avgPre_mpc;


%% kstest to find data distribution
findDataDistribution('Mean Phase Coherence',pre_mpc,pos_mpc,signifChannels)


%% find connections that present significant differences between the neutral and surprise conditions
close all;
clc;
findSignificantConnections('Mean Phase Coherence',pre_mpc,pos_mpc,signifChannels,avgPre_mpc,avgPos_mpc,'no')

%% find channels pairs with higher difference values between the connectivity of the neutral and surprise conditions
clc;
close all;
%Pre (MPC)
findHigherConnectivityConnections('Pre Stimulus','Mean Phase Coherence',avgPre_mpc,signifChannels,areasOfInterest)

%Pos (MPC)
findHigherConnectivityConnections('Pos Stimulus','Mean Phase Coherence',avgPos_mpc,signifChannels,areasOfInterest)

%Pos-Pre (MPC)
findHigherConnectivityConnections('Pos-Pre','Mean Phase Coherence',diffMPC,signifChannels,areasOfInterest)


%% Directed Transfer Function (eConnectome)
%NOTE:this connectivity method was implemented using the eConnectome toolbox
stim_data_cell=cell(length(participants),1);
for participant=1:length(participants)
    load('C:\Users\diogo\OneDrive - Universidade de Coimbra\Ambiente de Trabalho\dados_processados_2\p_'+string(participants(participant))+"_processed.mat",'fs','eeg_channels','each_processed_data') %load processed data
    stimData=(each_processed_data(stimulus_length,areasOfInterest))'; %temporal data corresponding to the surprise stimulus
    EEG.data=stimData;
    EEG.labels=signifChannels;
    EEG.type='EEG';
    EEG.nbchan=numRelevChann;
    EEG.points=length(stimData);
    EEG.srate=fs;
    EEG.labeltype='standard';
    EEG.unit='V';
    stim_data_cell{participant}=stimData;
    save('C:\Users\diogo\OneDrive - Universidade de Coimbra\Ambiente de Trabalho\DataStimulus\stim'+string(participant)+'.mat','EEG') %save in a mat variable compatible with the eConnectome toolbox
end



%% Average DTF matrix across participants
clc;
avgDTF=zeros(numRelevChann,numRelevChann,length(participants));

for part=1:length(participants)
    load('C:\Users\diogo\OneDrive - Universidade de Coimbra\Ambiente de Trabalho\DataDTF\dtf'+string(part)+'.mat','DTF') %load dtf file obtained from eConnectome
    avgDTF(:,:,part)=mean(DTF.matrix,3);
end

avgDTF=mean(avgDTF,3);

%Visualize connectivity matrix (heatmap)
figure;
imagesc(avgDTF);
colorbar;
xlabel('Channel'); ylabel('Channel');
xticklabels(signifChannels); yticklabels(signifChannels);
title('Directed Transfer Function (Pos Stimulus)');

