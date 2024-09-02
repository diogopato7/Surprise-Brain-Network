%% -------------------- Preprocessing --------------------

% Clear workspace and command window
clear all;
clc;
close all;

% Prompt the user for participant numbers
participantes_str = input('Select the numbers of the participants separated by commas (or "all" for all the 30 participants): ', 's');




if strcmp(participantes_str, "all")
    participants = 1:30;
else
    participants = str2num(strrep(participantes_str, ',', ' ')); % Convert string to numbers
end

% Initialize cell array to store file paths
selected_files = cell(1, numel(participants));

% Define folder path and list of files
folder_path = "C:\Users\diogo\OneDrive - Universidade de Coimbra\Ambiente de Trabalho\dados_eeg";
file_list = dir(fullfile(folder_path, '*.TRC'));

% Extract file paths for selected participants
for i = 1:numel(participants)
    selected_files{i} = fullfile(folder_path, file_list(participants(i)).name);
end


patients_data = cell(length(selected_files), 1); % cell initialization to store all data
ecg_data = cell(length(selected_files), 1); % cell initialization to store ECG data
filteredData=cell(length(patients_data),1); %cell initialization to store the filtered data
processedData_cell=cell(1, length(selected_files)); %cell initialization to store the processed data
Zfica_cell=cell(length(selected_files),1);
T_cell=cell(length(selected_files),1);

%video_start=[44,32]; %Time of the video initialization
video_start=[14,27,51,20,30,10,41,18,32,29,54,10,11,26,11,13,21,33,10,9,10,9,12,16,12,12,14,13,22,10]; %Time of the video initialization
video_start_selecionados = video_start(participants);
video_duration=235;%video duration in sec

for i=1:length(selected_files) %for all the files selected
    %Open and convert EEG file from TRC to mat file
    file_name = selected_files{i};
    Data = trc_file(file_name);
    Data.get_electrode_info();
    patients_data{i}=Data.def_data_access(Data.a_n_data_secs, 5, Data.a_file_elec_cell)';
    channels = Data.a_file_elec_cell; %channel order
    eeg_channels = channels(1:19);
    fs = Data.a_samp_freq; %sampling frequency

    patients_data{i}=patients_data{i}(video_start_selecionados(i)*fs:(video_start_selecionados(i)*fs+(video_duration*fs))-1,:);
    ecg_data{i}=patients_data{i}(:,20:21);
    patients_data{i}=patients_data{i}(:,1:19);


    % for j=1:length(eeg_channels)
    %     figure((i-1)*19+j)
    %     subplot(2,1,1);
    %     plot(0:1/fs:1/fs*(length(patients_data{i}(:,j))-1), patients_data{i}(:,j));
    %     xlabel('Time (s)'); ylabel('Voltage');
    %     title(string(eeg_channels(j))+'- Time Domain')
    % 
    %     subplot(2,1,2)
    %     f = [-numel(patients_data{i}(:,j))/2:numel(patients_data{i}(:,j))/2-1].*fs/numel(patients_data{i}(:,j));
    %     plot(f, abs(fftshift(fft(patients_data{i}(:,j)))), 'm');
    %     xlabel('Frequency (Hz)'); ylabel('Fourier Transform');
    %     title(string(eeg_channels(j))+'- Frequency Domain')
    % end

end


%% Filtering
clc;
close all;
%iirnotc filter
fc=[37 50]; %cut frequencies (remove 37Hz and 50Hz components)
wo=fc/(fs/2);
bw=wo/35; %bandwidth
[b, a]=iirnotch(wo(1), bw(1)); %notch filter- 37Hz
[c, d]=iirnotch(wo(2), bw(2)); %notch filter- 50Hz
fc2=[0.5 80];
for i=1:length(selected_files)
    filteredData{i}=filtfilt(b, a, patients_data{i}); %apply filter to the EEG data
    filteredData{i}=filtfilt(c, d, filteredData{i}); %apply filter to the EEG data

    %bandpass filter
    filteredData{i}=bandpass(filteredData{i}, fc2, fs, ImpulseResponse="fir");



    % for j=1:length(eeg_channels)
    %     figure((i-1)*19+j)
    %     subplot(2,1,1);
    %     plot(0:1/fs:1/fs*(length(filteredData{i}(:,j))-1), filteredData{i}(:,j));
    %     xlabel('Time (s)'); ylabel('Voltage');
    %     title(string(eeg_channels(j))+'- Time Domain')
    %     subplot(2,1,2)
    %     f = [-numel(filteredData{i}(:,j))/2:numel(filteredData{i}(:,j))/2-1].*fs/numel(filteredData{i}(:,j));
    %     plot(f, abs(fftshift(fft(filteredData{i}(:,j)))), 'm');
    %     xlabel('Frequency (Hz)'); ylabel('Fourier Transform');
    %     title(string(eeg_channels(j))+'- Frequency Domain')
    %
    % end

end

%% ECG Filtering
clc;
close all;

ecg_2_filtered_data=cell(1,length(selected_files));
locs_tempo_all = cell(1, length(selected_files));
intervals_all = cell(1, length(selected_files));
bpm_all = cell(1, length(selected_files));


for l=1:length(selected_files)
    each_ecg_data=ecg_data{l}(:,2);
    fs=512;
    % t=1:90;
    % plot(t,each_ecg_data(length(t)));


    timeWindow=2;
    N=fs*timeWindow;
    b=(1/N)*ones(1, N);
    ecg_1 = filtfilt (b, 1, each_ecg_data);
    ecg_1=each_ecg_data-ecg_1;

    t = 0:1/fs:30;
    fn=0.5*fs;
    order = 8;
    fc_ecg=30;
    wc_ecg = fc_ecg/fn;
    [b1_ecg,a1_ecg]=butter(8,wc_ecg);
    ecg_2=filtfilt(b1_ecg,a1_ecg,ecg_1);
    %plot(t,ecg_2(1:length(t)));
    if size(filteredData{l},2)~=20
        filteredData{l} = [filteredData{l}, ecg_2];
        ecg_2_filtered_data{l}=ecg_2;
    end
    % R Peak Detection
    wc=25;
    fc_2=wc/(0.5*fs);
    [b1_ecg,a1_ecg]=butter(order,fc_2);
    ecg_e2=filtfilt(b1_ecg,a1_ecg,ecg_2);
    wc=5;
    fc_2=wc/(0.5*fs);
    [b1_ecg_2,a1_ecg_2]=butter(order,fc_2,'High');
    ecg_2=filter(b1_ecg_2,a1_ecg_2,ecg_e2);


    ecg_3=diff(ecg_2);
    ecg_4=ecg_3.^2;
    timeWindow =0.02;
    N=round(fs*timeWindow);
    b2=(1/N)*ones(1, N);
    a2 = 1;
    ecg_5 = filter (b2, a2, ecg_4);
    ecg_threshold = 0.7 * mean (ecg_5);


    figure();
    plot(t,ecg_5(1:length(t)));
    yline(ecg_threshold)

    RR_intervals = cell(1, numel(ecg_data));

    [pk,locs] = findpeaks(ecg_5,'MinPeakHeight',ecg_threshold,'MinPeakDistance',250);
    intervals = cell(1,numel(locs)-1);
    for j=1:numel(locs)-1
        intervals{j} = locs(j):locs(j+1);
    end
    locs_tempo=[];
    for k=1:length(locs)
        locs_tempo(k)=locs(k)/512;
    end
    locs_tempo=locs_tempo';
    locs_tempo_all{l} = locs_tempo;
    intervals_all{l} = intervals;


    % figure()
    % plot(ecg_5(intervals{2}))
    % yline(ecg_threshold)
    % legend('ecg 5', 'threshold ecg')

    stimulus_time=60;
    [indice_antes_estimulo,~]=find(locs_tempo<stimulus_time,1,'last');
    indice_antes_estimulo=indice_antes_estimulo-2;

    [initial_index,~]=find(locs_tempo>45,1,'first');

    [indice_final,~]=find(locs_tempo<75,1,'last');
    valor_final=locs_tempo(indice_final);

    time_beats_before=locs_tempo(indice_antes_estimulo)-locs_tempo(initial_index);
    bpm_before=(indice_antes_estimulo-initial_index+1)*60/time_beats_before;

    time_beats_after=locs_tempo(indice_final)-locs_tempo(indice_antes_estimulo+2);
    bpm_after=(indice_final-indice_antes_estimulo-1)*60/time_beats_after;
    bpm_all{l} = [bpm_before, bpm_after];

    fprintf('Participant %d: Before the stimulus: %d After the stimulus: %d\n================================\n', participants(l), round(bpm_before), round(bpm_after));
    % Plotar o sinal de ECG
    plot(t, ecg_5(1:length(t)));
    hold on;

    locs_in_range = locs_tempo(locs_tempo <= t(end));

    plot(locs_in_range, ecg_5(locs(1:length(locs_in_range))), 'ro', 'MarkerSize', 10);

    title('ECG Peak Detection');
    xlabel('Time (s)');
    ylabel('ECG');

    yline(ecg_threshold)
    legend('ECG signal', 'Detected Peaks','ECG Threshold');

    hold off
end

%% Features espectrais do ECG

clc

vlf_power_all = cell(1, length(selected_files));
lf_power_all = cell(1, length(selected_files));
hf_power_all = cell(1, length(selected_files));
lf_hf_ratio_all = cell(1, length(selected_files));
total_power_all = cell(1, length(selected_files));

for l = 1:length(selected_files)
    stimulus_time = 60;
    post_stimulus_duration = 15;

    % Dados de ECG
    ecg_data = ecg_2_filtered_data{l}; % ECG filtrado para um participante
    t = (0:length(ecg_data)-1)/fs; % Vetor de tempo

    % Índices antes e depois do estímulo
    idx_before = (t > 45) & (t < 60);
    idx_after = (t >= stimulus_time) & (t <= stimulus_time + post_stimulus_duration);

    % Sinais antes e depois do estímulo
    ecg_before = ecg_data(idx_before);
    ecg_after = ecg_data(idx_after);

    % Cálculo da PSD usando o método de Welch
    [psd_before, f_before] = pwelch(ecg_before, [], [], [], fs);
    [psd_after, f_after] = pwelch(ecg_after, [], [], [], fs);

    % Frequency bands
    vlf_band = [0.0033 0.04];
    lf_band = [0.04 0.15];
    hf_band = [0.15 0.4];

    vlf_power_before = bandpower(psd_before, f_before, vlf_band,'psd');
    lf_power_before = bandpower(psd_before, f_before, lf_band,'psd');
    hf_power_before = bandpower(psd_before, f_before, hf_band,'psd');

    vlf_power_after = bandpower(psd_after, f_after, vlf_band,'psd');
    lf_power_after = bandpower(psd_after, f_after, lf_band,'psd');
    hf_power_after = bandpower(psd_after, f_after, hf_band,'psd');

    lf_hf_ratio_before = lf_power_before / hf_power_before;
    lf_hf_ratio_after = lf_power_after / hf_power_after;

    total_power_before = vlf_power_before + lf_power_before + hf_power_before;
    total_power_after = vlf_power_after + lf_power_after + hf_power_after;

    vlf_power_all{l} = [vlf_power_before, vlf_power_after];
    lf_power_all{l} = [lf_power_before, lf_power_after];
    hf_power_all{l} = [hf_power_before, hf_power_after];
    lf_hf_ratio_all{l} = [lf_hf_ratio_before, lf_hf_ratio_after];
    total_power_all{l} = [total_power_before, total_power_after];
   


    fprintf('================================\n');
    fprintf('Participant %d before the stimulus:\n', participants(l));
    fprintf('VLF Power: %.4f\n', vlf_power_before);
    fprintf('LF Power: %.4f\n', lf_power_before);
    fprintf('HF Power: %.4f\n', hf_power_before);
    fprintf('LF/HF Ratio: %.4f\n', lf_hf_ratio_before);
    fprintf('Total Power: %.4f\n', total_power_before);

    fprintf('================================\n');
    fprintf('Participant %d after the stimulus:\n', participants(l));
    fprintf('VLF Power: %.4f\n', vlf_power_after);
    fprintf('LF Power: %.4f\n', lf_power_after);
    fprintf('HF Power: %.4f\n', hf_power_after);
    fprintf('LF/HF Ratio: %.4f\n', lf_hf_ratio_after);
    fprintf('Total Power: %.4f\n', total_power_after);
end

%% HRV time domain
clc

RMSSD_all = cell(1, length(selected_files));
SDNN_all = cell(1, length(selected_files));
pNN50_all = cell(1, length(selected_files));

for l=1:length(selected_files)
    
    RMSSD_before = zeros(1, numel(RR_intervals));
    SDNN_before = zeros(1, numel(RR_intervals));
    pNN50_before = zeros(1, numel(RR_intervals));
    RMSSD_after = zeros(1, numel(RR_intervals));
    SDNN_after = zeros(1, numel(RR_intervals));
    pNN50_after = zeros(1, numel(RR_intervals));

    for i = 1:numel(intervals_all{l})
        interval_after_stimulus = locs_tempo_all{l} > stimulus_time & locs_tempo_all{l} <= (stimulus_time + 15);
        interval_before_stimulus = locs_tempo_all{l} <= (stimulus_time - 15) & locs_tempo_all{l} < stimulus_time ;

        RR_intervals_before = intervals_all{l}{i}(interval_before_stimulus);
        RMSSD_before(i) = sqrt(mean(diff(RR_intervals_before).^2));
        SDNN_before(i) = std(RR_intervals_before);
        pNN50_before(i) = sum(abs(diff(RR_intervals_before)) > 0.05) / length(RR_intervals_before) * 100;

        RR_intervals_after = intervals_all{l}{i}(interval_after_stimulus);
        RMSSD_after(i) = sqrt(mean(diff(RR_intervals_after).^2));
        SDNN_after(i) = std(RR_intervals_after);
        pNN50_after(i) = sum(abs(diff(RR_intervals_after)) > 0.05) / length(RR_intervals_after) * 100;
    end

    RMSSD_all{l} = [mean(RMSSD_before), mean(RMSSD_after)];
    SDNN_all{l} = [mean(SDNN_before), mean(SDNN_after)];
    pNN50_all{l} = [mean(pNN50_before), mean(pNN50_after)];

    

    fprintf('Participante %d antes do estímulo:\n', participants(l));
    fprintf('RMSSD: %.2f\n', RMSSD_before(l));
    fprintf('SDNN: %.2f\n', SDNN_before(l));
    fprintf('pNN50: %.2f\n', pNN50_before(l));
    fprintf('Participante %d depois do estímulo:\n', participants(l));
    fprintf('RMSSD: %.2f\n', RMSSD_after(l));
    fprintf('SDNN: %.2f\n', SDNN_after(l));
    fprintf('pNN50: %.2f\n', pNN50_after(l));
    fprintf('=======================\n');
end
    
    
%% ECG Tachogram

stimulus_time = 60; 
stimulus_time_2=180;
stimulus_time_3=218;

for i =1:length(participants)
RR_intervals = diff(locs_tempo_all{i});

RR_times = locs_tempo_all{i}(2:end);
%RR_times_mid = (locs_tempo(1:end-1) + locs_tempo(2:end)) / 2;
smooth_RR_intervals = smooth(RR_intervals, 0.1, 'rloess'); % Ajuste o fator de suavização conforme necessário



figure;
plot(RR_times, RR_intervals, '-o', 'LineWidth', 1);
hold on;


xline(stimulus_time, 'r--', 'Firts Stimulus', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5);
xline(stimulus_time_2, 'r--', 'Second Stimulus', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5);
xline(stimulus_time_3, 'r--', 'Third Stimulus', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5);



xlabel('Time (s)');
ylabel('RR Intervals (s)');
title('ECG Tachogram');
grid on;

%ylim([0, 1.2]); 

hold off;


legend('RR Intervals', 'Stimuli Times');

end

%% ICA
clc;
close all;

r = 20; % Number of independent components to be calculated

if length(participants)==1
    selected_participants=participants;
else
    fprintf(' Choose between these participants: ');
    fprintf('%d ', participants);
    fprintf('\n');
    participant_indexes = input('Select the participants numbers to do the ICA (separated by commas): ', 's');
    selected_participants = str2num(strrep(participant_indexes, ',', ' '));
end
lista = [];
for i=1:length(participants)
    for j=1:length(selected_participants)
        if participants(i) == selected_participants(j);
            lista = [lista i];
        end
    end
end

Zfica_cell = cell(size(selected_participants));
T_cell = cell(size(selected_participants));
W_cell=cell(size(selected_participants));



for i = 1:length(lista)
    [Zfica, W, T, mu] = fastICA(filteredData{lista(i)}', r,'negentropy');
    Zfica_cell{i} = Zfica;
    T_cell{i} = T;
    W_cell{i}=W;
end
%% ICA Components for each selected participant
close all;
fc2 = [0.5 80];
clc;


if length(selected_participants)==1
    selected_participants_ICA=selected_participants;
else
    fprintf('Selected Participants: ');
    disp(selected_participants)
    abc = sprintf(' Select the participant to see the ICA graphics  (separated by commas): ');
    user_choice = input(abc, 's');
    selected_participants_ICA = str2num(strrep(user_choice, ',', ' '));
end

lista2 = [];
for i=1:length(selected_participants)
    for j=1:length(selected_participants_ICA)
        if selected_participants(i) == selected_participants_ICA(j);
            lista2 = [lista2 i];
        end
    end
end

for i = 1:length(lista2)
    Zfica = Zfica_cell{lista2(i)};
    for j = 1:r
        figure; subplot(2,1,1);
        plot(0:1/fs:1/fs*(length(Zfica(j,:))-1), Zfica(j,:),'-'); xlim([0 200]); xlabel('Time (s)'); ylabel('Voltage'); title(['Component ' num2str(j) ' - Time Domain']);
        subplot(2,1,2);
        f = [-numel(Zfica(j,:))/2:numel(Zfica(j,:))/2-1].*fs/numel(Zfica(j,:));
        plot(f, abs(fftshift(fft(Zfica(j,:)))), 'm'); xlim(fc2); xlabel('Frequency (Hz)'); ylabel('Fourier Transform'); title(['Component ' num2str(j) ' - Frequency Domain']);
    end
end


correlations = zeros(20, 1);

for i = 1:20
    correlacao = corrcoef(ecg_2_filtered_data{lista}, Zfica(i,:));
    correlations(i) = correlacao(1,2);
end
correlations=abs(correlations);
[sorted_correlations, sorted_indexes] = sort(correlations, 'descend');
top_5_correlations = sorted_correlations(1:5);
top_5_indexes = sorted_indexes(1:5);

% Imprime as três maiores correlações e as componentes correspondentes
for i = 1:5
    fprintf("Correlation between the component %d and the ECG signal: %.4f\n", top_5_indexes(i), top_5_correlations(i));
end


%% Reconstruction of the EEG data without the noisy IC
clc;
close all;
C={};

ruidosos = cell(1, length(selected_participants_ICA));
fprintf('Participantes selecionados: ');
disp(selected_participants_ICA);

if length(selected_participants_ICA)==1
    selected_participants_noisy=selected_participants_ICA;
else

    prompt = sprintf(' Dos selected_participants quais voce quer escolher a componente ruidosa? (separados por vírgula): ');
    user_input = input(prompt, 's');



    selected_participants_noisy = str2num(strrep(user_input, ',', ' '));

end

for i=1:length(selected_participants_noisy)

    prompt2 = sprintf('Quais componentes ruidosos você quer eliminar do participante %d ? (separados por vírgula): ', selected_participants_noisy(i));
    user_input2 = input(prompt2, 's');
    componentes_ruidosos = str2num(user_input2);
    ruidosos{i} = componentes_ruidosos;
    C = ruidosos;
    

    C={[1,10],[6,8,11]}; % c= noisy components identified in the graphs
    T_cell{lista2(i)}(:,C{i})=0; %eliminate the noisy components

    processedData=T_cell{i}/W_cell{i}'*Zfica_cell{i}; %reconstruct the data without noise
    processedData_cell{i}=processedData;


    % % Plot of reconstructed data after ICA
    % for j=1:length(eeg_channels)
    %     %time
    %     figure; subplot(2,1,1);
    %     plot(0:1/fs:1/fs*(length(processedData_cell{i}(:,j))-1), processedData_cell{i}(:,j)); xlabel('Time (s)'); ylabel('Voltage'); title(string(eeg_channels(j))+'- Time Domain')
    %
    %     %frequency
    %     subplot(2,1,2);
    %     f = [-numel(processedData_cell{i}(:,j))/2:numel(processedData_cell{i}(:,j))/2-1].*fs/numel(processedData_cell{i}(:,j));
    %     plot(f, abs(fftshift(fft(processedData_cell{i}(:,j)))), 'm'); xlim(fc); xlabel('Frequency (Hz)'); ylabel('Fourier Transform'); title(string(eeg_channels(j))+'- Frequency Domain')
    % end
end
%%
clc;
close all;
%Save the preprocessed EEG data
% file_name = join(string(participants), "_") + "_processed.mat";
% save("C:\Users\diogo\OneDrive - Universidade de Coimbra\Ambiente de Trabalho\"+file_name,'fs','eeg_channels','processedData_cell')
for i = 1: length(selected_participants_noisy)
    each_processed_data=processedData_cell{i}';
    save("C:\Users\diogo\OneDrive - Universidade de Coimbra\Ambiente de Trabalho\dados_processados_2\p_"+selected_participants_noisy(i)+"_processed.mat",'fs','eeg_channels','each_processed_data')
end
%% Plot spectrogram using STFT
clc;
close all;
for participant = 14:14
load('C:\Users\diogo\OneDrive - Universidade de Coimbra\Ambiente de Trabalho\dados_processados_2\p_'+string(participant)+"_processed.mat",'fs','eeg_channels','each_processed_data') %load processed data
    processedData_cell{participant}=each_processed_data(:,1:19);
fc2=[0.5 50];
    for i=1:length(eeg_channels) %for all the eeg_channels
        figure;
        spectrogram(processedData_cell{participant}(:,i),1/(1/fs),0.5/(1/fs),5/(1/fs),fs,'yaxis');
        ylim(fc2);
        title("Short Time Fourier Transform - Channel "+eeg_channels(i))
    end
end
%% Plot spectrogram using CWT
clc;
clear all;
load("C:\Users\diogo\OneDrive - Universidade de Coimbra\Ambiente de Trabalho\dados_processados\p_7_processed.mat",'fs','eeg_channels','each_processed_data')
for i=1:19
    figure(i)
    cwt(each_processed_data(:,i),fs);
end
% figure;
% imagesc(0:1/fs:video_duration, 1:numel(oi), abs(oi));
% colorbar;