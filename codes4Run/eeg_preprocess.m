%% Analysis script for preprocess ERP  

%% It's always good to start with a clean sheet
clc, clear, close all, warning('off','all')

%% set project name
projectname = 'preprocess the EEG data';

%% prepare Toolbox
parent_folder = '/Users/baiweil/AnalysisDesk/ReAnalysisFreekData/EEG_eye_delay/';
package_folder = '/Users/baiweil/Documents/MATLAB_pacakges/';

cfg= [];
cfg.package_folder = package_folder;
cfg.futureTB_path = [package_folder 'futureToolBox'];
cd(cfg.futureTB_path);
futureTB_setup(cfg);

%% Set directions
write_dir_eeg = creatDir([parent_folder 'eeg_results']);
write_dir_fig = creatDir([write_dir_eeg filesep 'figures']);

read_dir1 = [parent_folder  'data' filesep 'EEG_raw/session1'];
read_dir2 = [parent_folder  'data' filesep 'EEG_raw/session2'];

write_dir_eye = creatDir([parent_folder 'eye_results']);

write_dir_tf = creatDir([write_dir_eeg filesep 'data_tf']);
write_dir_MVPA = creatDir([write_dir_eeg filesep 'results_MVPA']);
write_dir_epoch = creatDir([write_dir_eeg filesep 'data_epoch']);
write_dir_cleanData = creatDir([write_dir_eeg filesep 'cleaned_epoch_data']);

%%
trigs.left_en = [21 22];
trigs.right_en= [23 24];
trigs.left_hand = [21 23];
trigs.right_hand = [22 24];
trigs.epoch = [21:24];

%% channel and freqs selection
[channels,freqs] = eegFuture_selection();

%% epoch time 
epoch_time = [-3 1];

%% epoch 
%% clean data and merge data from 2 sessions
subList1 = get_subFiles(read_dir1,'*.cnt');
subList2 = get_subFiles(read_dir2,'*.cnt');

% run through subj loop
for subjInd = 1:length(subList1)
    
    %% subject specific information
    disp([' ']);  disp([' ']);      disp(['getting data from ', subList1{subjInd}(71:73), ' and ', subList2{subjInd}(71:74) ]);    disp([' ']);  disp([' ']);
    
    %% step 1 read in data in all at once and re-reference (for .cnt data works better to first read in all data, and later epoch it)
    cfg = [];
    cfg.dataset = subList1{subjInd};
    cfg.dataformat = 'ns_cnt32';
    cfg.headerformat = 'ns_cnt32';
    cfg.eventformat = 'ns_cnt32';
    cfg.reref = 'yes';
    cfg.implicitref = 'LM'; % left mastoid during recording
    cfg.refchannel = {'LM','RM'};
    cfg.dftfilter = 'yes';
    data1 = ft_preprocessing(cfg);
    
    % load session 2
    cfg.dataset = subList2{subjInd};
    data2 = ft_preprocessing(cfg);
    
    
    %% step 2 segment into trials of interest
    cfg = [];
    cfg.dataset = subList1{subjInd};
    cfg.trialdef.eventtype = 'trigger';
    cfg.trialdef.eventvalue = trigs.epoch;
    cfg.trialdef.prestim = -epoch_time(1);
    cfg.trialdef.poststim = epoch_time(2);
    cfg = ft_definetrial(cfg);
    cfg.trl(:,[1,2]) = cfg.trl(:,[1,2])+24; % to realign to photo-diode measurement: maximal deflection at 24 ms after trigger
    data1 = ft_redefinetrial(cfg, data1);
    
    
    cfg = [];
    cfg.dataset = subList2{subjInd};
    cfg.trialdef.eventtype = 'trigger';
    cfg.trialdef.eventvalue = trigs.epoch;
    cfg.trialdef.prestim = -epoch_time(1);
    cfg.trialdef.poststim = epoch_time(2);
    cfg = ft_definetrial(cfg);
    cfg.trl(:,[1,2]) = cfg.trl(:,[1,2])+24; % to realign to photo-diode measurement: maximal deflection at 24 ms after trigger
    data2 = ft_redefinetrial(cfg, data2);
        
    
    %% fix some noisy channels before appending the sessions
    
    % fix noisy channels in ds1 only
    if ismember(subjInd, 1)
        badchan = ismember(data1.label, {'FC5'});
        replacingchans = ismember(data1.label, {'FC3','FC1'});
        for trl = 1:length(data1.trial)
            data1.trial{trl}(badchan,:) = squeeze(mean(data1.trial{trl}(replacingchans,:)));
        end
    elseif ismember(subjInd, 2)
        badchan = ismember(data1.label, {'C1'});
        replacingchans = ismember(data1.label, {'C3','Cz'});
        for trl = 1:length(data2.trial)
            data1.trial{trl}(badchan,:) = squeeze(mean(data1.trial{trl}(replacingchans,:)));
        end
    end
    
    % fix noisy channels in ds2 only
    if ismember(subjInd, 18)
        badchan = ismember(data2.label, {'F8'});
        replacingchans = ismember(data2.label, {'AF3','F6'});
        for trl = 1:length(data2.trial)
            data2.trial{trl}(badchan,:) = squeeze(mean(data2.trial{trl}(replacingchans,:)));
        end
    end
    
    % fix noisy chan in both ds1 and ds2
    if ismember(subjInd, 25)
        badchan = ismember(data1.label, {'AF8'});
        replacingchans = ismember(data1.label, {'F8','Fp2'});
        for trl = 1:length(data1.trial)
            data1.trial{trl}(badchan,:) = squeeze(mean(data1.trial{trl}(replacingchans,:)));
        end
        for trl = 1:length(data2.trial)
            data2.trial{trl}(badchan,:) = squeeze(mean(data2.trial{trl}(replacingchans,:)));
        end
    end
    
    %% append datasets now  
    data2save = ft_appenddata(cfg, data1, data2);
    
    %% remove irrelevant channels
    cfg = [];
    cfg.channel = {'all','-LM','-RM'};
    data2save = ft_selectdata(cfg, data2save);
    
    %% correct EEG labels so fieltrip understand them properly
    lab2check = {'FP1','FP2','FPz','Cpz'};
    labcorrect = {'Fp1','Fp2','Fpz','CPz'};
    errchan = ismember(data2save.label, lab2check);
    whicherr = ismember(lab2check,data2save.label(errchan));
    data2save.label(errchan) = labcorrect(whicherr);
    data = data2save; clear data2save
    
    %% save appended data
    save([write_dir_epoch filesep 'pp' subList1{subjInd}(86:87) '.mat'], 'data');
end

%% clean data through ICA and visual remove 
sublist = get_subFiles(write_dir_epoch);


% run through subj loop
for subjInd = 1:length(sublist)
    
    cfg = [];
    
    % input and output file
    cfg.epoch_dir = sublist{subjInd};
    cfg.output_dir = write_dir_eeg;
    cfg.subj_name = sublist{subjInd}(end-7:end-4);
    
    % clean setting
    cfg.ica_here = true; % whether ICA has already been down before
    cfg.plot_out = false; % whether plot results
    cfg.ICA_onlyEEG = false; % whether run ICA including the EOG
    cfg.autoICA = true; % whether remove blink compoents auto 
    cfg.autoICA_co_remove = 0.4; % the correlation value used to remove blink compoents auto 
    cfg.i_channel = [channels.visualL channels.visualR channels.motionL channels.motionR];
    cfg.i_freq = [8 30];
    cfg.do_visualRemove = true; % whether clean the data by visual remove
    
    eegFuture_preprocess(cfg);  
end

%% get the cleaned data
cd([write_dir_eeg filesep 'data_clean_epoch/Epoch_data'])
sublist=dir('*.mat');
sublist={sublist.name}; 

% run through subj loop
for subjInd = 1 :length(sublist) 
    %% read data
    cfg.data_set = write_dir_eeg;
    cfg.subj_name = sublist{subjInd}(1:4);
    cfg.apply_ica = true;
    cfg.apply_VR = false;    
    cfg.apply_laplacian = true;
    cfg.ICA_onlyEEG = false;
    eeg_data = eegFuture_readData(cfg);
    
    save([write_dir_cleanData filesep sublist{subjInd}], 'eeg_data')
end