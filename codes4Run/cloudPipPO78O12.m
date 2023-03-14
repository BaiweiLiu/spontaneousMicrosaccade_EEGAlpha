%% Analysis scriptv for extract sacc events in delay phase in nn dataset
%% all piplinne 
%% It's always good to start wivth a clean sheet
clc, clear, close all, warning('off','all')

%% set project name
projectname = 'Extract sacc event left and right';

%% prepare Toolbox and data path 
platform = 'fgbp';

if strcmp(platform, 'hardDisk')
    plat_folder = '/Volumes/LaCie/EEG_data';
    parent_folder = [plat_folder filesep 'ReAnalysisFreek/EEG_eye_delay_NN/'];
    package_folder = '/Users/baiweil/Documents/MATLAB_pacakges/';
    
elseif strcmp(platform, 'laptop')
    plat_folder = '/Users/baiweil/AnalysisDesk';
    parent_folder = [plat_folder filesep 'ReAnalysisFreek/EEG_eye_delay_NN/'];
    package_folder = '/Users/baiweil/Documents/MATLAB_pacakges/';
    
elseif  strcmp(platform, 'fgbp')
    plat_folder = '/home/blu800/Desktop/analysisDesk';
    parent_folder = [plat_folder filesep 'ReAnalysisFreek/EEG_eye_delay_NN/'];
    package_folder = [plat_folder filesep 'MATLAB_pacakges/'];
end

cfg= [];
cfg.package_folder = package_folder;
cfg.futureTB_path = [package_folder 'futureToolBox'];
cd(cfg.futureTB_path);
futureTB_setup(cfg);

%% Set directions
write_dir_eeg = creatDir([parent_folder 'eeg_results']);
write_dir_fig = creatDir([write_dir_eeg filesep 'figures']);

write_dir_eye = creatDir([parent_folder 'eye_results']);
write_dir_be = creatDir([parent_folder 'behavior_results']);
write_dir_tf = creatDir([write_dir_eeg filesep 'data_tf']);
write_dir_MVPA = creatDir([write_dir_eeg filesep 'results_MVPA']);
write_dir_epoch = creatDir([write_dir_eeg filesep 'data_epoch']);
write_dir_cleanData = creatDir([write_dir_eeg filesep 'cleaned_epoch_data']);
write_dir_tf_reEpoch = creatDir([write_dir_eeg filesep 'data_tf_reEpoch']);

%%
trigs.left_en = [21 22];
trigs.right_en= [23 24];
trigs.left_hand = [21 23];
trigs.right_hand = [22 24];
trigs.epoch = [21:24];

%% channel and freqs selection
[channels,freqs] = eegFuture_selection();

%% 
threshold = 5;
minDis4Ms = 1;

%% epoch time 
epoch_time = [-3 1];

%% detect ms (toward vs. noMS vs. away)
data_file = [write_dir_eye filesep 'normalized_data/'];
sublist = get_subFiles(data_file);

for subInd = 1:length(sublist)

    cfg = [];
    cfg.read_file = sublist{subInd};
    cfg.smooth_raw = false; % important 
    cfg.smooth_der = true;
    cfg.smooth_der2 = false;
    cfg.smooth_step = 7;
    cfg.latency = [-3.0 1.0];
    cfg.ntrl = 5;
    cfg.threshold = threshold;
    cfg.winbef = [50 0];
    cfg.winaft = [50 100];
    cfg.minISI = 100;
    cfg.minDis = minDis4Ms;
    cfg.plotFig = false;
    cfg.output_dir = [write_dir_eye filesep 'saccDetect'];
    cfg.subjID = sublist{subInd}(end-7:end-4);
    eyeLab_mSaccDetect(cfg);
end

%%
cfg = [];
cfg.read_dir = [write_dir_eye filesep 'saccDetect'];
cfg.norm_dir = [write_dir_eye filesep 'normalized_data'];
cfg.time4fix = [-1.5 -0.5]; % -2000 ms - 0 ms 
cfg.output_dir = [write_dir_eye filesep 'saccDetect'];
eyeLab_mSaccAnalysis(cfg)

%% extract shift data

cfg=[];
cfg.time_i = [-1500, -500];
cfg.time = [-3000:1:1000];
cfg.input_dir = [write_dir_eye filesep 'saccDetect'];
cfg.output_dir = [write_dir_eye filesep 'group_results'];
cfg.outfile_name = ['GA_shift'];
cfg.condName = {'raw'};
eyeLab_GAshift(cfg)

%% re-epoch-data based on microsaccade
sublist = get_subFiles(write_dir_cleanData);

good_dir1 = [write_dir_eeg filesep 'data_clean_epoch/trial_keep'];
sublist_OKtl = get_subFiles(good_dir1);

sacc_dir = [write_dir_eye filesep 'saccDetect'];
sublist_shift = get_subFiles(sacc_dir);

gazeNorm_dir = [write_dir_eye filesep 'normalized_data'];
sublist_gaze = get_subFiles(gazeNorm_dir);

for subjInd = 1:length(sublist)
    
    infoDisp(['re epoch subj ' sublist{subjInd}(end-7:end-4)]);
    cfg = [];
    cfg.input_file = sublist{subjInd};
    cfg.goodness_file = {sublist_OKtl{subjInd}};    
    cfg.output_dir = [write_dir_eeg filesep 're_epoch'];
    cfg.fileName = [];
    cfg.time4sacc = [-1.5 -0.5];
    cfg.epochTime = [-1.2 1.2];
    cfg.trig_l = trigs.left_en;
    cfg.trig_r = trigs.right_en;
    v2struct(cfg); 
    
    sampleRate = 1000;
    eeg_time = [epochTime(1):0.001:epochTime(2)];
    epochTime = epochTime * sampleRate;
    
    out_file_corr = creatDir([output_dir filesep 'corr']);
    out_file_perm = creatDir([output_dir filesep 'perm']);
    
    % read and clean the data
    cfg = [];
    cfg.input_dir = input_file;
    cfg.goodness_file = goodness_file;
    [data,trialOK] = eegFuture_readER(cfg);
    
    rand_trial = data.trial(randperm(size(data.trial,1)),:,:);

    t4sacc = dsearchn(data.time',time4sacc')';
    
    % get shift info
    load(sublist_shift{subjInd});
    gazeL = gazeShift.gazeL(trialOK,:);
    gazeR = gazeShift.gazeR(trialOK,:);
    gazeStart = gazeShift.gaze_startPos(trialOK,:);
    gazeEnd = gazeShift.gaze_endPos(trialOK,:);
    thre_gazeVol_old = gazeShift.thre_gazeVol(trialOK,:);
    dwell_time_old = gazeShift.dwell_time(trialOK,:);
    dwell_time_return_old = gazeShift.dwell_time_return(trialOK,:);
    start_sacc_old = gazeShift.gaze_start(trialOK,:);
    return_sacc_old = gazeShift.gaze_return(trialOK,:);
    
    % get gaze postion 
    load(sublist_gaze{subjInd})
    gazePos = squeeze(eye_data.trial(trialOK,1,t4sacc(1):t4sacc(2)));
    
    numTrials = size(gazeL,1);
    
    % initial data
    %% get trials for left Sacc
    trialSize = 0;
    for tlInd = 1:numTrials
        infoDisp(['processing trial ' num2str(tlInd) '/' num2str(numTrials)],'loop',tlInd,  numTrials)
        shiftPos = [];
        % get shift time
        shiftTimeL = find(gazeL(tlInd,:));
        
        % check shift in timewindow
        shiftInTimeL = shiftTimeL((shiftTimeL>=t4sacc(1)) & (shiftTimeL<=t4sacc(2)));
        shiftPos =[gazeStart(tlInd,shiftInTimeL)',gazeEnd(tlInd,shiftInTimeL)'];
        
        for newtlInd = 1: length(shiftInTimeL)
            trialSize = trialSize +1;
        end
        
        %% get trials for right Sacc
        % get shift time
        shiftTimeR = find(gazeR(tlInd,:));
        
        % check shift in timewindow
        shiftInTimeR = shiftTimeR((shiftTimeR>=t4sacc(1)) & (shiftTimeR<=t4sacc(2)));
        shiftPos =[gazeStart(tlInd,shiftInTimeR)',gazeEnd(tlInd,shiftInTimeR)'];
        
        for newtlInd = 1: length(shiftInTimeR)
           trialSize = trialSize +1;
        end 
    end
    
    infoDisp(['will have '  num2str(trialSize)  ' trials'], 'line')
    
    %% get trials for left Sacc
    newTrial = zeros(trialSize,61,epochTime(2) - epochTime(1) +1); 
    thresholdInfo = zeros(1,trialSize);
    start_return = cell(1,trialSize);
    dwell_time = zeros(1,trialSize);
    dwell_time_return = zeros(1,trialSize);
    displace = zeros(1,trialSize);
    trialInfo =  zeros(1,trialSize);
    OldtrialNumber = zeros(1,trialSize);
    timetoEncoding=zeros(1,trialSize);
    delayTime =zeros(1,trialSize);
    medianGaze = zeros(1,trialSize);
    shiftStartEnd =zeros(trialSize,2);
    saccEvent=zeros(1,trialSize);
    toward = zeros(1,trialSize);
    fakeTrial=zeros(trialSize,61,epochTime(2) - epochTime(1) +1); 
    new_trial_ind = 0;
    
    for tlInd = 1:numTrials
        infoDisp(['processing trial ' num2str(tlInd) '/' num2str(numTrials)],'loop',tlInd,  numTrials)
        shiftPos = [];
        % get shift time
        shiftTimeL = find(gazeL(tlInd,:));
        
        % check shift in timewindow
        shiftInTimeL = shiftTimeL((shiftTimeL>=t4sacc(1)) & (shiftTimeL<=t4sacc(2)));
        shiftPos =[gazeStart(tlInd,shiftInTimeL)',gazeEnd(tlInd,shiftInTimeL)'];

        
        for newtlInd = 1: length(shiftInTimeL)
            new_trial_ind = new_trial_ind +1;
            % Epoch data
            starEpoch = shiftInTimeL(newtlInd) + epochTime(1); 
            endEpoch = shiftInTimeL(newtlInd) + epochTime(2);
            dwell_time_tl = dwell_time_old(tlInd,shiftInTimeL(newtlInd));
            dwell_time_return_tl = dwell_time_return_old(tlInd,shiftInTimeL(newtlInd));
            
            newTrial(new_trial_ind,:,:) = data.trial(tlInd,:,starEpoch:endEpoch);
            fakeTrial(new_trial_ind,:,:) = rand_trial(tlInd,:,starEpoch:endEpoch);
            
            saccEvent(new_trial_ind) = 1; % left
            thresholdInfo(new_trial_ind) = thre_gazeVol_old(tlInd,shiftInTimeL(newtlInd));
            if start_sacc_old(tlInd,shiftInTimeL(newtlInd))
                start_return{new_trial_ind} = 'start';
            elseif return_sacc_old(tlInd,shiftInTimeL(newtlInd))
                start_return{new_trial_ind} = 'return';
            end
            dwell_time(new_trial_ind) = dwell_time_tl;
            dwell_time_return(new_trial_ind) = dwell_time_return_tl;
            trialInfo(new_trial_ind) =  data.trialinfo(tlInd);
            
            OldtrialNumber(new_trial_ind) = tlInd;
            
            timetoEncoding(new_trial_ind) = shiftInTimeL(newtlInd) - 1000;
            
            medianGaze(new_trial_ind) = nanmedian(gazePos(tlInd,:));
            
            shiftStartEnd(new_trial_ind,:) = shiftPos(newtlInd,:);
            displace(new_trial_ind) = shiftPos(newtlInd,2) -  shiftPos(newtlInd,1);
            
            if  ismember(data.trialinfo(tlInd), trig_l) 
                toward(new_trial_ind) = 1;
            else
                toward(new_trial_ind) = 2;
            end
        end
        
        %% get trials for right Sacc
        % get shift time
        shiftTimeR = find(gazeR(tlInd,:));
        
        % check shift in timewindow
        shiftInTimeR = shiftTimeR((shiftTimeR>=t4sacc(1)) & (shiftTimeR<=t4sacc(2)));
        shiftPos =[gazeStart(tlInd,shiftInTimeR)',gazeEnd(tlInd,shiftInTimeR)'];
        
        for newtlInd = 1: length(shiftInTimeR)
           new_trial_ind = new_trial_ind +1;
            % Epoch data
            starEpoch = shiftInTimeR(newtlInd) + epochTime(1); 
            endEpoch = shiftInTimeR(newtlInd) + epochTime(2);
            dwell_time_tl = dwell_time_old(tlInd,shiftInTimeR(newtlInd));
            dwell_time_return_tl = dwell_time_return_old(tlInd,shiftInTimeR(newtlInd));
            
            if start_sacc_old(tlInd,shiftInTimeR(newtlInd))
                start_return{new_trial_ind} = 'start';
            elseif return_sacc_old(tlInd,shiftInTimeR(newtlInd))
                start_return{new_trial_ind} = 'return';
            end
            
            newTrial(new_trial_ind,:,:) = data.trial(tlInd,:,starEpoch:endEpoch);
            fakeTrial(new_trial_ind,:,:) = rand_trial(tlInd,:,starEpoch:endEpoch);
            thresholdInfo(new_trial_ind) = thre_gazeVol_old(tlInd,shiftInTimeR(newtlInd));
            
            dwell_time(new_trial_ind) = dwell_time_tl;
            dwell_time_return(new_trial_ind) = dwell_time_return_tl;
            
            saccEvent(new_trial_ind) = 2; 
            
            trialInfo(new_trial_ind) =  data.trialinfo(tlInd);
            
            OldtrialNumber(new_trial_ind) = tlInd;
            
            timetoEncoding(new_trial_ind) = shiftInTimeR(newtlInd) - 1000;
            
            medianGaze(new_trial_ind) = nanmedian(gazePos(tlInd,:));
            
            shiftStartEnd(new_trial_ind,:) = shiftPos(newtlInd,:);
            
            displace(new_trial_ind) = shiftPos(newtlInd,2) -  shiftPos(newtlInd,1);
            
            if  ismember(data.trialinfo(tlInd), trig_r) 
                toward(new_trial_ind) = 1;
            else
                toward(new_trial_ind) = 2;
            end
        end 
        

    end
    
    eeg_data = [];
    eeg_data.time = eeg_time;
    eeg_data.label = data.label;
    eeg_data.elec = data.elec;
    eeg_data.trial = newTrial;
    eeg_data.trialinfo = trialInfo';
    eeg_data.dimord = data.dimord;
    newTrial = [];
    
    cfg = [];
    cfg.resamplefs = 1000;
    eeg_data = ft_resampledata(cfg, eeg_data);
    
    eeg_data.saccEvent = saccEvent';
    eeg_data.OldtrialNumber = OldtrialNumber';
    eeg_data.timetoEncoding = timetoEncoding';
    eeg_data.medianGaze = medianGaze';
    eeg_data.shiftStartEnd = shiftStartEnd;
    eeg_data.toward= toward';
    eeg_data.dwell_time = dwell_time';
    eeg_data.dwell_time_return = dwell_time_return';
    eeg_data.start_return = start_return';
    eeg_data.thre_gazeVol = thresholdInfo';
    eeg_data.displace = displace';
    save([out_file_corr filesep input_file(end-7:end)], 'eeg_data');
    
    eeg_data = [];
    eeg_data.time = eeg_time;
    eeg_data.label = data.label;
    eeg_data.elec = data.elec;
    eeg_data.trial = fakeTrial;
    eeg_data.trialinfo = trialInfo';
    eeg_data.dimord = data.dimord;
    fakeTrial = [];
    
    cfg = [];
    cfg.resamplefs = 1000;
    eeg_data = ft_resampledata(cfg, eeg_data);
    
    eeg_data.saccEvent = saccEvent';
    eeg_data.OldtrialNumber = OldtrialNumber';
    eeg_data.timetoEncoding = timetoEncoding';
    eeg_data.medianGaze = medianGaze';
    eeg_data.shiftStartEnd = shiftStartEnd;
    eeg_data.toward= toward';
    eeg_data.dwell_time = dwell_time';
    eeg_data.dwell_time_return = dwell_time_return';
    eeg_data.start_return = start_return';
    eeg_data.thre_gazeVol = thresholdInfo';
    eeg_data.displace = displace';
    save([out_file_perm filesep input_file(end-7:end)], 'eeg_data');
    
end

%% do a time-frequency analysis
write_dir_tf_reEpoch_corr = creatDir([write_dir_eeg filesep 'data_tf_reEpoch_corr']);

output_dir = [write_dir_eeg filesep 're_epoch'];
out_file_corr = creatDir([output_dir filesep 'corr']);
sublist_corr = get_subFiles(out_file_corr);

% run through subj loop
for subjInd = 1: length(sublist_corr) 
    
    %% load cleaned data
    infoDisp(['load cleaned correct data for subject ' sublist_corr{subjInd}(end-7: end-4)]);
    load(sublist_corr{subjInd});


    %% Do TF analysis
    windowsize = 0.3; % 300 ms sliding time window

    cfg = [];
    cfg.taper = 'hanning'; % default hanning taper again
    cfg.method = 'mtmconvol'; % time-resolved frequency analysis (= short-time fourier transform)
    cfg.foi = [1:1:50]; % frequencies of interest (foi)
    cfg.toi = [eeg_data.time(1):0.02:eeg_data.time(end)]; % times of interest (i.e. where to centre the time-windows for each Fourier analysis on)
    cfg.t_ftimwin = ones(1, length(cfg.foi))*windowsize; % same window for each frequency.
    cfg.pad = 5; % padding the data to obtains the desired frequency resolution
    cfg.keeptrials = 'yes';
    cfg.output = 'fourier'; % get  to 'fourier' to get phase and phpower 
    data_tfr = ft_freqanalysis(cfg, eeg_data); % tfr.powspctrm = 4D (trials, electrodes, frequencies, timepoints => see also "tfr.dimord");
    
    data_tfr.saccEvent = eeg_data.saccEvent ;
    data_tfr.OldtrialNumber = eeg_data.OldtrialNumber ;
    data_tfr.timetoEncoding = eeg_data.timetoEncoding ;
    data_tfr.medianGaze = eeg_data.medianGaze;
    data_tfr.shiftStartEnd = eeg_data.shiftStartEnd;
    data_tfr.toward = eeg_data.toward;
    data_tfr.dwell_time = eeg_data.dwell_time;
    data_tfr.dwell_time_return = eeg_data.dwell_time_return;
    data_tfr.start_return = eeg_data.start_return;
    data_tfr.thre_gazeVol = eeg_data.thre_gazeVol;
    data_tfr.displace = eeg_data.displace;
    
    
    save([write_dir_tf_reEpoch_corr filesep sublist{subjInd}(end-7:end)], 'data_tfr');
    
    clear data_tfr
    clear eeg_data
end

%% Grand average alpha lateraztion data (PO7 and PO8)
sublist = get_subFiles(write_dir_tf_reEpoch_corr);
sublist_erp = get_subFiles([write_dir_eeg filesep 're_epoch' filesep 'corr']);

valueChannel = {['PO7'] ['PO8']};
outfile_dir = creatDir([write_dir_eeg filesep 'group_results']);
outfile_name = [outfile_dir filesep 'GA_saccSize_cvsi_TF-PO78' '.mat'];

% inti GA 
GA_struct = [];
GA_struct.power_cvsi = zeros(length(sublist), 3, 50, 121);
GA_struct.power_co = zeros(length(sublist), 3, 50, 121);
GA_struct.power_ip = zeros(length(sublist), 3, 50, 121);
GA_struct.power_co_log = zeros(length(sublist), 3, 50, 121);
GA_struct.power_ip_log = zeros(length(sublist), 3, 50, 121);

GA_struct.phaseLock_co = zeros(length(sublist), 3, 50, 121);
GA_struct.phaseLock_ip = zeros(length(sublist), 3, 50, 121);
GA_struct.phase_co = zeros(length(sublist), 3, 50, 121);
GA_struct.phase_ip = zeros(length(sublist), 3, 50, 121);

GA_struct.erp_co = zeros(length(sublist), 3, 2401);
GA_struct.erp_ip = zeros(length(sublist), 3, 2401);

% run through subj loop
for subjInd = 1 : length(sublist) 
    infoDisp(['processing subj ' sublist{subjInd}(end-7:end-4)])
    
     %%
    infoDisp(['step1: load dataset'], 'line')
    load(sublist{subjInd})

    %% trust trial
    tl_ok = ones(length(data_tfr.trialinfo),1);
    
    tl_ok = logical(tl_ok);
    %% select useful trials
    infoDisp(['step3: select the useful trials' '(' num2str(sum(tl_ok)) ')'], 'line')    
    
    data_tfr.fourierspctrm = data_tfr.fourierspctrm(tl_ok,:,:,:);
    
     %% 
    infoDisp(['step4: caculated power from fourierspectrm'], 'line')
    data_tfr.powspctrm = abs(data_tfr.fourierspctrm).^2;
    
    %% sels 
    
    channel_left = match_str(data_tfr.label, valueChannel{1});
    channel_right = match_str(data_tfr.label, valueChannel{2});
    
    left = ismember(data_tfr.saccEvent, 1);
    right = ismember(data_tfr.saccEvent, 2);
    
    left = left(tl_ok);
    right = right(tl_ok);
    
    startDis = abs(data_tfr.shiftStartEnd(tl_ok,1) - data_tfr.medianGaze(tl_ok));
    endDis = abs(data_tfr.shiftStartEnd(tl_ok,2) - data_tfr.medianGaze(tl_ok));
    away_sel  = startDis < endDis;
    return_sel = startDis > endDis;
    
    to_target_sel = data_tfr.toward ==1;
    away_target_sel = data_tfr.toward ==2;    
    to_target_sel = to_target_sel;
    away_target_sel = away_target_sel;   
    
    to_target_sel = to_target_sel(tl_ok);
    away_target_sel = away_target_sel(tl_ok);
    
    saccSize = abs(data_tfr.shiftStartEnd(tl_ok,1) - data_tfr.shiftStartEnd(tl_ok,2));
    cut_point = nanmedian(saccSize);

    smallSacc_sel = saccSize < cut_point;
    largeSacc_sel = saccSize > cut_point;
    
    % only use trial whose size < 1 visual degree 
    %all_tl = logical(ones(1,sum(tl_ok)))';
    all_tl = saccSize <= (100/5.7);
    
    %% tf lateration
    infoDisp(['step5: get group level cvsi'], 'line')
    
    %conds1 = {all_tl  smallSacc_sel largeSacc_sel};
    conds1 = {all_tl}
    conds2 = {all_tl away_sel  return_sel}; % raw, start MS, return MS
    %conds3 = {all_tl to_target_sel away_target_sel}; % raw, toward MS, away MS
    conds3 = {all_tl}

    for cond1_ind = 1:length(conds1)
        for cond2_ind = 1:length(conds2)
            for cond3_ind = 1:length(conds3)
                cond_sel =conds2{cond2_ind} & conds3{cond3_ind};
               
                a = mean(data_tfr.powspctrm(cond_sel & left,channel_right,:,:)); % contra-chR
                b = mean(data_tfr.powspctrm(cond_sel & right,channel_right,:,:)); % ipsi-chR
                c = mean(data_tfr.powspctrm(cond_sel & right,channel_left,:,:)); % contra-chL
                d = mean(data_tfr.powspctrm(cond_sel & left,channel_left,:,:)); % ipsi-chL

                cvsi_chR = squeeze(((a-b) ./ (a+b)) * 100);
                cvsi_chL = squeeze(((c-d) ./ (c+d)) * 100);
                GA_struct.power_cvsi(subjInd,cond2_ind,:,:) = (cvsi_chR + cvsi_chL) ./ 2;
                GA_struct.power_co(subjInd,cond2_ind,:,:) = squeeze((a + c) ./ 2);
                GA_struct.power_ip(subjInd,cond2_ind,:,:) = squeeze((b + d) ./ 2);
            end
        end
    end
    
    clear a b c d cvsi_chR cvsi_chL 
    %% tf lateration (log data)
    infoDisp(['step5: get group level loged c and i'], 'line')
    data_tfr.powspctrm = log(data_tfr.powspctrm);
    data_tfr.powspctrm =[]; % release memory space
    for cond1_ind = 1:length(conds1)
        for cond2_ind = 1:length(conds2)
            for cond3_ind = 1:length(conds3)
                cond_sel =conds2{cond2_ind} & conds3{cond3_ind};
               
                a = mean(data_tfr.powspctrm(cond_sel & left,channel_right,:,:)); % contra-chR
                b = mean(data_tfr.powspctrm(cond_sel & right,channel_right,:,:)); % ipsi-chR
                c = mean(data_tfr.powspctrm(cond_sel & right,channel_left,:,:)); % contra-chL
                d = mean(data_tfr.powspctrm(cond_sel & left,channel_left,:,:)); % ipsi-chL

                GA_struct.power_co_log(subjInd,cond2_ind,:,:) = squeeze((a + c) ./ 2);
                GA_struct.power_ip_log(subjInd,cond2_ind,:,:) = squeeze((b + d) ./ 2);
            end
        end
    end
    clear a b c d
    %% itpc extract
    infoDisp(['get phase coherernce from subj ' sublist{subjInd}(end-7:end-4)],'line')
    data_tfr.powspctrm = [];
    d = squeeze(data_tfr.fourierspctrm); % data of interest
    data_tfr.fourierspctrm = [];
    d_absnorm = d ./ abs(d); % each vector [in matrix] normalised to length 1
    clear d
    for cond1_ind = 1:length(conds1)
        for cond2_ind = 1:length(conds2)
            for cond3_ind = 1:length(conds3)
                cond_sel =conds2{cond2_ind} & conds3{cond3_ind};
                
                a = d_absnorm(cond_sel & left,channel_right,:,:); % contra-chR
                b = d_absnorm(cond_sel & right,channel_right,:,:); % ipsi-chR
                c = d_absnorm(cond_sel & right,channel_left,:,:); % contra-chL
                d = d_absnorm(cond_sel & left,channel_left,:,:); % ipsi-chL
                
                contra = squeeze(mean([a; c],1));
                ipsi = squeeze(mean([b; d],1));
                
                GA_struct.phaseLock_co(subjInd,cond2_ind,:,:) = abs(contra);
                GA_struct.phaseLock_ip(subjInd,cond2_ind,:,:) = abs(ipsi);
                
                GA_struct.phase_co(subjInd,cond2_ind,:,:) = contra;
                GA_struct.phase_ip(subjInd,cond2_ind,:,:) = ipsi;
                
            end
        end
    end
    
    if subjInd == length(sublist)
        GA_struct.time = data_tfr.time;
        GA_struct.freq = data_tfr.freq;
        GA_struct.eleLabel = data_tfr.label;
    end
    clear a b c d contra ipsi data_tfr
    %% ERP 
    infoDisp(['get ERP from subj ' sublist_erp{subjInd}(end-7:end-4)],'line')
    load(sublist_erp{subjInd})
    
    cfg = [];
    cfg.baseline     = [-0.5 -0.15];
    eeg_data_baseCorr = ft_timelockbaseline(cfg, eeg_data)
    
    clear eeg_data
    
    for cond1_ind = 1:length(conds1)
        for cond2_ind = 1:length(conds2)
            for cond3_ind = 1:length(conds3)
                cond_sel =conds2{cond2_ind} & conds3{cond3_ind};
                
                a = squeeze(mean(eeg_data_baseCorr.trial(cond_sel & left,channel_right,:))); % contra-chR
                b = squeeze(mean(eeg_data_baseCorr.trial(cond_sel & right,channel_right,:))); % ipsi-chR
                c = squeeze(mean(eeg_data_baseCorr.trial(cond_sel & right,channel_left,:))); % contra-chL
                d = squeeze(mean(eeg_data_baseCorr.trial(cond_sel & left,channel_left,:))); % ipsi-chL
                
                contra = (a+c)./2;
                ipsi = (b+d)./2;
                GA_struct.erp_co(subjInd,cond2_ind,:) = contra;
                GA_struct.erp_ip(subjInd,cond2_ind,:) = ipsi;
            end
        end
    end
    
    if subjInd == length(sublist)
        GA_struct.time_erp = eeg_data_baseCorr.time;
    end
    clear eeg_data_baseCorr
end

%
GA_struct.dimName = ['subjID_', 'cond1', 'cond2','cond3','chan','freq', 'time_'];
GA_struct.cond1Label = {'raw'};
GA_struct.cond2Label = {'raw' 'away' 'return'};
GA_struct.cond3Label = {'raw'};
save(outfile_name, 'GA_struct');

%% Grand average alpha lateraztion data (o1 vs. o2)

sublist = get_subFiles(write_dir_tf_reEpoch_corr);
sublist_erp = get_subFiles([write_dir_eeg filesep 're_epoch' filesep 'corr']);

valueChannel = {['O1'] ['O2']};
outfile_dir = creatDir([write_dir_eeg filesep 'group_results']);
outfile_name = [outfile_dir filesep 'GA_saccSize_cvsi_TF-O12' '.mat'];

% inti GA 
GA_struct = [];
GA_struct.power_cvsi = zeros(length(sublist), 3, 50, 121);
GA_struct.power_co = zeros(length(sublist), 3, 50, 121);
GA_struct.power_ip = zeros(length(sublist), 3, 50, 121);

GA_struct.power_co_log = zeros(length(sublist), 3, 50, 121);
GA_struct.power_ip_log = zeros(length(sublist), 3, 50, 121);

GA_struct.phaseLock_co = zeros(length(sublist), 3, 50, 121);
GA_struct.phaseLock_ip = zeros(length(sublist), 3, 50, 121);
GA_struct.phase_co = zeros(length(sublist), 3, 50, 121);
GA_struct.phase_ip = zeros(length(sublist), 3, 50, 121);

GA_struct.erp_co = zeros(length(sublist), 3, 2401);
GA_struct.erp_ip = zeros(length(sublist), 3, 2401);
% run through subj loop
for subjInd = 1 : length(sublist) 
    infoDisp(['processing subj ' sublist{subjInd}(end-7:end-4)])
    
     %%
    infoDisp(['step1: load dataset'], 'line')
    load(sublist{subjInd})

    %% trust trial
    tl_ok = ones(length(data_tfr.trialinfo),1);
    
    tl_ok = logical(tl_ok);
    %% select useful trials
    infoDisp(['step3: select the useful trials' '(' num2str(sum(tl_ok)) ')'], 'line')    
    
    data_tfr.fourierspctrm = data_tfr.fourierspctrm(tl_ok,:,:,:);
    
     %% 
    infoDisp(['step4: caculated power from fourierspectrm'], 'line')
    data_tfr.powspctrm = abs(data_tfr.fourierspctrm).^2;
    
    %% sels 
    
    channel_left = match_str(data_tfr.label, valueChannel{1});
    channel_right = match_str(data_tfr.label, valueChannel{2});
    
    left = ismember(data_tfr.saccEvent, 1);
    right = ismember(data_tfr.saccEvent, 2);
    
    left = left(tl_ok);
    right = right(tl_ok);
    
    startDis = abs(data_tfr.shiftStartEnd(tl_ok,1) - data_tfr.medianGaze(tl_ok));
    endDis = abs(data_tfr.shiftStartEnd(tl_ok,2) - data_tfr.medianGaze(tl_ok));
    away_sel  = startDis < endDis;
    return_sel = startDis > endDis;
    
    to_target_sel = data_tfr.toward ==1;
    away_target_sel = data_tfr.toward ==2;    
    to_target_sel = to_target_sel;
    away_target_sel = away_target_sel;   
    
    to_target_sel = to_target_sel(tl_ok);
    away_target_sel = away_target_sel(tl_ok);
    
    saccSize = abs(data_tfr.shiftStartEnd(tl_ok,1) - data_tfr.shiftStartEnd(tl_ok,2));
    cut_point = nanmedian(saccSize);

    smallSacc_sel = saccSize < cut_point;
    largeSacc_sel = saccSize > cut_point;
    
    % only use trial whose size < 1 visual degree 
    %all_tl = logical(ones(1,sum(tl_ok)))';
    all_tl = saccSize <= (100/5.7);
    
    %% tf lateration
    infoDisp(['step5: get group level cvsi'], 'line')
    
    %conds1 = {all_tl  smallSacc_sel largeSacc_sel};
    conds1 = {all_tl}
    conds2 = {all_tl away_sel  return_sel}; % raw, start MS, return MS
    %conds3 = {all_tl to_target_sel away_target_sel}; % raw, toward MS, away MS
    conds3 = {all_tl}

    for cond1_ind = 1:length(conds1)
        for cond2_ind = 1:length(conds2)
            for cond3_ind = 1:length(conds3)
                cond_sel =conds2{cond2_ind} & conds3{cond3_ind};
               
                a = mean(data_tfr.powspctrm(cond_sel & left,channel_right,:,:)); % contra-chR
                b = mean(data_tfr.powspctrm(cond_sel & right,channel_right,:,:)); % ipsi-chR
                c = mean(data_tfr.powspctrm(cond_sel & right,channel_left,:,:)); % contra-chL
                d = mean(data_tfr.powspctrm(cond_sel & left,channel_left,:,:)); % ipsi-chL

                cvsi_chR = squeeze(((a-b) ./ (a+b)) * 100);
                cvsi_chL = squeeze(((c-d) ./ (c+d)) * 100);
                GA_struct.power_cvsi(subjInd,cond2_ind,:,:) = (cvsi_chR + cvsi_chL) ./ 2;
                GA_struct.power_co(subjInd,cond2_ind,:,:) = squeeze((a + c) ./ 2);
                GA_struct.power_ip(subjInd,cond2_ind,:,:) = squeeze((b + d) ./ 2);
            end
        end
    end
    
    clear a b c d cvsi_chR cvsi_chL 
        %% tf lateration (log data)
    infoDisp(['step5: get group level loged c and i'], 'line')
    data_tfr.powspctrm = log(data_tfr.powspctrm);
    
    for cond1_ind = 1:length(conds1)
        for cond2_ind = 1:length(conds2)
            for cond3_ind = 1:length(conds3)
                cond_sel =conds2{cond2_ind} & conds3{cond3_ind};
               
                a = mean(data_tfr.powspctrm(cond_sel & left,channel_right,:,:)); % contra-chR
                b = mean(data_tfr.powspctrm(cond_sel & right,channel_right,:,:)); % ipsi-chR
                c = mean(data_tfr.powspctrm(cond_sel & right,channel_left,:,:)); % contra-chL
                d = mean(data_tfr.powspctrm(cond_sel & left,channel_left,:,:)); % ipsi-chL

                GA_struct.power_co_log(subjInd,cond2_ind,:,:) = squeeze((a + c) ./ 2);
                GA_struct.power_ip_log(subjInd,cond2_ind,:,:) = squeeze((b + d) ./ 2);
            end
        end
    end
    clear a b c d
    %% itpc extract
    infoDisp(['get phase coherernce from subj ' sublist{subjInd}(end-7:end-4)],'line')
    data_tfr.powspctrm = [];
    d = squeeze(data_tfr.fourierspctrm); % data of interest
    data_tfr.fourierspctrm = [];
    d_absnorm = d ./ abs(d); % each vector [in matrix] normalised to length 1
    clear d 
    for cond1_ind = 1:length(conds1)
        for cond2_ind = 1:length(conds2)
            for cond3_ind = 1:length(conds3)
                cond_sel =conds2{cond2_ind} & conds3{cond3_ind};
                
                a = d_absnorm(cond_sel & left,channel_right,:,:); % contra-chR
                b = d_absnorm(cond_sel & right,channel_right,:,:); % ipsi-chR
                c = d_absnorm(cond_sel & right,channel_left,:,:); % contra-chL
                d = d_absnorm(cond_sel & left,channel_left,:,:); % ipsi-chL
                
                contra = squeeze(mean([a; c],1));
                ipsi = squeeze(mean([b; d],1));
                
                GA_struct.phaseLock_co(subjInd,cond2_ind,:,:) = abs(contra);
                GA_struct.phaseLock_ip(subjInd,cond2_ind,:,:) = abs(ipsi);
                
                GA_struct.phase_co(subjInd,cond2_ind,:,:) = contra;
                GA_struct.phase_ip(subjInd,cond2_ind,:,:) = ipsi;
                
            end
        end
    end
    
    if subjInd == length(sublist)
        GA_struct.time = data_tfr.time;
        GA_struct.freq = data_tfr.freq;
        GA_struct.eleLabel = data_tfr.label;
    end
    clear a b c d contra ipsi data_tfr
     %% ERP 
    infoDisp(['get ERP from subj ' sublist_erp{subjInd}(end-7:end-4)],'line')
    load(sublist_erp{subjInd})
    
    cfg = [];
    cfg.baseline     = [-0.5 -0.15];
    eeg_data_baseCorr = ft_timelockbaseline(cfg, eeg_data)
    
    clear eeg_data
    
    for cond1_ind = 1:length(conds1)
        for cond2_ind = 1:length(conds2)
            for cond3_ind = 1:length(conds3)
                cond_sel =conds2{cond2_ind} & conds3{cond3_ind};
                
                a = squeeze(mean(eeg_data_baseCorr.trial(cond_sel & left,channel_right,:))); % contra-chR
                b = squeeze(mean(eeg_data_baseCorr.trial(cond_sel & right,channel_right,:))); % ipsi-chR
                c = squeeze(mean(eeg_data_baseCorr.trial(cond_sel & right,channel_left,:))); % contra-chL
                d = squeeze(mean(eeg_data_baseCorr.trial(cond_sel & left,channel_left,:))); % ipsi-chL
                
                contra = (a+c)./2;
                ipsi = (b+d)./2;
                GA_struct.erp_co(subjInd,cond2_ind,:) = contra;
                GA_struct.erp_ip(subjInd,cond2_ind,:) = ipsi;
            end
        end
    end
    
    if subjInd == length(sublist)
        GA_struct.time_erp = eeg_data_baseCorr.time;
    end
    clear eeg_data_baseCorr
end

%
GA_struct.dimName = ['subjID_', 'cond1', 'cond2','cond3','chan','freq', 'time_'];
GA_struct.cond1Label = {'raw'};
GA_struct.cond2Label = {'raw' 'away' 'return'};
GA_struct.cond3Label = {'raw'};
save(outfile_name, 'GA_struct');

%% Grand average TF data (topography)
sublist = get_subFiles(write_dir_tf_reEpoch_corr);
valueChannel = {channels.visualL channels.visualR};
outfile_dir = creatDir([write_dir_eeg filesep 'group_results']);
outfile_name = [outfile_dir filesep 'GA_saccSize_TF' '.mat'];

% inti GA 
GA_struct = [];
GA_struct.data_itpc = zeros(length(sublist), 3, 3, 3, 61, 50 , 120);

% run through subj loop
for subjInd = 1 : length(sublist) 
    infoDisp(['processing subj ' sublist{subjInd}(end-7:end-4)])
    
     %%
    infoDisp(['step1: load dataset'], 'line')
    load(sublist{subjInd})
    
    %% short 2.25 s long 4.25 s
    infoDisp(['step2: Get trust trials based on time window 3000'], 'line')
    
    %% trust trial
    tl_ok = ones(length(data_tfr.trialinfo),1);
    
    tl_ok = logical(tl_ok);
    
    %% select useful trials
    infoDisp(['step3: select the useful trials' '(' num2str(sum(tl_ok)) ')'], 'line')    
    
    data_tfr.fourierspctrm = data_tfr.fourierspctrm(tl_ok,:,:,:);
    
    %% 
    infoDisp(['step4: caculated power from fourierspectrm'], 'line')
    data_tfr.powspctrm = abs(data_tfr.fourierspctrm).^2;
    
    %% presacc baseline correct 
    infoDisp(['step5: Baseline correct'], 'line')
    cfg =[];
    cfg.data_tfr = data_tfr;
    cfg.windowsize = 0.3;
    cfg.baseWin = [-1 -0.5];

    data_tfr = eegFuture_baseLineTF(cfg)
    
    %% sels 
    
    left = ismember(data_tfr.saccEvent, 1);
    right = ismember(data_tfr.saccEvent, 2);
    
    left = left(tl_ok);
    right = right(tl_ok);
    
    startDis = abs(data_tfr.shiftStartEnd(tl_ok,1) - data_tfr.medianGaze(tl_ok));
    endDis = abs(data_tfr.shiftStartEnd(tl_ok,2) - data_tfr.medianGaze(tl_ok));
    away_sel  = startDis < endDis;
    return_sel = startDis > endDis;
    
    to_target_sel = data_tfr.toward ==1;
    away_target_sel = data_tfr.toward ==2;    
    to_target_sel = to_target_sel;
    away_target_sel = away_target_sel;   
    
    to_target_sel = to_target_sel(tl_ok);
    away_target_sel = away_target_sel(tl_ok);
    
    saccSize = abs(data_tfr.shiftStartEnd(tl_ok,1) - data_tfr.shiftStartEnd(tl_ok,2));
    cut_point = nanmedian(saccSize);
    smallSacc_sel = saccSize < cut_point;
    largeSacc_sel = saccSize > cut_point;
    
    all_tl = logical(ones(1,sum(tl_ok)))';
    
    %% power
    infoDisp(['step4: get group level power'], 'line')
    
    conds1 = {all_tl smallSacc_sel largeSacc_sel};
    conds2 = {all_tl away_sel  return_sel}; % raw, start MS, return MS
    conds3 = {all_tl to_target_sel  away_target_sel}; % raw, toward MS, away MS
    
    for cond1_ind = 1:length(conds1)
        for cond2_ind = 1:length(conds2)
            for cond3_ind = 1:length(conds3)
    
            cond_sel = conds1{cond1_ind} & conds2{cond2_ind} & conds3{cond3_ind};
        
            GA_struct.raw_pow(subjInd,cond1_ind,cond2_ind,cond3_ind,:,:,:) = squeeze(mean(data_tfr.powspctrm(cond_sel,:,:,:)));
            GA_struct.left_pow(subjInd,cond1_ind,cond2_ind,cond3_ind,:,:,:) = squeeze(mean(data_tfr.powspctrm(cond_sel&left,:,:,:)));
            GA_struct.right_pow(subjInd,cond1_ind,cond2_ind,cond3_ind,:,:,:) = squeeze(mean(data_tfr.powspctrm(cond_sel&right,:,:,:)));
            GA_struct.diff_pow(subjInd,cond1_ind,cond2_ind,cond3_ind,:,:,:) = GA_struct.left_pow(subjInd,cond1_ind,cond2_ind,cond3_ind,:,:,:) - ...
                                                                    GA_struct.right_pow(subjInd,cond1_ind,cond2_ind,cond3_ind,:,:,:);
            end
        end
    end
    
end

%
GA_struct.time = data_tfr.time;
GA_struct.freq = data_tfr.freq;
GA_struct.eleLabel = data_tfr.label;
GA_struct.dimName = ['subjID_', 'cond1', 'cond2','cond3','chan','freq', 'time_'];
GA_struct.cond2Label = {'raw' 'smallSacc' 'largeSacc'};
GA_struct.cond2Label = {'raw' 'away' 'return'};
GA_struct.cond3Label = {'raw' 'toward' 'away'};
save(outfile_name, 'GA_struct');

%% Grand average ITPC topography 

sublist = get_subFiles(write_dir_tf_reEpoch_corr);
sublist_erp = get_subFiles([write_dir_eeg filesep 're_epoch' filesep 'corr']);

outfile_name = [write_dir_GA filesep 'GA_saccSize_itpc_topo' '.mat'];

% inti GA 
GA_struct = [];

GA_struct.phaseLock_left = zeros(length(sublist), 3, 61, 50, 121);
GA_struct.phaseLock_right = zeros(length(sublist), 3, 61, 50, 121);
GA_struct.phaseLock_all = zeros(length(sublist), 3, 61, 50, 121);

% run through subj loop
for subjInd = 1 : length(sublist) 
    infoDisp(['processing subj ' sublist{subjInd}(end-7:end-4)])
    
     %%
    infoDisp(['step1: load dataset'], 'line')
    load(sublist{subjInd})

    %% trust trial
    tl_ok = ones(length(data_tfr.trialinfo),1);
    
    tl_ok = logical(tl_ok);
    %% select useful trials
    infoDisp(['step3: select the useful trials' '(' num2str(sum(tl_ok)) ')'], 'line')    
    
    data_tfr.fourierspctrm = data_tfr.fourierspctrm(tl_ok,:,:,:);
    
    %% sels
    
    left = ismember(data_tfr.saccEvent, 1);
    right = ismember(data_tfr.saccEvent, 2);
    
    left_sel = left(tl_ok);
    right_sel = right(tl_ok);
    
    startDis = abs(data_tfr.shiftStartEnd(tl_ok,1) - data_tfr.medianGaze(tl_ok));
    endDis = abs(data_tfr.shiftStartEnd(tl_ok,2) - data_tfr.medianGaze(tl_ok));
    away_sel  = startDis < endDis;
    return_sel = startDis > endDis;
    
    saccSize = abs(data_tfr.shiftStartEnd(tl_ok,1) - data_tfr.shiftStartEnd(tl_ok,2));
    
    % only use trial whose size < 1 visual degree 
    %all_tl = logical(ones(1,sum(tl_ok)))';
    all_tl = saccSize <= (100/5.7);
    
    %% get ITPC topo 
    infoDisp(['step5: get group level itpc for electrodes'], 'line')

    data_coh = squeeze(data_tfr.fourierspctrm); % data of interest
    data_tfr.fourierspctrm = [];
    d_absnorm = data_coh ./ abs(data_coh); % each vector [in matrix] normalised to length 1
    clear data_coh
    
    conds1 = {all_tl away_sel&all_tl  return_sel&all_tl}; % raw, start MS, return MS
    
    for cond1_ind = 1:length(conds1)
            cond_sel =conds1{cond1_ind} & all_tl;
            
            a = squeeze(mean(d_absnorm(cond_sel & left,:,:,:),1)); % contra-chR
            b = squeeze(mean(d_absnorm(cond_sel & right,:,:,:),1)); % ipsi-chR
            c = squeeze(mean(d_absnorm(cond_sel,:,:,:),1)); % ipsi-chR
            
            GA_struct.phaseLock_left(subjInd,cond1_ind,:,:,:) = abs(a);
            GA_struct.phaseLock_right(subjInd,cond1_ind,:,:,:) = abs(b);
            GA_struct.phaseLock_all(subjInd,cond1_ind,:,:,:) = abs(c);

    end
    
    if subjInd == length(sublist)
        GA_struct.time = data_tfr.time;
        GA_struct.freq = data_tfr.freq;
        GA_struct.eleLabel = data_tfr.label;
        GA_struct.elec = data_tfr.elec;
    end
    clear a b c data_coh d_absnorm data_tfr
end

%
GA_struct.dimName = ['subjID_', 'cond1','chan','freq', 'time_'];

GA_struct.cond2Label = {'raw' 'away' 'return'};

save(outfile_name, 'GA_struct');
