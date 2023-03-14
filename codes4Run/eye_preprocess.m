%% Analysis script for preprocess eye gaze data 

%% It's always good to start with a clean sheet
clc, clear, close all, warning('off','all')

%% set project name
projectname = 'Re analysis eye data in delay phase';

%% prepare Toolbox
parent_folder = '/Users/baiweil/AnalysisDesk/ReAnalysisFreekData/EEG_eye_delay/';
package_folder = '/Users/baiweil/Documents/MATLAB_pacakges/';

cfg= [];
cfg.package_folder = package_folder;
cfg.futureTB_path = [package_folder 'futureToolBox'];
cd(cfg.futureTB_path);
futureTB_setup(cfg);

%% Set directions
read_dir1 = [parent_folder  'data/eye_raw/exp1_sess1_eye'];
read_dir2 = [parent_folder  'data/eye_raw/exp1_sess2_eye'];

write_dir_eye = creatDir([parent_folder 'eye_results']);
write_dir_fig = creatDir([write_dir_eye filesep 'figures']);

%%
trigs.left_en = [21 22];
trigs.right_en= [23 24];
trigs.left_hand = [21 23];
trigs.right_hand = [22 24];
trigs.epoch = [21:24];

%% set parameter
epoch_time = [-3 1.0]; % the first is prestim, the second is poststim
localiser_trig = [201:209];
localiser_epoch_time = [-0.5 1.5];

%% step 1: epoch trial and localizer data without interp
sublist1 = get_subFiles(read_dir1,'*.asc');
sublist2 = get_subFiles(read_dir2,'*.asc');

for subjInd = 1:length(sublist1)
    cfg =[];
    cfg.input_dirs = {[sublist1{subjInd}] [sublist2{subjInd}]};
    cfg.output_dir = write_dir_eye;
    cfg.epoch_trig = trigs.epoch;
    cfg.epoch_time = epoch_time;
    cfg.localiser_trig = localiser_trig;
    cfg.localiser_epoch_time = localiser_epoch_time;
    cfg.subjID = sublist1{subjInd}(end-8:end-7);
    cfg.maxBinkDuration = 0; % default is 500, 0 means no interp
    eyeLab_preprocess(cfg)
end

%% step 2: mark the NaN cluster in interesting time window
% readDir = [write_dir_eye filesep 'epoched_data']; 
% sublist = get_subFiles(readDir);
% 
% for subjInd = 1:length(sublist)
%     cfg =[];
%     cfg.input_dir = sublist{subjInd};
%     cfg.output_dir = write_dir_eye;
%     cfg.time_i = [0.2 0.6]; % 200 - 600 ms after probe
%     eyeLab_markNaN(cfg)
% end

%% step 3: localiser
readDir = [write_dir_eye filesep 'localiser_data']; 
wr_dir_locl = creatDir([write_dir_eye filesep 'normalized_cal']); 
sublist = get_subFiles(readDir);

for subjInd = 1:length(sublist)
    
    cfg = [];
    cfg.input_dir = sublist{subjInd};
    cfg.localiser_triggers =[201 204 207 205 203 206 209];
    cfg.plotTitle = sublist{subjInd}(end-7:end-4);
    cfg.averg_time = [0.5 1];
    cfg.sc_size = [1680,1050];
    cfg.output_dir = wr_dir_locl;
    eyeLab_localiser(cfg);
end

%% step4: normalize the data
rd_dir_epoch = [write_dir_eye filesep 'epoched_data']; 
rd_dir_locl = [write_dir_eye filesep 'normalized_cal/data']; 

sublist_epoch=get_subFiles(rd_dir_epoch);
sublist_locl=get_subFiles(rd_dir_locl);

for subjInd = 1:length(sublist_epoch)
    
    cfg=[];
    cfg.epoch_dir = sublist_epoch{subjInd};
    cfg.locl_dir = sublist_locl{subjInd};
    cfg.output_dir = [write_dir_eye filesep 'normalized_data'];
    cfg.plotTitle = sublist_epoch{subjInd}(end-7:end-4);
    data = eyeLab_normalize(cfg);
    
    cfg=[];
    cfg.eye_data = data;
    cfg.noisyLine = [50 50]; % one value for one channel
    cfg.channel = {'eyeX'};
    cfg.time_i = [0 0.6]; % the window that used to detect ms 
    cfg.plotTitle =  sublist_epoch{subjInd}(end-7:end-4);
    cfg.output_dir = [write_dir_eye filesep 'trialClean_autoDetect'];
    eyeLab_detectNoisy(cfg);
end


