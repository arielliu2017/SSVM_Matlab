clc;
clear;
close all;

load SERIES_DATA_WITHOUT_COORDS.mat

%% EXAMPLE # 1
% HAVE A LOOK AT README
% DEE2 H3N2

INDICES.study = 2;
INDICES.subject = ':';
INDICES.shedding = 1;
INDICES.symp = 1;
INDICES.early_treatment = 3; % NA
INDICES.sham = 1;
INDICES.gender = 1;
INDICES.age = ':';
INDICES.log_symp_score = ':';
INDICES.time = ':';

sample_indx = get_sample_indx_wrapper(INDICES, ANNOTATION_DATA, ALL_IDS);

pw_it = 5;
pw_mat_1 = cell2mat(SERIES_DATA_WITHOUT_COORDS(pw_it, sample_indx))';
assert(size(PW_DATA{pw_it},1) == size(pw_mat_1,2));
assert(size(sample_indx,1) == size(pw_mat_1,1));


%% EXAMPLE # 2

INDICES.study = 2;
INDICES.subject = 122;
INDICES.shedding = 1;
INDICES.symp = 1;
INDICES.early_treatment = 3; % NA
INDICES.sham = 1;
INDICES.gender = 1;
INDICES.age = ':';
INDICES.log_symp_score = ':';
INDICES.time = 13;
sample_indx = get_sample_indx_wrapper(INDICES, ANNOTATION_DATA, ALL_IDS);
assert(sample_indx == 498);

PW_IT = 5;
pw_mat_2 = cell2mat(SERIES_DATA_WITHOUT_COORDS(PW_IT, sample_indx))';