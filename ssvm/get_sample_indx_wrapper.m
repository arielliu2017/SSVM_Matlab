function [common_sample_indx] = get_sample_indx_wrapper(INDICES, ANNOTATION_DATA, ALL_IDS)

total_indx_set = (1:size(ANNOTATION_DATA.SAMPLEID, 1))';

% CAN BE CONVERTED TO A LOOP, BUT IS THIS EASY FOR DEBUGGING???

common_sample_indx = total_indx_set;

study_sample_indx = get_sample_indx(cellstr(ANNOTATION_DATA.STUDYID), ALL_IDS.study_ids, INDICES.study);
common_sample_indx = intersect(common_sample_indx, study_sample_indx);

subject_sample_indx = get_sample_indx(cellstr(ANNOTATION_DATA.SUBJECTID), ALL_IDS.subject_ids, INDICES.subject);
common_sample_indx = intersect(common_sample_indx, subject_sample_indx);

shedding_sample_indx = get_sample_indx(ANNOTATION_DATA.SHEDDING_SC1, ALL_IDS.shedding_ids, INDICES.shedding, 'INT');
common_sample_indx = intersect(common_sample_indx, shedding_sample_indx);

symp_sample_indx = get_sample_indx(ANNOTATION_DATA.SYMPTOMATIC_SC2, ALL_IDS.symp_ids, INDICES.symp, 'INT');
common_sample_indx = intersect(common_sample_indx, symp_sample_indx);

early_treatment_sample_indx = get_sample_indx(cellstr(ANNOTATION_DATA.EARLYTX), ALL_IDS.early_treatment_ids, INDICES.early_treatment);
common_sample_indx = intersect(common_sample_indx, early_treatment_sample_indx);

sham_sample_indx = get_sample_indx(cellstr(ANNOTATION_DATA.SHAM), ALL_IDS.sham_ids, INDICES.sham);
common_sample_indx = intersect(common_sample_indx, sham_sample_indx);

gender_sample_indx = get_sample_indx(cellstr(ANNOTATION_DATA.GENDER), ALL_IDS.gender_ids, INDICES.gender);
common_sample_indx = intersect(common_sample_indx, gender_sample_indx);

age_sample_indx = get_sample_indx(ANNOTATION_DATA.AGE, ALL_IDS.age_ids, INDICES.age, 'INT');
common_sample_indx = intersect(common_sample_indx, age_sample_indx);

log_symp_score_sample_indx = get_sample_indx(ANNOTATION_DATA.LOGSYMPTSCORE_SC3, ALL_IDS.log_symp_score_ids, INDICES.log_symp_score, 'INT');
common_sample_indx = intersect(common_sample_indx, log_symp_score_sample_indx);

time_sample_indx = get_sample_indx(ANNOTATION_DATA.TIMEHOURS, ALL_IDS.time_ids, INDICES.time, 'INT');
common_sample_indx = intersect(common_sample_indx, time_sample_indx);

end