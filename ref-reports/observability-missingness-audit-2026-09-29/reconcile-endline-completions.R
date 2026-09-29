# Read-only reconstruction of the saved endline sample, using production filters.
.libPaths(c(.libPaths(), '/home/ed/R/x86_64-pc-linux-gnu-library/4.5'))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(lubridate))
out <- 'ref-reports/observability-missingness-audit-2026-09-29'
raw <- read_csv('data/raw-data/Endline Survey.csv', guess_max=10000, show_col_types=FALSE) %>%
  mutate(SubmissionDate=parse_datetime(SubmissionDate, '%b %d, %Y %I:%M:%S %p',
    locale=locale(tz='America/New_York')) %>% format(tz='Africa/Nairobi') %>%
    parse_datetime(locale=locale(tz='Africa/Nairobi')))
saved <- readRDS('data/clean-data/clean-endline-data.rds')
stages <- list(raw=raw)
stages$present_nonweb_after_start <- raw %>% filter(deviceid!='(web)', present==1, SubmissionDate>='2016-10-18')
stages$gps_available <- stages$present_nonweb_after_start %>% filter(!is.na(`gps-Longitude`),!is.na(`gps-Latitude`))
stages$valid_fieldwork <- stages$gps_available %>% filter(cluster_id!=1163 | SubmissionDate>='2016-11-14',enumerator!=111)
stages$eligible_contact_records <- stages$valid_fieldwork %>% filter(date(SubmissionDate)!='2016-11-7' | test==1)
reached <- stages$eligible_contact_records
stages$completed_records <- reached %>% filter(if_all(c(present,interview,consent),~!is.na(.x)&.x==1))
stages$unique_completed_respondents <- stages$completed_records %>% arrange(person,SubmissionDate) %>% distinct(person,.keep_all=TRUE)
final <- stages$unique_completed_respondents
stopifnot(nrow(final)==3678L, setequal(final$person,saved$KEY.individ),setequal(final$KEY,saved$KEY))
# Verify each person's selected survey key, not merely membership of each set.
key_check <- final %>% select(person,KEY) %>% inner_join(saved %>% select(KEY.individ,KEY),by=c('person'='KEY.individ'),suffix=c('.rebuilt','.saved'))
stopifnot(nrow(key_check)==3678L,all(key_check$KEY.rebuilt==key_check$KEY.saved))
counts <- imap_dfr(stages,function(x,stage) bind_rows(
 tibble(stage=stage,sms='all',records=nrow(x),unique_people=n_distinct(x$person)),
 x %>% group_by(sms) %>% summarise(records=n(),unique_people=n_distinct(person),.groups='drop') %>% mutate(stage=stage,.before=1)))
write_csv(counts,file.path(out,'completion-flow.csv'))
status <- reached %>% count(interview,consent,name='records')
write_csv(status,file.path(out,'contact-status-counts.csv'))
# A record-level audit using hashes instead of respondent IDs.
record_audit <- raw %>% transmute(record_hash=map_chr(KEY,~digest::digest(.x,algo='sha256',serialize=FALSE)),
 person_hash=map_chr(person,~digest::digest(.x,algo='sha256',serialize=FALSE)),
 sms, outcome=case_when(
 !KEY %in% stages$present_nonweb_after_start$KEY ~ 'not_present_nonweb_after_start',
 !KEY %in% stages$gps_available$KEY ~ 'missing_gps',
 !KEY %in% stages$valid_fieldwork$KEY ~ 'fieldwork_validity_exclusion',
 !KEY %in% reached$KEY ~ 'november_7_test_filter',
 !KEY %in% stages$completed_records$KEY ~ 'interview_or_consent_not_affirmative',
 !KEY %in% final$KEY ~ 'repeat_completed_submission',
 TRUE ~ 'retained'))
write_csv(record_audit,file.path(out,'completion-record-audit.csv'))
# Deduplicated contact accounting: someone with any eligible completed interview counts as completed.
people <- reached %>% group_by(person) %>% summarise(has_completion=any(KEY%in%final$KEY),.groups='drop')
cat('Verified exact saved person/submission-key pairs:',nrow(final),'\n')
cat('Eligible contact records:',nrow(reached),'unique people:',nrow(people),'\n')
cat('Records without affirmative interview and consent:',nrow(reached)-nrow(stages$completed_records),'\n')
cat('Completed records:',nrow(stages$completed_records),'repeat completed records:',nrow(stages$completed_records)-nrow(final),'\n')
cat('Unique contacted people without an eligible completion:',sum(!people$has_completion),'\n')
print(counts %>% filter(stage%in%c('eligible_contact_records','completed_records','unique_completed_respondents')))
print(status)
