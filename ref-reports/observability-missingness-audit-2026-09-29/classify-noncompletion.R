# Uses existing reconstruction; emits aggregates and hashed IDs only.
sink(tempfile('endline-reconstruction-'))
source('ref-reports/observability-missingness-audit-2026-09-29/reconcile-endline-completions.R')
sink()
nc <- reached %>% filter(!person %in% final$person)
people <- nc %>% group_by(person) %>% summarise(
 records=n(), recorded_nonconsent=any(consent==0,na.rm=TRUE),
 recorded_unable=any(interview==0,na.rm=TRUE),
 .groups='drop') %>% mutate(category=case_when(
 recorded_nonconsent~'Recorded non-consent',
 recorded_unable~'Unable to interview (no recorded non-consent)',
 TRUE~'Incomplete status fields'))
stopifnot(nrow(nc)==103L,nrow(people)==98L)
counts<-people%>%count(category)
write_csv(counts,file.path(out,'noncompletion-person-counts.csv'))
write_csv(people%>%mutate(person_hash=map_chr(person,~digest::digest(.x,algo='sha256',serialize=FALSE)))%>%select(-person),file.path(out,'noncompletion-person-audit.csv'))
unknown <- nc %>% filter(person%in%people$person[people$category=='Incomplete status fields'])
write_csv(unknown%>%count(interview,consent),file.path(out,'incomplete-person-status-counts.csv'))
cat('Survey records:',nrow(nc),'unique noncompleters:',nrow(people),'\n')
print(counts)
cat('Overlap of non-consent and unable flags:\n');print(people%>%count(recorded_nonconsent,recorded_unable))
cat('Incomplete status combinations:\n');print(unknown%>%count(interview,consent))
