.libPaths(c(.libPaths(), '/home/ed/R/x86_64-pc-linux-gnu-library/4.5'))
d <- readRDS('data/clean-data/clean-endline-data.rds')
k <- readRDS('data/clean-data/clean-endline-know-table-data.rds')
cat('Clean endline respondents:', nrow(d), '\n')
print(table(SMS=d$sms.treatment, module_record=d$in_know_table))
print(table(module_type=d$survey.type, module_record=d$in_know_table,useNA='always'))
no <- d[d$sms.treatment=='sms.control',]
j <- match(no$KEY.individ,k$KEY.individ)
cat('No-SMS module types:\n');print(table(k$know.table.type[j],useNA='always'))
a <- !is.na(j) & k$know.table.type[j]=='table.A'
cat('No-SMS A with record:',sum(a),'with positive recognition:',sum(a & k$obs_know_person[j]>0,na.rm=TRUE),'\n')
cat('All knowledge summary records:\n');print(table(k$know.table.type))
cat('Records matched to clean endline:\n');print(table(k$know.table.type,k$KEY.individ%in%d$KEY.individ))
raw <- read.csv('data/raw-data/Endline Survey.csv',check.names=FALSE)
ra <- read.csv('data/raw-data/Endline Survey-survey-sec_D-tableA.csv',check.names=FALSE)
rb <- read.csv('data/raw-data/Endline Survey-survey-sec_D-tableB.csv',check.names=FALSE)
r <- raw[raw$person %in% no$KEY.individ & raw$present==1 & raw$interview==1 & raw$consent==1 & !is.na(raw$consent) & !is.na(raw$interview) & !is.na(raw$present),]
r <- r[!is.na(r$KEY),]
cat('Completed raw submissions linked to final no-SMS IDs:',nrow(r),'\n')
missing <- !r$KEY %in% c(ra$PARENT_KEY,rb$PARENT_KEY)
cat('Of these, raw submissions without knowledge records:',sum(missing),'unique individuals:',length(unique(r$person[missing])),'\n')
stopifnot(nrow(no)==2659,sum(!no$in_know_table)==252,sum(!d$in_know_table)==379,sum(a)==1204,sum(a & k$obs_know_person[j]>0,na.rm=TRUE)==1141)
