stopifnot(getRversion()=='4.5.3')
db<-read.delim('/tmp/r-packages.tsv',stringsAsFactors=FALSE)
# All sources are local, exact-version archives. Retry only installation order,
# never alternate package versions; the final version check is mandatory.
remaining<-db$package
for(pass in seq_len(nrow(db))){
 before<-remaining
 for(p in before){
  file<-file.path('/tmp/r-sources',paste0(p,'_',db$version[db$package==p],'.tar.gz'))
  desc<-tempfile();dir.create(desc);untar(file,files=paste0(p,'/DESCRIPTION'),exdir=desc)
  d<-read.dcf(file.path(desc,p,'DESCRIPTION'))
  fields<-intersect(c('Depends','Imports','LinkingTo'),colnames(d))
  deps<-unique(trimws(gsub('\\s*\\([^)]*\\)','',unlist(strsplit(paste(d[1,fields],collapse=','),',')))))
  deps<-setdiff(deps,c('R',''))
  installed<-rownames(installed.packages())
  if(!all(deps%in%installed))next
  install.packages(file,repos=NULL,type='source',Ncpus=2)
  stopifnot(requireNamespace(p,quietly=TRUE),
            utils::packageDescription(p,fields='Version')==db$version[db$package==p])
  remaining<-setdiff(remaining,p)
 }
 if(!length(remaining))break
 stopifnot(length(remaining)<length(before))
}
stopifnot(!length(remaining))
writeLines(capture.output(sessionInfo()),'/opt/derived-r-session.txt')
