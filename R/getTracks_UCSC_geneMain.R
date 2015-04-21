#library('yaml')
library('RMySQL')
library('RCurl')
library('rjson')
library('latticeExtra')

xmlDir<-'http://beta.chromozoom.org/tracks/UCSC/'

### some parameters
track.directory     <-'~/temp/chromozoom/geneMain/'

### access to UCSC mysql database
db.host<-'genome-mysql.cse.ucsc.edu'
db.user<-'genome'
db.pwd <-NULL

### list of our genomes of interest: only the ones with geneMain.xml
selected.genome.name<-gsub('.+"','',gsub('\\.geneMain\\.trakList\\.xml.*','',grep('geneMain.trakList.xml',readLines(xmlDir),value=TRUE)))

##### obtain the main gene track 
for(db.name in selected.genome.name){ #[75:length(selected.genome.name)]

  ### red the table name that has the main gene data from the xml
  track.name<-sub(".+>","",sub("</tableName>","",grep("<tableName",readLines(paste(xmlDir,db.name,'.geneMain.trakList.xml',sep='')),value=TRUE)))
  file.name<-paste(track.directory,db.name,'.',track.name,'.bed',sep='')

  conn <- dbConnect(RMySQL::MySQL(), host=db.host, user=db.user,dbname=db.name)
  my.sql<-paste('select * from ',track.name,sep='')
  track.data<-dbGetQuery(conn,my.sql)
  dbDisconnect(conn)

  ### set score and add the itemRGB, exonSizes, exonStartsZero field
  if(!( 'score' %in% colnames(track.data)))  track.data$score<-0;
  track.data$itemRGB<-'0,0,250'
  track.data$exonSizes<-0
  track.data$exonStartsZero<-''
  for(i in 1:nrow(track.data)){
  	 # simplify if the exonCount (ie blockCount)==1
  	 if(track.data$exonCount[i]==1) {
  	 	track.data$exonSizes[i]<-''; 
  	 } else {
  	   e=track.data$exonEnds[i]  ; en<-as.numeric(unlist(strsplit(e,',')));
       s=track.data$exonStarts[i]; sn<-as.numeric(unlist(strsplit(s,',')));
       track.data$exonSizes[i]<-sub(',$','',paste(en-sn+1,collapse=','))
       track.data$exonStartsZero[i]<-sub(',$','',paste( ifelse( (sn-sn[1]-1)<0,0, sn-sn[1]-1) ,collapse=','))
     }  
  }
  
  if(nrow(track.data)<=1){print(sprintf('Warning! Probably a bigWig file (no content, just location):%s',file.name));next} 
  file.headers<-paste('track name=',track.name,' description=',track.name,' itemRgb="On"',sep='')
  expected.fields  <-c("chrom","txStart"   ,"txEnd"   ,if("name2" %in% colnames(track.data)) "name2" else "name","score","strand","cdsStart"  ,"cdsEnd"  ,'itemRGB',"exonCount" ,"exonSizes" , "exonStartsZero")   ### DB
  destination.cols <-c("chrom","chromStart","chromEnd",                                                   "name","score","strand","thickStart","thickEnd",'itemRGB',"blockCount","blockSizes", "blockStarts")   ### file
  file.headers<-paste(file.headers,paste(paste('# ',destination.cols),collapse="\n"),sep="\n")
  writeLines(file.headers,file.name)
  write.table(track.data[order(track.data$chrom,track.data$txStart),match(expected.fields,colnames(track.data))],file.name,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)

           
}
