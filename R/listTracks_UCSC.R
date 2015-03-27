#library('yaml')
library('RMySQL')
library('RCurl')
library('rjson')
library('latticeExtra')

ajaxDir<-switch(
Sys.info()['nodename']
                ,'Adrians.local'          ='http://localhost/~pasculescu/chromozoom/php/'
                ,'llama.mshri.on.ca'      ='http://beta.chromozoom.org/php/'
                ,'pawson-imac.local'      ='http://localhost/~pasculescu/chromozoom/php/'  
)

### some parameters
track.directory     <-'~/temp/chromozoom/tracks/'

### access to UCSC mysql database
db.host<-'genome-mysql.cse.ucsc.edu'
db.user<-'genome'
db.pwd <-NULL

### list of our genomes of interest
conn<-paste(ajaxDir,'chromsizes.php',sep='')           
genome.json<-getURL(conn) ### the response is obtained here

##### parse the JSON response in an R structure data.frame
genome.as.list<-fromJSON(genome.json)
genome.data<-Reduce(rbind,Map(data.frame,genome.as.list[lapply(genome.as.list,length)==3])) ## accept only records with all fileds (name,species,assemblyDate)

### connect to  UCSC database
conn <- dbConnect(RMySQL::MySQL(), host=db.host, user=db.user)
#### the full list of accepted  genome names can be obtained by listing all avalable databases
available.genome.name.UCSC<-dbGetQuery(conn,'show databases')
genomes.not.in.UCSC.db<-setdiff(genome.data$name,available.genome.name.UCSC[,1])
print(paste('WARNING! genomes in chromozoom and not in UCSC db:',paste(sort(genomes.not.in.UCSC.db),collapse=', ')))

latest.genome.name   <-tapply(genome.data$name,sub('[0-9]+$','',genome.data$name)
                             ,function(g){ pre<-unique(sub('[0-9]+$','',g)); post<-as.numeric(sub(pre,'',g,fixed=T)); paste(pre,max(post,na.rm=T),sep='') } )
latest.genome.name   <-setdiff(latest.genome.name,genomes.not.in.UCSC.db)
print(paste('Latest genomes UCSC:\n',paste(latest.genome.name,collapse=', ')) )                             
dbDisconnect(conn) ### no need for a connection

#selected.genome.name <-c('dm6','sacCer3','hg19')  
selected.genome.name <-intersect(as.character(genome.data$name),available.genome.name.UCSC[,1]) 

##### obtain information about wanted tracks
for(db.name in selected.genome.name[68:length(selected.genome.name)]){
  if(! (db.name %in% genome.data$name )) {print(paste('ERROR! selected genome not in UCSC list:',db.name))}

  ### connect to the desired database
  conn <- dbConnect(RMySQL::MySQL(), host=db.host, user=db.user,dbname=db.name)

  if(!dbExistsTable(conn,'trackDb')){dbDisconnect(conn); next}
  
  ### select only the track tables in the desired format
  my.sql<-paste('select tableName,shortLabel,`type`,longLabel,html,grp,settings from '
               ,'trackDb where grp="genes" and (`type` like "bed%" or `type` like "big%" or `type` like "wig%" )',sep='')            
  ### execute sql and obtain all table names that contain desired tracks
  track.table<-dbGetQuery(conn,my.sql)

  file.name<-paste(track.directory,db.name,'.genes.trakList.xml',sep='')
  #### now obtain information for each track and store it as a file in the proper format
  ### use a differential from existing and from the ones that are in link list
  xml.content<-''
  track.table<-track.table[gsub(' .*','',track.table$type) %in% c('bed','wig','bigWig'),]
  if(nrow(track.table)<=0) { dbDisconnect(conn); next }
  for(i in 1:nrow(track.table)){
     track.type.info<-track.table$type[i]
     track.type <-gsub(' .*','',track.type.info)
     track.name<-track.table$tableName[i]
     ### check if track type is the desired one
     if(!(track.type %in% c('bed','wig','bigWig'))) next
     
     ### check if table exists
     if(!dbExistsTable(conn,track.name)) next

     my.sql<-paste('select count(*) from ',track.name,sep='')
     track.nmbRows<-dbGetQuery(conn,my.sql)

     if(track.nmbRows<=1) notes=sprintf('Warning! Probably a bigWig file (no content, just location)') else notes=""; 

     xml.row<-paste(
       '<tableName type="',track.table$type[i],'" grp="',track.table$grp[i],'"  nmbRows="',track.nmbRows,'">',track.name,'</tableName>\n'
      ,'<shortLabel><![CDATA[',track.table$shortLabel[i],']]></shortLabel>\n'
      ,'<longLabel><![CDATA[',track.table$longLabel[i],']]></longLabel>\n'
#      ,'<html>',gsub('<br>','<br/>',gsub(' > ',' &gt; ',gsub(' < ',' &lt; ',gsub('=_BLANK','="_blank"',track.table$html[i],ignore.case=TRUE))),ignore.case=TRUE),'</html>\n'
      ,'<html><![CDATA[',track.table$html[i],']]></html>\n'
      ,'<settings><![CDATA[',track.table$settings[i],']]></settings>\n'
      ,'<notes>',notes,'</notes>\n'
      ,sep='')


     xml.content<-paste(xml.content,'<trackInfo>',xml.row,'</trackInfo>\n',sep='')
  }
  xml.content<-paste(
    '<?xml version="1.0" encoding="iso-8859-1" standalone="no"?>\n'
   ,'<?xml-stylesheet type="text/xsl" href="trackList.xsl"?>\n'
   ,'<trackList source="UCSC" organism="',db.name,'" extractionDate="',date(),'">\n'
   ,xml.content
   ,'</trackList>\n'
   ,sep='')
  writeLines(xml.content,file.name)
  
  dbDisconnect(conn)
            
}
