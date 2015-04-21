library('RMySQL')
### some parameters
track.directory    <-'/Users/pasculescu/Sites/chromozoom/tracks/UCSC/'  # ehere the bed files are
#track.directory    <-'/Users/administrator/temp/chromozoom/geneMain/'  # ehere the bed files are
bin.directory      <-'/Users/pasculescu/Sites/bin/'  # ehere the bin files are
chromsize.file.name <-paste(track.directory,'chrom.sizes',sep='')
sorted.bed.file.name<-paste(track.directory,'sorted.bed',sep='')

### access to UCSC mysql database
db.host<-'genome-mysql.cse.ucsc.edu'
db.user<-'genome'
db.pwd <-NULL

bed.files<-list.files(track.directory,pattern='\\.bed$',full.names=TRUE)

##### obtain the main gene track 
for(f in bed.files[1]){

  db.name<-sub('.+/','',sub('\\..+','',f)) 

  conn <- dbConnect(RMySQL::MySQL(), host=db.host, user=db.user,dbname=db.name)
  my.sql<-paste('select chrom,size from chromInfo',sep='')
  chrom.data<-dbGetQuery(conn,my.sql)
  dbDisconnect(conn)
  if(nrow(chrom.data)<=1){print(sprintf('Warning! missing chromInfo table:%s',db.name));next} 

  ### create a chromsize file
  write.table(chrom.data[order(chrom.data[,1]),],chromsize.file.name,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
  f.out<-sub('\\.bed$','.bb',f)

  #sort `sort -k1,1 -k2,2n`
  #and get rid of comment rows and the '^track' header and fill in the length, (end-start) and 0, for genes with one exon only !
  bed<-readLines(f)
  bed<-bed[grep('^#',bed,invert=TRUE)]
  bed<-bed[grep('^track',bed,invert=TRUE)]
  filter<-grep('\t\t$',bed)
  bed[filter]<-sapply(bed[filter],function(x) {e<-as.numeric(unlist(strsplit(x,'\t'))[2:3]); sub('\t$',paste(e[2]-e[1],',\t0,',sep=''),x)})
  order.1<-sub('\t.+','',bed)
  order.2<-as.numeric(sub('^[^\t]+\t([^\t]+).+','\\1',bed))
  bed<-bed[order(order.1,order.2)]  ### ??? does not work
  write(bed,sorted.bed.file.name)
  
  
  cmd<-paste(bin.directory,'bedToBigBed ',sorted.bed.file.name,' ',chromsize.file.name,' ',f.out,sep='' )
  ### error: Error line 245 of /Users/pasculescu/Sites/chromozoom/tracks/UCSC/sorted_2.bed: BED blocks must be in ascending order without overlap. Blocks 0 and 1 overlap.

         
}
