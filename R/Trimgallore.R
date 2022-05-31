#Install Trim Gallore
#'InstallTrimGalore
#'@description 
#'Install [Trimgalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) priori to check the availability of [cutapat](/cutadapt.readthedocs.io/en/stable/installation.html)
#'it will install Trimgalore in the Software directory handle by [GenomeDB] library
#'@export
#'@usage
#'\dontrun{
#'InstallTrimGalore
#'}
InstallTrimGalore <- function(){
  cut.path <- unlist(stringr::str_split(system2("whereis","cutadapt",stdout=T),":"))
  if(cut.path[2]==""){
    stop("Please verify if cutadapt seems not to be installed, please visit `https://cutadapt.readthedocs.io/en/stable/installation.html`")
  }
  
  if(require(GenomeDB)==FALSE){
    stop("Pls you should install GenomeDB package")
  }
  software <- GenomeDB:::.OpenConfigFile()
  if(length(software)<1){
    stop("\nPlease install GenomeDB package See https://github.com/elmerfer/GenomeDB")
  }
  
  if(is.null(software$Software$main)){
    software$Software$main <- file.path(software$main,"Softwares")
    stopifnot(dir.create(software$Software$main))
    # software$Software$TrimGalore$main <-  file.path(software$Software$main,"Trimgalore")#to install software versions and index files
    # stopifnot(dir.create(software$Software$STrimGalore$main))
  }
  
  
  tmp.file <- tempfile()
  download.file(url="https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz",method = "wget",
                destfile = tmp.file, extra = "--quiet")
  if(!file.exists(tmp.file)){
    stop("Trim Galore download failed")
  }
  zfiles <- untar(tmp.file, list=T)
  untar(tarfile = tmp.file,
        exdir = software$Software$main)
  
  software$Software$TrimGalore$main <- file.path(software$Software$main, basename(zfiles[1]))
  software$Software$TrimGalore$command <- file.path(software$Software$main, zfiles[stringr::str_detect(zfiles,"/trim_galore")])
  version<-system2(software$Software$TrimGalore$command,"--version", stdout = T)
  if(all(stringr::str_detect(version,"version 0.6.6")==FALSE)==TRUE){
    
    file.remove(tmp.file)
    stop("TrimGalore  installation FAILED")
  }else{
    cat("TrimGalore  installation SUCCESS")
  }
  file.remove(tmp.file)
  GenomeDB:::.OpenConfigFile(software)
}

#' RunTrimgalore
#' RunTrimgalore(fileR1)
#' @description RunTrimgalores is intended to perform adapters and trimming of read from paired end sequencing data
#' @param fileR1 the fastq file read 1 (the read 2 will be automatically called) accepted files are XXX_1.fastq/_r1.fastq or fq.gz. Their paired files should replace the _1 (_R1) by _2 (_R2)
#' @outdir if missing the current directory will be used. It will be created if does not exists
#' @param replace (boolean) If TRUE the original files will be reomved and replaced by the trimmed ones under the same original name but gziped (.fq.gz), otherwise it will created two new files 
#' as XXX_1_val_1.fq.gz and XXX_2_val_2,fq.gz
#' @export
#' @usage 
#' \dontrun{
#' file1 <- "my.file_1.fastq #it is expected to have their pair as ny.file_2.fastq
#' res <- RunTrimgalore(file1)
#' res ##this will hold the file1 name and the elapsed time as well as the file sizes prior and after trimming
#' }
RunTrimgalore <- function(fileR1, outdir, replace=FALSE){
  software <- GenomeDB:::.OpenConfigFile()
  fileR2 <- stringr::str_replace_all(fileR1, "_R1.","_R2.")
  fileR2 <- stringr::str_replace_all(fileR1, "_1.","_2.")
  if(file.exists(software$Software$TrimGalore$command)==FALSE){
    stop("Trimgalore is not installed, please see InstallTtrimgalore()")
  }
  
  gziped <- ifelse(stringr::str_detect(fileR1,".gz"),"--gzip","--dont_gzip")
  outdir <- ifelse(missing(outdir), paste0("--output_dir ",dirname(fileR1)),paste0("--output_dir ",file.path(dirname(fileR1),outdir)))
  
  t1<- Sys.time()
  system2(command = software$Software$TrimGalore$command,
          args =c("--paired",
                  gziped,
                  outdir,
                  fileR1,
                  fileR2))
  
  etime <- as.numeric(Sys.time()-t1,units="hours")
  
  of.size1 <-file.info(fileR1)$size
  of.size2 <-file.info(fileR2)$size
  if(gziped=="--gzip"){
    ofile <- paste0(unlist(stringr::str_split(basename(fileR1),".f"))[1],"_val_1.fq.gz")
    file.rename(ofile,stringr::str_replace_all(ofile,"_1_|_R1_",""))
    ofile <-stringr::str_replace_all(ofile,"_1_|_R1_","")
    ot.size1 <-file.info(file.path(outdir,ofile))$size  
    of2 <- paste0(unlist(stringr::str_split(basename(fileR2),".f"))[1],"_val_2.fq.gz")
    file.rename(of2,stringr::str_replace_all(of2,"_2_|_R2_",""))
    ot.size2 <-file.info(file.path(outdir,of2))$size
  }else{
    ofile <- paste0(unlist(stringr::str_split(basename(fileR1),".f"))[1],"_val_1.fastq")
    file.rename(ofile,stringr::str_replace_all(ofile,"_1_|_R1_",""))
    ofile <-stringr::str_replace_all(ofile,"_1_|_R1_","")
    ot.size1 <-file.info(file.path(outdir,paste0(unlist(stringr::str_split(basename(fileR1),".f"))[1],"_val_1.fastq")))$size  
    of2 <- paste0(unlist(stringr::str_split(basename(fileR2),".f"))[1],"_val_2.fastq")
    file.rename(of2,stringr::str_replace_all(of2,"_2_|_R2_",""))
    ot.size2 <-file.info(file.path(outdir,of2))$size    
  }
  
  
  
  
  return(c(ofile=fileR1,tfile = ofile ,Time= etime, osize1=of.size1,osize2=of.size1,tsize1=ot.size1,tsize2=ot.size2))
}
