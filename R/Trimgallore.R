#Install Trim Gallore
#'InstallTrimGalore
#'@description 
#'Install [Trimgalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) priori to check the availability of [cutapat](/cutadapt.readthedocs.io/en/stable/installation.html)
#'it will install Trimgalore in the Software directory handle by [GenomeDB] library
#'@export
#'@uses
#'/dontrun{
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
#' @param fileR1 the fastq file read 1 (the read 2 will be automatically called)
#' @export
RunTrimgalore <- function(fileR1){
  software <- GenomeDB:::.OpenConfigFile()
  fileR2 <- stringr::str_replace_all(fileR1, "_R1.","_R2.")
  fileR2 <- stringr::str_replace_all(fileR1, "_1.","_2.")
  if(file.exists(software$Software$TrimGalore$command)==FALSE){
    stop("Trimgalore is not installed, please see InstallTtrimgalore()")
  }
  t1<- Sys.time()
  system2(command = software$Software$TrimGalore$command,
          args =c("--paired",paste0("--output_dir ",dirname(fileR1)),fileR1,fileR2))
  etime <- as.numeric(Sys.time()-t1,units="hours")
  of.size1 <-file.info(fileR1)$size
  ot.size1 <-file.info(paste0(unlist(stringr::str_split(fileR1,".f"))[1],"_val_1.fq.gz"))$size
  of.size2 <-file.info(fileR2)$size
  ot.size2 <-file.info(paste0(unlist(stringr::str_split(fileR2,".f"))[1],"_val_2.fq.gz"))$size
  return(c(Time= etime, osize1=of.size1,osize2=of.size1,tsize1=ot.size1,tsize2=ot.size2))
}
