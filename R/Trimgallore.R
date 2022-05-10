#Install Trim Gallore

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
  if(all(stringr::str_detect(version,"version 0.6.6"))==FALSE){
    
    file.remove(tmp.file)
    stop("TrimGalore  installation FAILED")
  }else{
    cat("TrimGalore  installation SUCCESS")
  }
  file.remove(tmp.file)
  GenomeDB:::.OpenConfigFile(software)
}

