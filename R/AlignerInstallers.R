
#' InstallAligners
#' This function will install and configure the use fo the STAR sofware version 2.7.6a compatible with the Arriba software, and the Rsubread bioconductor package
#' @export
#' it will build the a directory based on \code{\link[GenomeDB]{SetDbdirectory}} function
InstallAligners <- function(){
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
    software$Software$STAR$main <-  file.path(software$Software$main,"STAR")#to install softwar versions and index files
    stopifnot(dir.create(software$Software$STAR$main))
    software$Software$Rsubread$main <-  file.path(software$Software$main,"Rsubread")  #to install index files
    stopifnot(dir.create(software$Software$Rsubread$main))
  }
  if(require("Rsubread")==FALSE){
    message("\nInstalling Rsubread")
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")

    BiocManager::install("Rsubread")
  }

  cat("\nDownloading the STAR aligner software version 2.7.6a compatible with Arriba version 2.0.1")
  tmp.destfile <- tempfile()
  download.file(url = "https://github.com/alexdobin/STAR/archive/2.7.6a.tar.gz",
                method = "wget",
                destfile = tmp.destfile, extra = "--quiet")

  files <- untar(tarfile = tmp.destfile,list = T )
  star.f <- files[which(stringr::str_detect(files, "Linux"))]
  star.f[stringr::str_detect(star.f, "static")]
  str.comp <- which(unlist(stringr::str_split(star.f[stringr::str_detect(star.f, "static")][2],"/") )=="STAR")-1

  software$Software$STAR$command <- file.path(software$Software$STAR$main,"STAR")
  software$Software$STAR$version <- "2.7.6a"

  untar(tarfile = tmp.destfile,
        files = star.f[stringr::str_detect(star.f, "static")][-1],
        exdir = software$Software$STAR$main,
        extras = paste0("--strip-components ",str.comp))
  file.remove(tmp.destfile)

  if(file.exists(software$Software$STAR$command)){
    cat("\nSTAR 2.7.6a installed")
  }
  system2(command = software$Software$STAR$command)

  software$Software$STAR$alignmentPrefix <- "_STAR_Aligned_"
  software$Software$STAR$alignmentPrefixSorted <- "_STAR_AlignedSortedByCoordinates_"

  GenomeDB:::.OpenConfigFile(config = software)
}
