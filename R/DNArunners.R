#' RunAlign
#' a wrapper to the \code{\link[Rsubread]{align}} function
#'@param sbjFile string full path to the R1.fastq or or gz file
#'@param species the species you whant to use, it should be stored in your database pls see \code{\link[GenomeDB]{ShowSpecies}}
#'@param version the code indicating the genome+annotation pair
#'@param nThreads integer, number of CPUs to be used
#'@export
#'
RunAlign <- function(sbjFile, species, version, nThreads){
  if(require(GenomeDB)==FALSE){
    stop("Pls you should install GenomeDB package")
  }
  software <- GenomeDB:::.OpenConfigFile()
  available.species <- names(software$GenomesDB)[-1]
  
  if(all(stringr::str_detect(available.species,species)==FALSE) ){
    stop(paste0("\n",species, " NOT in Database\nAvailable specias are :", available.species))
  }
  
  file1 <- sbjFile
  file2 <- .IdentifyRead2(file1)
  
  cat(paste0("\nSubject files\n",file1,"\n",file2,"\n"))
  if(!all(file.exists(file1,file2))){
    stop("ERROR some files not found")
  }
  
  if(all(stringr::str_detect(names(software$GenomesDB)[-1], species)==FALSE)){
    stop(paste0("\nThe ", species, "NOT found"))
  }
  
  
  if(missing(nThreads)){
    nThreads <- max(1,parallel::detectCores()-2)
  }
  
  if(nThreads < 4){
    message(paste0("Subjunct running on ",nThreads, " may be not optimal"))
  }else{
    message(paste0("Subjunct running on ",nThreads))
  }
  
  genomeDir <- software$Software$Rsubread$main[[version]]
  
  
  
  Rsubread::align(index = genomeDir,
                    readfile1 = file1,
                    readfile2 = file2,
                  type = "dna",
                    nthreads = nThreads,
                    detectSV = TRUE,
                  sortReadsByCoordinates = TRUE,
                  annot.ext = software$GenomesDB[[species]]$version[version],
                  isGTF = TRUE)
}


.IdentifyRead2 <- function(file1){
  ##is trimmed?
  if(stringr::str_detect(file1,"1.trim.fastq")){
    if(stringr::str_detect(file1, ".gz")){
      ##esta gzip
      file2 <- stringr::str_replace(file1,"1.trim.fastq.gz","2.trim.fastq.gz")
    }else{
      file2 <- stringr::str_replace(file1,"1.trim.fastq","2.trim.fastq")
    }
  }else{
    if(stringr::str_detect(file1, ".gz")){
      ##esta gzip
      file2 <- stringr::str_replace(file1,"1.fastq.gz","2.fastq.gz")
    }else{
      file2 <- stringr::str_replace(file1,"1.fastq","2.fastq")
    }
  }
  return(file2)
}

