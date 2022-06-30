
#'RunSTAR
#'@param sbjFile string full path to the R1.fastq or or gz file
#'@param species the species you whant to use, it should be stored in your database pls see \code{\link[GenomeDB]{ShowSpecies}}
#'@param version the code indicating the genome+annotation pair
#'@param twoPass "None" or "Basic". If set to "Basic", Annotated junctions will be included in both the 1st and 2nd passes. To run STAR 2-pass mapping for
#'each sample separately. STAR will perform the 1st pass mapping,
#'then it will automatically extract junctions, insert them into the genome index, and, finally, re-map
#'all reads in the 2nd mapping pass. This option can be used with annotations, which can be included
#'either at the run-time or at the genome generation step.
#'@param nThreads integer, number of CPUs to be used
#'@export
#'@return the file name of the resulted bam file with an attribute "assemblyVersion" that shows the reference-annotation verision 
#'see \code{\link[GenomeDB]{GetGenomesInDB}}
#'
RunSTAR <- function(sbjFile, species, version, twoPass = c("None","Basic"), nThreads){
  if(require(GenomeDB)==FALSE){
    stop("Pls you should install GenomeDB package")
  }
  software <- GenomeDB:::.OpenConfigFile()
  available.species <- names(software$GenomesDB)[-1]
  if(is.null(software$Software$STAR$main)){
    stop("STAR not installed")
  }
  if(all(stringr::str_detect(available.species,species)==FALSE) ){
    stop(paste0("\n",species, " NOT in Database\nAvailable specias are :", available.species))
  }

  twoPass <- match.arg(twoPass[1], choices = c("None","Basic"))
  file1 <- sbjFile
  if(stringr::str_detect(file1, "_R1.fastq|_R1.fq")){
    file2 <- stringr::str_replace(file1,"_R1.f","_R2.f")
  }else{
    if(stringr::str_detect(file1, "_1.fastq|_1.fq")){
      file2 <- stringr::str_replace(file1,"_1.f","_2.f")
    }else{
      stop("unrecognizable file, it should be sbj_R1.fastq/.gz or sbj_1.fastq/.gz or sbj_R1.fq/.gz or sbj_1.fq/.gz")
    }
  }
  
  
  cat(paste0("\nSubject files\n",file1,"\n",file2,"\n"))
  if(!all(file.exists(file1,file2))){
    stop("ERROR some files not found")
  }

  if(all(c(version %in% names(software$GenomesDB[[species]]$version))==FALSE)){
    stop(paste0("\nThe ", version, "NOT found"))
  }


  out.file.prefix <- file.path(dirname(file1),paste0(basename(dirname(file1)),software$Software$STAR$alignmentPrefix))
  out.file <- paste0(out.file.prefix,"Aligned.out.bam")
  if(missing(nThreads)){
    nThreads <- max(1,parallel::detectCores()-2)
  }

  if(nThreads < 4){
    message(paste0("STAR running on ",nThreads, " may be not optimal"))
  }else{
    message(paste0("STAR running on ",nThreads))
  }

  genomeDir <- software$Software$STAR$main[[version]]

  fasta.gtf.files <- GenomeDB::GetGenome(species, version)

  log <- system2(command = software$Software$STAR$command,
          args = c(paste0("--runThreadN ",nThreads),
                   paste0("--genomeDir " ,genomeDir),
                   paste0("--readFilesIn ",file1, " ", file2),
                   paste0("--readFilesCommand ", ifelse(stringr::str_detect(file1,".gz"),"zcat","-")),
                   # "--outStd BAM_Unsorted",
                   "--outSAMtype BAM Unsorted",
                   "--outSAMunmapped Within",
                   "--outBAMcompression 0",
                   "--outFilterMultimapNmax 50",
                   "--peOverlapNbasesMin 10",
                   "--alignSplicedMateMapLminOverLmate 0.5",
                   "--alignSJstitchMismatchNmax 5 -1 5 5",
                   "--chimSegmentMin 10",
                   "--chimOutType WithinBAM HardClip",
                   "--chimJunctionOverhangMin 10",
                   "--chimScoreDropMax 30",
                   "--chimScoreJunctionNonGTAG 0",
                   "--chimScoreSeparation 1",
                   "--chimSegmentReadGapMax 3",
                   "--chimMultimapNmax 50",
                   paste0("--outFileNamePrefix ",out.file.prefix),
                   paste0("--twopassMode ", twoPass),
                   paste0("--sjdbGTFfile ",fasta.gtf.files$gtf)
          ), stdout = TRUE)#out.file)

  if(!file.exists(out.file)){
    stop("ERROR aligment")
  }
  attr(out.file,"assemblyVersion") <- version
  return(out.file)
}

#' RunSubjunct
#' a wrapper to the \code{\link[Rsubread]{subjuct}} function
#'@param sbjFile string full path to the R1.fastq or or gz file
#'@param species the species you whant to use, it should be stored in your database pls see \code{\link[GenomeDB]{ShowSpecies}}
#'@param version the code indicating the genome+annotation pair
#'@param nThreads integer, number of CPUs to be used
#'@export
#'@return the file name of the resulted bam file with an attribute "assemblyVersion" that shows the reference-annotation verision 
#'see \code{\link[GenomeDB]{GetGenomesInDB}}
#'
RunSubjunct <- function(sbjFile, species, version, nThreads){
  if(require(GenomeDB)==FALSE){
    stop("Pls you should install GenomeDB package")
  }
  software <- GenomeDB:::.OpenConfigFile()
  available.species <- names(software$GenomesDB)[-1]
  
  if(all(stringr::str_detect(available.species,species)==FALSE) ){
    stop(paste0("\n",species, " NOT in Database\nAvailable specias are :", available.species))
  }
  
  file1 <- sbjFile
  if(stringr::str_detect(file1, "_R1.fastq|_R1.fq")){
    file2 <- stringr::str_replace(file1,"_R1.f","_R2.f")
    outf <- stringr::str_remove_all(file1,"_R1.fastq|_R1.fastq.gz|_R1.fq|_R1.fq.gz")
    out.file <- paste0(outf,"_Rsubread_alignment.bam")
  }else{
    if(stringr::str_detect(file1, "_1.fastq|_1.fq")){
      file2 <- stringr::str_replace(file1,"_1.f","_2.f")
      outf <- stringr::str_remove_all(file1,"_1.fastq|_1.fastq.gz|_1.fq|_1.fq.gz")
      out.file <- paste0(outf,"_Rsubread_alignment.bam")
    }else{
      stop("unrecognizable file, it should be sbj_R1.fastq/.gz or sbj_R1.fq/.gz or sbj_1.fastq/.gz or sbj_1.fq/.gz")
    }
  }
  
  cat(paste0("\nSubject files\n",file1,"\n",file2,"\n"))
  if(!all(file.exists(file1,file2))){
    stop("ERROR some files not found")
  }
  
  if(all(c(version %in% names(software$GenomesDB[[species]]$version))==FALSE)){
    stop(paste0("\nThe ", version, "NOT found"))
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
  
  
  
  stats <- Rsubread::subjunc(index = genomeDir,
                    readfile1 = file1,
                    readfile2 = file2,
                    nthreads = nThreads,
                    reportAllJunctions = TRUE,
                    output_file = out.file)
  if(!file.exists(out.file)){
    stop(paste0("ERROR: file not created " ,out.file ))
  }
  attr(out.file,"assemblyVersion") <- version
  return(out.file)
}
