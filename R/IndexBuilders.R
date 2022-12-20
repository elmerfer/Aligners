##Aligners
UpdateSTARindex <- function(){
  if(require(GenomeDB)==FALSE){
    stop("Pls you should install GenomeDB package")
  }
  software <- GenomeDB:::.OpenConfigFile()
  species <- names(software$GenomesDB)[-1]

    if(is.null(software$Software$STAR$main)){
    stop("STAR not installed")
    }
  cat("\nBuilding STAR 2.7.10b - ARRIBA Genome Index...This may take some time")
  thr <- ifelse(parallel::detectCores() > 3, parallel::detectCores()-1, parallel::detectCores())
  cat(paste0("\nRunning STAR genomeGenerate with ",thr," CPU cores"))
  for(sp in species){

    versiones <- basename(software$GenomesDB[[sp]]$version)
    for(vv in versiones){
      if(!(vv %in% names(software$Software$STAR$main))){
        ##check if the genome version exists
        if(dir.exists(file.path(software$Software$STAR$main[1],vv))==FALSE){
          if(dir.create(file.path(software$Software$STAR$main[1],vv))==FALSE){
            stop("error creating index directory")
          }
          ##se creó el directorio
          software$Software$STAR$main[[vv]] <- file.path(software$Software$STAR$main[1],vv)
          fasta.gtf.files <- GenomeDB::GetGenome(sp, vv)
          genome.fasta <- file.path(software$GenomesDB[[sp]],paste0(vv,"/"))
          system2(command = software$Software$STAR$command,
                  args = c("--runMode genomeGenerate",
                           paste0("--genomeDir ",software$Software$STAR$main[[vv]]),
                           paste0("--genomeFastaFiles ",fasta.gtf.files$fasta),
                           paste0("--sjdbGTFfile ",fasta.gtf.files$gtf),
                           paste0("--runThreadN ",thr),
                           "--sjdbOverhang 250") )  
      }
      
      }else{
        message(paste0("\nThe species:", sp," and version: ",vv, "has already being created for STAR aligner"))
      }
      
    }
  }
  GenomeDB:::.OpenConfigFile(software)

}

#' UpdateRsubreadindex
#' It creates the index files for \code{\link[Rsubread]{buildindex}} for each species in the database of \code{\link[GenomeDB]{SetDBdirectory}}
#' @export
UpdateRsubreadindex <- function(){
  if(require(GenomeDB)==FALSE){
    stop("Pls you should install GenomeDB package")
  }
  software <- GenomeDB:::.OpenConfigFile()
  if(is.null(software$Software$Rsubread)){
    stop("Pls run InstallAligners() first to set up directories")
  }
  species <- names(software$GenomesDB)[-1]


  for(sp in species){
    versiones <- names(software$GenomesDB[[sp]]$version)
    for(vv in versiones){
      if(!(vv %in% names(software$Software$Rsubread$main))){
        ##no esta aun el indice
        if( dir.exists(file.path(software$Software$Rsubread$main[1],vv))==FALSE){
          if(dir.create(file.path(software$Software$Rsubread$main[1],vv))==FALSE){
            stop("error creating Rsubread index directory")
          }
        }
        ##se creo el directorio
        software$Software$Rsubread$main <- c(software$Software$Rsubread$main,
                                             file.path(file.path(software$Software$Rsubread$main[1],vv)))
        names(software$Software$Rsubread$main)[length(software$Software$Rsubread$main)] <- vv
        
        fasta.gtf.files <- GenomeDB::GetGenome(sp, vv)
        
        Rsubread::buildindex(basename = software$Software$Rsubread$main[[vv]],
                             reference = fasta.gtf.files$fasta)
      }else{
        message(paste0("\nThe species:", sp," and version: ",vv, "has already being created for Rsubread"))
      }
      

    }
  }
  GenomeDB:::.OpenConfigFile(software)
}

