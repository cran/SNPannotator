# load gene name file

# gene data file had to be put outside the package due to size limit from CRAN
getGeneFile <- function(server)
{

  grch37 <- grepl(pattern = '37',x = server)

  if(grch37)
    gene.File <- system.file("extdata", "Gene_Names_Ensembl_104_GRCh37.rds", package = "SNPannotator")
  else
    gene.File <- system.file("extdata", "Gene_Names_Ensembl_104_GRCh38.rds", package = "SNPannotator")


  if (!file.exists(gene.File))
    stop('Gene file not found in the package!', call. = FALSE)
  else
    return(readRDS(gene.File))
}

getGeneFile.external <- function(geneFilePath)
{

  gene.File <- NULL

  if(file.exists(geneFilePath))
    gene.File <- tryCatch({
      readRDS(geneFilePath)
    },
    error = function(er){
      message(paste0('Error reading file. (',geneFilePath,')',er$message))
    }
    )
  else
    message(paste0('File not found! (',geneFilePath,')'))

  return(gene.File)

}

# load regulatory file
# regulatory data file had to be put outside the package due to size limit from CRAN
getRegulatoryFile <- function(server)
{

  reg.File <- system.file("extdata", "homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20201218.rds", package = "SNPannotator")


  if (!file.exists(reg.File))
    stop('Regulatory file not found in the package!', call. = FALSE)
  else
    return(readRDS(reg.File))
}

getRegulatoryFile.external <- function(regFilePath)
{
  reg.File <- NULL

  if(file.exists(regFilePath))
    reg.File <- tryCatch({
      readRDS(regFilePath)
    },
    error = function(er){
      message(paste0('Error reading file. (',regFilePath,')',er$message))
    }
    )
  else
    message(paste0('File not found! (',regFilePath,')'))


  return(reg.File)
}

# load cytoband file
getCytobandFile <- function(server)
{

  c.File <- system.file("extdata", "cytoband_GRCh37.rds", package = "SNPannotator")


  if (!file.exists(c.File))
    stop('Cytoband information file not found in the package!', call. = FALSE)
  else
    return(readRDS(c.File))
}
