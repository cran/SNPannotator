#' Run the annotation pipeline on a list of variants
#'
#' This function receives a list of variants and checks their information on Ensembl website via
#' the Ensembl REST API server.
#'
#' @param rslist A vector of rs numbers.
#' @param server Name of the server. "https://rest.ensembl.org" can be used for GRCh38
#' and "https://grch37.rest.ensembl.org" for GRCh37.
#' @param db The population database for calculating LD scores. This can be found using `listDatabases` function.
#' @param outputPath The report file will be saved in this path as an Excel file (.xlsx)
#' @param window_size Number of base pairs around the variant for checking LD scores (max = 500kb)
#' @param r2 The LD threshold for selecting variants around the target SNP.
#' @param LDlist If set to TRUE, variants in high LD will be found and added to the output.
#' @param cadd If set to TRUE, the CADD scores will be added to variant information.
#' @param geneNames.file path the gene information file (*.rds). Default value is NULL and ENSEMBL website will be checked if no file is provided.
#' @param regulatoryType.file path the variants regulatory type information file (*.rds). Default value is NULL and this step will be skipped if no file is provided.
#' @param cores set to a value above 0 for parallel processing.
#' @return a data table with all variant information is returned.
#' @examples
#' \dontrun{
#' # select the required server
#' server <- "https://grch37.rest.ensembl.org"
#'
#' # select the database for population data
#' # this can be selected from listDatabases() function
#' db <- "1000GENOMES:phase_3:EUR"
#'
#' # create a vector of required SNPs
#' rslist=c('rs236349')
#'
#' output <- annotate(rslist,server,db,
#'   outputPath = paste(tempdir(),'sampleOutput.xlsx',sep="/"),
#'   window_size = 500,
#'   r2 = .9,
#'   cadd = FALSE)
#'   }
#'
annotate <- function(rslist,
                              server,
                              db,
                              outputPath,
                              window_size=500,
                              r2=0.5,
                              LDlist = TRUE,
                              cadd = FALSE,
                              geneNames.file=NULL,
                              regulatoryType.file=NULL,
                              cores=0) {


  if(file.exists(outputPath))
  {
    outputPath = sprintf('%s_%s',createRandString(),outputPath)
  }


  stopifnot('RS list must be a vector.' = is.vector(rslist),
            'Range for r2 is 0-1' = (r2 >0 & r2<=1),
            'Range for window size is 100-500' = (window_size > 100 & window_size<=500)
  )

  # check if system has the number of selected cores
  check.core.count(cores)


  # load gene name file - to add gene information for intergenic variations
  # based on grch37 or 38
  g.set <- NULL
  get.gene.info.from.API <- TRUE

  if(!is.null(geneNames.file))
  {
    g.set <- getGeneFile.external(geneNames.file)
    get.gene.info.from.API <- FALSE
  }



  # load regulatory file
  # based on grch37
  r.set <- NULL
  if(!is.null(regulatoryType.file))
    r.set <- getRegulatoryFile.external(regulatoryType.file)

  # load cytoband file
  # based on grch37
  c.set <- getCytobandFile(server)

  # set start time
  start.time <-  proc.time()

  on.exit({
    if(!is.null(start.time))
      cat(sprintf("\nRun time: %s",timetaken(start.time)),fill = TRUE)# END LOG
  })

  ######################
  ######################

  pingServer(server)

  #wb <- createWorkbook()
  output <- data.table()

  cat("Fetching variant information ...", fill = TRUE)

  for(i in 1:length(rslist))
  {
    rs= rslist[i]

    # get variant info + variants in high LD
    varInfo <- tryCatch(
      {
        inSilicoSeqPipeline(rs, server, db , window_size , r2 , i, length(rslist), LDlist, get.gene.info.from.API)
      },
      error = function(err) {
        message(paste0("Error in fetching variant information: ", err$message))
        NULL
      }
    )

    # rsid is checked vs synonym id
    if(!is.null(varInfo$LDList))
      varInfo$LDList <- checkReturnedVariantData(varInfo)

    if(!is.null(varInfo))
    {

      tab <- returnVariantDatatable(i, varInfo,db)

      # fill missing genes
      if(!is.null(g.set))
      {
        cat("Checking gene information ...", fill = TRUE)

        tab[Gene=='',
            c('GeneId','Gene') := find.nearest.gene(Chr,Pos_37,g.set),
            by=list(Chr,Pos_37)]

        tab[Cytoband=='',
            Cytoband := find.band(Chr,Pos_37,c.set),
            by=list(Chr,Pos_37)]
      }
      # add regulatory features

      if(!is.null(r.set))
      {

        cat("Checking regulatory feature information ...", fill = TRUE)

        tab[,c('RegId','RegType') := find.regulatory(Chr,Pos_37,r.set),
            by=list(Chr,Pos_37)]
      }

      ## check if this variants is already present in the output
      ## maybe the provided variants are not independent with each other
      checkifOutputIncludesRS(rs,output)

      output <- rbind(output,tab)

      # add a sheet for each variant
      # appendXLSXfile(tab,sheetName = rs,fileName = outputPath)
    }
    else
    {
      cat("Variant not found.\n")
    }
  }
  cat("done\n", fill = TRUE)


  if(cadd & nrow(output) > 0)
  {
    # cat("=============================================\n")
    cat('Fetching CADD scores ... ', fill = TRUE)
    # cat("=============================================\n")

    ## for progress bar 2
    pb <- progress::progress_bar$new(format = "[:bar] :current/:total (:percent)", total = length(output$Linked_SNP))
    pb$tick(0)
    ##

    conseq.dt <- tryCatch(
      {
        if(cores == 0)
          getConseqOnList(output$Linked_SNP,server,pb)
        else
          getConseqOnList.parallel(output$Linked_SNP,server,cores)
      },
      error =function(err)
      {
        message(paste0("Error in fetching CADD information: ",err$message))
        NULL
      }
    )

    if(!is.null(conseq.dt))
    {
     # output = merge(x=output,y=conseq.dt,by.x = 'Linked_SNP',by.y = 'rs',all.x = TRUE,sort = FALSE)
      output[is.element(Linked_SNP,conseq.dt$rs) , cadd := conseq.dt[conseq.dt$rs == Linked_SNP]$cadd ]

      if(any(grepl(x = output$cadd,pattern = ',')))
        output[grepl(x = cadd,pattern = ','), cadd := select.CADD.scores(Alleles,cadd),by = Linked_SNP]
    }

  }



  # save excel file
  if(nrow(output) > 0)
    appendXLSXfile(output,sheetName = 'All Variants',fileName = outputPath, addFirst = FALSE)


  # save html file
  if(rmarkdown::pandoc_available())
  {
    if(nrow(output) > 0 & LDlist== TRUE)
    {
      html.file = gsub(x = outputPath, pattern = 'xlsx',replacement = 'html')

      if(generate.report.file(output, html.file,server,db,window_size,r2))
        cat(sprintf('File saved: %s',html.file),fill = TRUE)
    }
  }else {
    cat('Pandoc is not available. Skipping HTML report.',fill = TRUE)
  }



  invisible(output)
}
