fetch <- function(ext,server)
{
  counter <- 1
  notFound <- TRUE

  while(counter < 6 && notFound)
  {
    #cat(sprintf("Try number %s",counter),fill = TRUE)

    r <- GET(paste(server, ext, sep = ""), content_type("application/json"))

    if(checkStatusCode(r))
      notFound <- FALSE
    else
      counter <- counter + 1

  }

  if(!notFound)
    return(fromJSON(toJSON(content(r))))
  else
  {
    #cat("Not found after 5 tries.",fill = TRUE)
    return(NULL)
  }
}


#' Checks if the service is alive
#'
#' This function test whether the Ensembl server is accessible or not
#'
#' @param server name of the server. "https://rest.ensembl.org" can be used for GRCh38
#' and "https://grch37.rest.ensembl.org" for GRCh37.
#' @return a message is displayed to the user
#' @examples
#' # select the required Ensembl server
#' server = "https://rest.ensembl.org"
#'
#' # check if server is accessible
#' pingServer(server)
#'
pingServer <- function(server)
{
  tryCatch(
    {

      ext <- "/info/ping?"
      r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
      stop_for_status(r)

      response <- fromJSON(toJSON(content(r)))

      cat("\nServer ping test ... ")
      if(response$ping == 1)
        cat("accessible\n\n")
      else
        stop("Not accessible\n")

    },
    error = function(err) stop(err$message,call. = FALSE)
  )
}

#' List population from human database (1000 Genomes project)
#'
#' This function list the name, description and size of the available populations
#' in 1000 Genomes project database. This database will be used for returning variables in high LD
#' with the target SNP.
#'
#' @param server name of the server. "https://rest.ensembl.org" can be used for GRCh38
#' and "https://grch37.rest.ensembl.org" for GRCh37.
#' @return A data table is returned which includes the name, description and size of the available populations
#' in 1000 Genomes project database.
#' @examples
#' # select the required Ensembl server
#' server = "https://rest.ensembl.org"
#'
#' # check the available population data for the selected server
#' listDatabases("https://rest.ensembl.org")
#'
listDatabases <- function(server)
{
  ext <- "/info/variation/populations/human?"
  r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
  stop_for_status(r)

  response <- fromJSON(toJSON(content(r)))

  setDT(response)

  return(response[grepl(x=name,pattern = '1000'),list(name,size,description)])

}

#' Shows the data releases available on this REST server.
#'
#' Shows the data releases available on this REST server.
#' May return more than one release (unfrequent non-standard Ensembl configuration).
#'
#' @param server name of the server.
#' @return a message is displayed to the user
#' @examples
#' # select the required Ensembl server
#' server = "https://rest.ensembl.org"
#'
#' # check the data releases of the selected server
#' releaseVersion(server)
#'
releaseVersion <- function(server)
{
  ext <- "/info/data?"
  r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
  stop_for_status(r)

  return( fromJSON(toJSON(content(r))))

}

getVariantInfo <- function(rsID,server,pb=NULL,get.gene.info.from.API)
{
  on.exit({
    if(!is.null(pb))
      pb$tick(1)
  })

  ext <- sprintf("/variation/human/%s?phenotypes=1;pops=1",rsID)
  #cat(sprintf("Fetching variant info: ... "))

  var <- fetch(ext,server)

  if(!is.null(var))
  {
    if(length(var$mappings) == 0)
      #  cat("done",fill = TRUE)
      # else
    {
      cat(sprintf("no data found"))

      return(NULL)
    }

    if(get.gene.info.from.API)
      var <- getVariantGene(var, server)

    return(var)
  }
  else
    return(NULL)

}

getVariantConsequence <- function(rsID,server)
{
  ext <- sprintf("/vep/human/id/%s?CADD=1",rsID)
  # cat(sprintf(" Fetching variant consequence: ... "))

  response <- fetch(ext,server)

  if(!is.null(response))
  {


    # if(length(response$most_severe_consequence) == 0)
    #   cat("done",fill = TRUE)
    # else
    # {
    #   cat(sprintf("%s- no data found",rsID))
    #   return(NULL)
    # }

    if(!is.null(response$transcript_consequences))
      dt= list(most_severe_consequence=response$most_severe_consequence,
               consequences = response$transcript_consequences)
    else
      dt= list(most_severe_consequence=response$most_severe_consequence,
               consequences = response$intergenic_consequences)

    return(dt)
  }
  else
    return(NULL)

}


getVariantGene <- function(varInfo,server)
{
  ext <- sprintf("/overlap/region/human/%s?feature=gene;feature=band",varInfo$mappings$location[[1]])
  # cat(sprintf("Fetching gene info: ... "))

  response <- fetch(ext,server)

  if(!is.null(response))
  {
    #cat(sprintf("done (%s).", varInfo$mappings$location[[1]]), fill = TRUE)
    # cat("done",fill = TRUE)

    varInfo$geneInfo <- response
  }

  return(varInfo)
}


checkStatusCode <- function(response)
{
  code <- response$status_code

  if(code == 200L)
    return(TRUE)
  else if(code == 400L)
  {
    #cat('Bad Request.', fill=TRUE)
    return(FALSE)
  }
  else if (code== 403L) {
    cat('You are submitting far too many requests.', fill=TRUE)
    return(FALSE)
  }
  else if (code== 404L)
  {
    cat('Indicates a badly formatted request. Check your URL.', fill=TRUE)
    return(FALSE)
  }
  else if (code== 408L)
  {
    cat('Timeout.', fill=TRUE)
    return(FALSE)
  }
  else if(code == 429L)
  {
    cat('You have been rate-limited.', fill=TRUE)
    return(FALSE)
  }
  else if (code== 503L){
    cat('The service is temporarily down.', fill=TRUE)
    return(FALSE)
  }
  else
  {
    cat('unknown error.')
    return(FALSE)
  }


}

checkRemainingLimit <- function(response)
{
  if(r$headers$`x-ratelimit-remaining` < 1)
    stop(sprintf('your requests are limited. Try again in %s seconds.'),r$headers$`x-ratelimit-reset`)
}

correctSynonymIDs <- function(varInfo)
{
  if(is.null(varInfo$LDList))
    return(list(0))
  LDTable <- varInfo$LDList # table including targetSNP and linkedSNPs , names and LDs
  LDVariants <- varInfo$LDlistFull # list of found data for LD variants
  LDVariantsNames <- sapply(LDVariants,function(x) return(x$name)) # name of found variants in high LD
  # check if the same name or a synonym is used for the found variant info

  LDTable[,checkSynonyms := ifelse(is.element(variation2 ,LDVariantsNames),0,1)]

  if(LDTable[checkSynonyms ==1,.N] > 0)
  {

    LDTable[checkSynonyms == 1,
            variation3 := returnSynonym(variation2,LDVariants), by = variation2]

    # some notifications
    changedIDsTable <- LDTable[checkSynonyms == 1,]
    if(nrow(changedIDsTable) > 0)
    {
      cat("\n== Some synonym rsIDs were changed in the original LD table.",fill = TRUE)
      for(i in seq_len(nrow(changedIDsTable)))
        cat(sprintf("%s - %s changed to %s",
                    i,
                    changedIDsTable[i,variation2],
                    changedIDsTable[i,variation3]),fill = TRUE)
    }

    LDTable[,variation2 := ifelse(is.na(variation3),variation2,variation3)]
  }
  return(LDTable)
}

returnSynonym <- function(rsID, LDVariants)
{
  newID <- rsID

  if(is.null(LDVariants) || length(LDVariants) < 1)
    return(newID)
  else
    for(i in seq_len(length(LDVariants)))
      if(is.element(rsID,LDVariants[[i]]$synonyms))
        newID <- LDVariants[[i]]$name

    return(newID)
}

saveOutputData <- function(varInfoTbl, wb)
{
  # wb <- createWorkbook()
  sheetName=varInfoTbl[1,gSNP]
  addWorksheet(wb, sheetName)
  writeDataTable(wb, sheetName, x=varInfoTbl, tableStyle = "TableStyleMedium9")
  freezePane(wb, sheetName, firstRow = TRUE)
  setColWidths(wb, sheetName, cols = 1:ncol(varInfoTbl), widths = "auto")
  #saveWorkbook(wb, path, overwrite = overwriteExistingFile)
  return(wb)
}

appendXLSXfile <- function(output,sheetName,fileName,addFirst = FALSE)
{
  if(!file.exists(fileName))
  {
    write.xlsx(x = output,file = fileName,sheetName= sheetName)
  }
  else
  {

    tryCatch({

      wb <- loadWorkbook(fileName)
      addWorksheet(wb = wb, sheetName = sheetName)
      writeData(wb = wb, sheet = sheetName, x = output, colNames = TRUE, rowNames = FALSE)

      if(addFirst)
      {
        sheet.count <- length(openxlsx::getSheetNames(fileName))
        worksheetOrder(wb) <- c(sheet.count+1,seq(sheet.count))
      }

      saveWorkbook(wb = wb, file = fileName, overwrite = TRUE)

    },
    error = function(err){
      print(err)
    })
  }
  cat(sprintf('File saved: %s',fileName),fill = TRUE)
}


getConseqOnList <- function(rslist, server, pb) {

  out.dt <- data.table(rs=character(),cadd=character())

  for(i in seq_len(length(rslist)))
  {
    rs= rslist[i]
    #cat(i , "-", rs , '... ')

    # store the returned scores
    outString <- ""

    conseq <- getVariantConsequence(rs,server)

    if(all(is.element(c('cadd_phred','variant_allele') , names(conseq$consequences[[1]]))))
    {
      conseq.dt <- unique.data.frame(data.table(conseq$consequences[[1]][,'variant_allele'],
                                                conseq$consequences[[1]][,'cadd_phred']))
      if(length(conseq.dt$V2) > 1)
        conseq.dt <- sapply(conseq.dt, function(x) ifelse(x == "NULL", NA, x))

      conseq.dt <- data.table(matrix(unlist(conseq.dt),nrow = nrow(conseq.dt),byrow = FALSE))
      outString <-  paste(paste(conseq.dt$V1,collapse = ","),
                          paste(conseq.dt$V2,collapse = ","),sep=" = ")
    }

    dt <- data.table(rs=rs,cadd=outString)
    out.dt <- rbind(out.dt,dt)

    pb$tick(1)
  }
  cat("done\n",fill = TRUE)
  return(out.dt)
}

getConseqOnList.parallel <- function(rslist, server,cores) {

  cat("")
  cluster <- parallel::makeCluster(cores)
  registerDoParallel(cluster)

  data <- foreach(i = seq_len(length(rslist)), .combine=rbind,  .packages=c('jsonlite','httr','xml2')) %dopar%
    {
      rs= rslist[i]

      # store the returned scores
      outString <- ""

      conseq <- getVariantConsequence(rs,server)

      if(all(is.element(c('cadd_phred','variant_allele') , names(conseq$consequences[[1]]))))
      {
        conseq.dt <- unique.data.frame(data.table(conseq$consequences[[1]][,'variant_allele'],
                                                  conseq$consequences[[1]][,'cadd_phred']))
        if(length(conseq.dt$V2) > 1)
          conseq.dt <- sapply(conseq.dt, function(x) ifelse(x == "NULL", NA, x))

        conseq.dt <- data.table(matrix(unlist(conseq.dt),nrow = nrow(conseq.dt),byrow = FALSE))
        outString <-  paste(paste(conseq.dt$V1,collapse = ","),
                            paste(conseq.dt$V2,collapse = ","),sep=" = ")
      }

      data.table(rs=rs,cadd=outString)
    }

  stopImplicitCluster()

  cat("done\n",fill = TRUE)

  return(data)
}

find.nearest.gene <- function(ch, pos, gene.set)
{
  data <- gene.set[chr==ch,]
  data[, on.gene := ifelse(chr==ch & pos>start & pos<end, 1, 0)]

  found <- data[on.gene == 1 ,]
  if(nrow(found)>0)
    return(list(paste(found$id,collapse=','),paste(found$name,collapse=',')))

  ##
  data <- gene.set[chr==ch & type=='protein_coding',]

  data[, dist1 := pos-start]
  data[, dist2 := pos-end]

  g1 <- data[,which(dist1<0 & dist2<0)][1]

  v1 <- data[(g1-1)]
  v2 <- data[g1]

  #data <- data[(g1-1):g1]

  # return(list(paste(data$id,collapse = ','),
  #             paste(c(sprintf('%s(dist=%s)',v1$name,v1$dist2),
  #                     sprintf('%s(dist=%s)',v2$name,abs(v2$dist1))),
  #                   collapse = ',')
  #             )
  #        )

  val <- list(sprintf('%s,%s',v1$id,v2$id),
              sprintf('%s(dist=%s), %s(dist=%s)',v1$name,v1$dist2,
                      v2$name,abs(v2$dist1)))
  return(val)


}

find.regulatory<- function(Chr,Pos_37,r.set){
  regulatory <- r.set[chr == Chr & start< Pos_37 & end > Pos_37,]
  if(nrow(regulatory) > 0)
    return(list(paste(regulatory$id,collapse = ','),
                paste(regulatory$type,collapse = ',')))
  else
    return(list('',''))
}

find.band<- function(Chr,Pos_37,c.set){
  band <- c.set[chr == Chr & start< Pos_37 & end > Pos_37, ]
  if(nrow(band) > 0)
    return(band$band)
  else
    return('')
}


check.core.count <- function(cores) {
  avail.cores <- parallel::detectCores()

  if (cores > avail.cores)
    stop('Selected number of cores is not available.',call. = FALSE)
}


createRandString<- function() {
  v = c(sample(0:9, 4, replace = TRUE),
        sample(letters, 1, replace = TRUE))
  return(paste0(v,collapse = ""))
}

checkifOutputIncludesRS <- function(rs,output)
{
  if(is.element(rs,output$Linked_SNP))
    cat(sprintf('== WARNING: variant is in LD with previously checked variant(s): %s',
                paste(output[Linked_SNP == rs, gSNP],collapse = ',')),
        fill = TRUE)
}


generate.report.file <- function(output,html.file,server,db,window_size,r2)
{
  tryCatch({
    options(bitmapType='cairo')
    rmarkdown::render(system.file("rmd", 'variantReport.Rmd', package = "SNPannotator"),
                      output_dir = getwd(),
                      output_file =  html.file,
                      quiet = TRUE)

    return(TRUE)
  },
  error = function(x) {
    cat(x$message,fill = TRUE)
    return(FALSE)
  }
  )
}

select.CADD.scores <- function(alleles,cadd.scores)
{
  alleles <- unlist(strsplit(alleles,'/'))
  cadd.scores <- unlist(strsplit(cadd.scores,' = '))
  cadd.scores1 <- unlist(strsplit(cadd.scores[1],','))
  cadd.scores2 <- unlist(strsplit(cadd.scores[2],','))

  i <- which(is.element(cadd.scores1,alleles))

  return(paste(paste(cadd.scores1[i],collapse = ','),'=',paste(cadd.scores2[i],collapse = ',')))
}
