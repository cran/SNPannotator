inSilicoSeqPipeline <- function(rsID,
                                server,
                                database,
                                window_size,
                                r2,
                                i,
                                variantCount,
                                addLinkedVars,
                                get.gene.info.from.API)
{

  rsID <- tolower(rsID)
  if(!startsWith(rsID,"rs"))
    stop("rsID is not correct.",call. = FALSE)

  #cat("=============================================", fill = TRUE)
  cat(sprintf("Fetching data for target SNP (#%s/%s): %s ",i,variantCount,rsID), fill = TRUE)
  # cat("=============================================", fill = TRUE)

  varInfo <- getVariantInfo(rsID,server,NULL,get.gene.info.from.API)

  if(is.null(varInfo))
    return(NULL)

  if(!addLinkedVars)
    return(varInfo)
  else
    varInfo$LDList <- getVariantLDs(rsID,server,database, window_size, r2)

  if(is.null(varInfo$LDList))
  {
    varInfo$LDcount <- 0
    varInfo$LDlistFull <- NULL
  }else
  {
    varInfo$LDcount <- nrow(varInfo$LDList)
    varInfo$LDlistFull <- list()

    pb <- progress::progress_bar$new(format = "[:bar] :current/:total (:percent)", total = varInfo$LDcount)
    pb$tick(0)

    for(j in seq_len(varInfo$LDcount))
    {
      #cat(sprintf("\n%s- %s (%s/%s) - rsID: %s",i, rsID ,j,varInfo$LDcount,varInfo$LDList[j,variation2]), fill=TRUE)
      varInfo$LDlistFull[[j]] <- getVariantInfo(varInfo$LDList[j,variation2] , server, pb,get.gene.info.from.API)
    }
  }

  #varInfo$LDList <- correctSynonymIDs(varInfo)

  return(varInfo)
}
