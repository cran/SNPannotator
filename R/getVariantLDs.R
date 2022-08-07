#' Find variants in high LD with the lead SNP.
#'
#' This function returns a list of variables that are in high LD with the lead variant.
#'
#' @param rslist A vector of rs numbers.
#' @param server Name of the server. "https://rest.ensembl.org" can be used for GRCh38
#' and "https://grch37.rest.ensembl.org" for GRCh37.
#' @param db The population database for calculating LD scores. This can be found using `listDatabases` function.
#' @param window_size Number of base pairs around the variant for checking LD scores (max = 500kb)
#' @param r2 The LD threshold for selecting variants around the target SNP.
#' @return a data table with variant information.
#'
LDlist <- function(rslist,server,db, window_size, r2)
{

  start.time <-  proc.time()

  on.exit({
    if(!is.null(start.time))
      cat(sprintf("\nRun time: %s",timetaken(start.time)))# END LOG
  })


  pb <- progress::progress_bar$new(format = "[:bar] :current/:total (:percent)", total = length(rslist))
  pb$tick(0)

  output <- data.frame()
  cat("Fetching variants in high LD",fill = TRUE)

  for(i in 1:length(rslist))
  {
    rs= rslist[i]

    # get variants in high LD
    ldlist <- getVariantLDs(rs,server,db, window_size, r2)


     if(!is.null(ldlist))
    {
       # add variant number in the first column
       ldlist$number = i

       output <- rbind(output,ldlist)
    }
    pb$tick(1)
  }

  if(nrow(output)>0)
  {
    setDT(output)
    output <- output[order(-r2),
                     list(variation2,r2),
                     keyby=list(number,variation1)]
    names(output) <- c('number','variant1','variant2','R2')
    return(output)
  }
  else
    return(NULL)
}

getVariantLDs <- function(rsID,server,db, window_size, r2)
{

  ext <- sprintf("/ld/human/%s/%s?window_size=%s;r2=%s",
                 rsID,
                 db,
                 window_size,
                 r2)
  resp <- fetch(ext,server)
  LDlist <- setDT(resp)

  if(!is.null(LDlist) && nrow(LDlist) > 0)

  {
    LDlist <- LDlist[, lapply(.SD, as.character)]
    non.rs <- LDlist[!startsWith(variation2 ,prefix = "rs"), .N]

    # cat(sprintf("\n== Found %s variants in high LD",
    #             nrow(LDlist)),
    #     fill=TRUE)

    # remove variants that don't start with rs , e.g. esv569117
    if(LDlist[!startsWith(variation2 ,prefix = "rs"), .N] > 0)
    {
      cat(sprintf("Inappropriate IDs: %s",
                  paste(LDlist[!startsWith(variation2 ,prefix = "rs"),variation2],
                        collapse ="|")))

      LDlist <- LDlist[startsWith(variation2 ,prefix = "rs"),]


    }

    LDlist <- LDlist[ ,r2:= as.numeric(r2)]
    LDlist <- LDlist[ ,d_prime:= as.numeric(d_prime)]
    LDlist[order(-r2)]
    return(LDlist)
  }
  else
  {
    #cat("No variants in high LD found.",fill = TRUE)
    return(NULL)
  }

}
