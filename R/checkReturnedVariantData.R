checkReturnedVariantData <- function(varInfo)
{
  #cat('\nData verification for target SNP ... ')

  #correct <- TRUE

  # ldVars <- unlist(varInfo$LDList$variation2)

  l1 <- as.data.table(varInfo$LDList)

  l2 <- data.frame(matrix(NA,length(varInfo$LDlistFull),2))
  names(l2) <- c('rs_id','syn_id')

  for(i in seq_len(length(varInfo$LDlistFull)))
  {
    l2[i,] = c(varInfo$LDlistFull[[i]]$name,
               paste(varInfo$LDlistFull[[i]]$synonym,collapse = ','))
  }

  setDT(l2)


  l1[,variation2 :=ifelse(variation2 %in% l2$rs_id,
                          variation2,
                          l2[grepl(pattern = paste0(variation2,'\\b'), x=syn_id),rs_id]),
     by = variation2]


  return(l1)

  # if(length(ldVars) != length(returnedVars))
  # {
  #
  #   cat(sprintf("RSID count mismatch! in LD list = %s , in information list = %s",
  #               length(ldVars),
  #               length(returnedVars)
  #   ),fill = TRUE)
  #   correct <- FALSE
  # }




  # wrongs <- which(ldVars != returnedVars)
  #
  # if(length(wrongs) > 0)
  # {
  #   cat(sprintf("RSID mismatch! in LD list = %s , in information list = %s",
  #               paste(ldVars[wrongs],collapse = "/"),
  #               paste(returnedVars[wrongs],collapse = "/")
  #   ),fill = TRUE)
  #
  #   correct <- FALSE
  # }

  # if(correct)
  #   cat('done\n\n')
  # else
  #   cat('problem found\n\n')

}
