returnVariantDatatable <- function(varNum , varInfo, database )
{

  # create an empty data table where the first row is the target SNP
  tab <-  getVariantData(varNum,varInfo,varInfo,database, LD = 1)

  # fill the above data table if any variant found in high LD
  if(!is.null(varInfo$LDList) & !is.null(ncol(varInfo$LDList)))
  {
    LDTable <- setDT(varInfo$LDList)

    for(i in seq_len(length(varInfo$LDlistFull)))
    {
      LD = LDTable[variation2 == varInfo$LDlistFull[[i]]$name, r2]
      tab <- rbind(tab, getVariantData(varNum,varInfo,varInfo$LDlistFull[[i]],database,LD))
    }
  }

  tab <- tab[order(-LD),]

  return(tab)
}
