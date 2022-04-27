convert.CADD <- function(cadd)
{
  if(grepl(x = cadd,pattern = '='))
    cadd =gsub(pattern = '(.+=\\s)(.+)',replacement = '\\2',x = cadd)

  cadd = as.numeric(cadd)

  if(is.na(cadd) | grepl(x = cadd, pattern = ','))
    return('')
  else
    return(sprintf('top %#.1f%%',(10^(-cadd/10))*100))

}


select.CADD.scores <- function(alleles,cadd.scores)
{
  alleles <- unlist(strsplit(alleles,'/'))
  cadd.scores <- unlist(strsplit(cadd.scores,' = '))
  cadd.scores1 <- unlist(strsplit(cadd.scores[1],','))
  cadd.scores2 <- unlist(strsplit(cadd.scores[2],','))

  i <- which(is.element(cadd.scores1,alleles))[1]

  return(paste(paste(cadd.scores1[i],collapse = ','),'=',paste(cadd.scores2[i],collapse = ',')))
}
