getVariantData <- function(varNum, gVar, linkedVar, database , LD)
{

  ## alleles and frequency
  pops <- linkedVar$populations
  setDT(pops)
  # pops <- pops[grep(x= population,pattern = population),]
  pops <- pops[population == database,]
  alleles <- paste(pops$allele,collapse = "/")
  freqs <- paste(pops$frequency,collapse = "/")

  ancestral_Allele <- ifelse(inherits(linkedVar$mappings$ancestral_allele, 'list'),
                              linkedVar$mappings$ancestral_allele[[1]],
                              "")
  ## create GWASCATALOG phenotype string
  phenotypes <- linkedVar$phenotypes

  if(!is.null(nrow(phenotypes)))
    phenotypes <- paste(unique(phenotypes$trait),collapse = ", ") else
      phenotypes <- ""


  ## create Gene table
  geneTab <- data.table()

  # replace if data exists
  if(!is.null(linkedVar$geneInfo))
  {

    geneTab <- linkedVar$geneInfo
  # some times description is data.frame instead of list
   if(inherits(geneTab$description, 'list'))
     geneTab$description <- NULL

   setDT(geneTab)
  }

  geneNames <- ""
  geneIDs <- ""
  #geneDesc <- ""
  band <- ""

  if(nrow(geneTab) > 0 && geneTab[feature_type == 'gene', .N] > 0)
  {
    #geneNames <- paste(geneTab[feature_type == 'gene', external_name], collapse = "/")
    #geneIDs <- paste(geneTab[feature_type == 'gene', gene_id], collapse = "/")

    geneNames <- paste(geneTab[feature_type == 'gene',]$external_name, collapse = "/")
    geneIDs <- paste(geneTab[feature_type == 'gene', ]$gene_id, collapse = "/")
    #geneDesc <- paste(geneTab[feature_type == 'gene', description], collapse = "/")
  }

  if(nrow(geneTab) > 0 && geneTab[feature_type == 'band', .N] > 0)
  {
    band <- paste0(geneTab[feature_type == "band",]$seq_region_name,
                   geneTab[feature_type == "band",]$id)
  }

  ## create table
  tab <- data.table(`#gSNP` = varNum,
                    gSNP = gVar$name,
                    Linked_SNP = linkedVar$name,
                    Chr = linkedVar$mappings$seq_region_name[[1]],
                    Pos_37 = linkedVar$mappings$start[[1]],
                    # Alleles = linkedVar$mappings$allele_string[[1]],
                    # `Ancestral allele` = ancestral_Allele,
                    # `Minor allele`  = linkedVar$minor_allele,
                    # MAF = linkedVar$MAF,
                    Alleles = alleles,
                    `Allele frequencies` = freqs,
                    LD = LD ,
                    Cytoband = band,
                    Type = linkedVar$most_severe_consequence,
                    Gene = geneNames,
                    GeneId = geneIDs,
                    #  GeneDescription = geneDesc,
                    Associations = phenotypes
  )

  return(tab)
}
