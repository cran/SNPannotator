.onAttach <-
  function(libname, pkgname) {

    welcomeString = 'initializing SNPannotator package ...'

    packageStartupMessage(welcomeString)

  }


.onLoad <- function(libname, pkgname) {

}


utils::globalVariables(c('r','variation2','variation3','syn_id','rs_id','checkSynonyms','seq_region_name',
                         'variation1','chr','start','Chr',
                         'end','on.gene','type','dist1',
                         'dist2','i','population','feature_type',
                         'external_name','gene_id','id','d_prime',
                         'name','size','description','r2','Gene',
                         'Pos_37','Cytoband','gSNP','Linked_SNP','number','Alleles'))
