% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/codes.R
\name{listDatabases}
\alias{listDatabases}
\title{List population from human database (1000 Genomes project)}
\usage{
listDatabases(server)
}
\arguments{
\item{server}{name of the server. "https://rest.ensembl.org" can be used for GRCh38
and "https://grch37.rest.ensembl.org" for GRCh37.}
}
\value{
A data table is returned which includes the name, description and size of the available populations
in 1000 Genomes project database.
}
\description{
This function list the name, description and size of the available populations
in 1000 Genomes project database. This database will be used for returning variables in high LD
with the target SNP.
}
