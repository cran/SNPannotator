% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getVariantLDs.R
\name{LDlist}
\alias{LDlist}
\title{Find variants in high LD with the lead SNP.}
\usage{
LDlist(rslist, server, db, window_size, r2)
}
\arguments{
\item{rslist}{A vector of rs numbers.}

\item{server}{Name of the server. "https://rest.ensembl.org" can be used for GRCh38
and "https://grch37.rest.ensembl.org" for GRCh37.}

\item{db}{The population database for calculating LD scores. This can be found using \code{listDatabases} function.}

\item{window_size}{Number of base pairs around the variant for checking LD scores (max = 500kb)}

\item{r2}{The LD threshold for selecting variants around the target SNP.}
}
\value{
a data table with variant information.
}
\description{
This function returns a list of variables that are in high LD with the lead variant.
}
