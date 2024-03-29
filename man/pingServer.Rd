% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/codes.R
\name{pingServer}
\alias{pingServer}
\title{Checks if the service is alive}
\usage{
pingServer(server)
}
\arguments{
\item{server}{name of the server. "https://rest.ensembl.org" can be used for GRCh38
and "https://grch37.rest.ensembl.org" for GRCh37.}
}
\value{
a message is displayed to the user
}
\description{
This function test whether the Ensembl server is accessible or not
}
