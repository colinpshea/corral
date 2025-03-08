% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/determineAlleleMatches.R
\name{determineAllAlleleMatches}
\alias{determineAllAlleleMatches}
\title{Assess matching of alleles, comparing to selves and other colonies.}
\usage{
determineAllAlleleMatches(dataset)
}
\arguments{
\item{dataset}{A data frame with \code{Coral_ID} as the first column and single-letter \code{IUPAC} allele data for \code{nrow(dataset)} colonies at each of \code{2:ncol(dataset)} loci.}
}
\value{
This function returns a data frame with six columns: \code{coral1} (name of coral 1), \code{coral2} (name of coral 2), \code{allele1} (allele for coral 1 at that locus), \code{allele2} (allele for coral 2 at that locus), \code{locus} (the locus name), and match (are the \code{IUPAC} alleles the same: \code{TRUE} or \code{FALSE}). All possible pairwise comparisons with selves and others are included here. The resulting file is used by \code{groupByGenets()}, which assigns colonies to genets.
}
\description{
For each pairwise combination of individuals, this function classifies allele pairs at each locus as either matching or not. This function applies \code{classifyAllelePairs()}, one column or locus at a time, by looping sequentially through all of \code{2:ncol(dataset)} loci. This function includes allele matches at each locus for colonies compared to themselves and compared to other colonies and is specifically used for genet assignment, which requires both types of comparisons.
}
