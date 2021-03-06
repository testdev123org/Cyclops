% Generated by roxygen2 (4.0.1): do not edit by hand
\name{isInitialized}
\alias{isInitialized}
\title{isInitialized}
\usage{
isInitialized(object)
}
\arguments{
\item{object}{OHDSI data object to test}
}
\description{
\code{isInitialized} determines if an OHDSI data object is properly
initialized and remains in memory.  OHSDI data objects do not
serialized/deserialize their back-end memory across R sessions.
}

