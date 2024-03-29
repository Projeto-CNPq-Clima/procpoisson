% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MainWeibullNewSA.R
\name{MainWeibullNewSA}
\alias{MainWeibullNewSA}
\title{Algoritmo MCMC Poisson Process}
\usage{
MainWeibullNewSA(data, sites, iter = 1000, bar = 900)
}
\arguments{
\item{data}{Matriz de dados simulados de um processo de Poisson}

\item{sites}{Malha de pontos}

\item{iter}{número inteiro positivo, interações do algoritmo MCMC}

\item{bar}{número inteiro positivo, período de burn-in do algoritmo MCMC}
}
\description{
Esta é a função para realizar o algoritmo Monte Carlo via Cadeias de Markov para os parâmetros de um processo de Poisson.
}
