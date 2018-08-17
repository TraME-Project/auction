\name{auction}

\alias{auction}

\title{Auction Algorithms for Transportation Problems}

\description{
    Several variations on the auction algorithms for real-valued transportation problems.
}

\usage{
auction(Phi,demand,supply,algorithm=1,max_prob=TRUE)
}

\arguments{
    \item{Phi}{
    Costs}
	\item{demand}{
	Demand vector.}
	\item{supply}{
	Supply vector.}
    \item{algorithm=1}{
	1 for general auction; 2 for assignment auction; 3 for auction with similar objects; 4 for similar objects and persons.}
    \item{max_prob=TRUE}{
    Maximization problem or minimization problem.
    }
    \item{eps_init=-1.0}{
    Initial step size.
    }
    \item{eps_min=-1.0}{
    Minimum step size.
    }
    \item{run_bellman_ford=FALSE}{
    Run Bellman-Ford algorithm in the event of a duality gap.
    }
    \item{bf_upper_prices=TRUE}{
    Return the upper (TRUE) or lower bounds (FALSE) on prices.
    }
    \item{verbose=TRUE}{
    Verbose printing.
    }

}

\details{
    See the references below.
}

\value{
The function returns a list of objects...
}

\author{Keith O'Hara}

\references{
    Walsh, J.D. and Luca Dieci, \emph{General auction method for real-valued optimal transport,} \emph{Mimeo}, 2010.
}

\examples{
\dontrun{
n <- 3
m <- 2

cost <- matrix(c(3,6,7,2,8,2), nrow = n)

p  <-  rep(1,n)
q  <-  rep(1,m)

q[1] <- 2

auction(cost,q,p,1,FALSE)
}
}

\keyword{auction}
\concept{auction}