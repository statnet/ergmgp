\docType{package}
\name{ergmgp-package}
\alias{ergmgp-package}
\alias{ergmgp}
\alias{EGP_init}

\title{Tools for Modeling ERGM Generating Processes}

\description{
Tools for simulation and analysis of continuous time graph processes with equilibria that can be described in exponential family random graph (ERGM) form.
}

\details{
A random graph \eqn{G} on support \eqn{\mathcal{G}}{Supp(G)} is said to be expressed in exponential family random graph (ERGM) form when its probability mass function (pmf) is written as

\deqn{\Pr(G=g|\theta,X) = \frac{\exp\left(\theta^T w(g,X)\right) h(g)}{\sum_{g' \in \mathcal{G}} \exp\left(\theta^T w(g',X)\right) h(g')}}{%
      Pr(G=g|theta,X) = exp[theta^T w(g,X)] h(g) / 
                     [ sum_{g' in Supp(G)} exp[theta^T w(g',X)] h(g')]
}

where \eqn{\theta}{theta} is a parameter vector, \eqn{w} is a vector of sufficient statistics, \eqn{X} is a covariate set, and \eqn{h} is a reference measure.  This form is quite general; in particular, any pmf on finite support can be written in ERGM form (albeit not always elegantly), making it a natural language for expressing graph distributions.

Now, consider a continuous time process whose state space is \eqn{\mathcal{G}}{Supp(G)}.  A process of this type having an equilibrium distribution that can be written in (known) ERGM form is said to be an \emph{ERGM generating process} or EGP.  Although there are many types of EGPs, perhaps the most natural are continuous time Markov chains (CTMCs) whose transitions involve the addition or removal of individual edges from the current graph state.  The transition rates of such CTMCs have the topology of the Hamming adjacency on \eqn{\mathcal{G}}{Supp(G)}; this is only sensible when considering graphs on a fixed vertex set, which is the typical use case.  We can think of this class of EGPs as continuous time analogs of the Markov chains used in ERGM simulation (see \code{\link[ergm]{ergm}}), with equilibrium obtained in the limit of time rather than simulation steps.  EGPs are potentially useful as dynamic interpretations of empirically obtained ERGMs, or as \emph{a priori} models in their own right.  Since many Markovian EGPs are identified by their equilibrium ERGM together with a pacing constant, they are also natural choices when dynamics must be inferred from limited data (e.g., a cross-sectional network observation together with pacing or duration information).

The \code{ergmgp} package supports a number of different EGPs, all of which are currently Markovian with support on graphs or digraphs of fixed order.  The following EGPs are currently supported; all definitions and notation follow Butts (2023).  Define  \eqn{q(g) = \theta^T w(g,X) + \ln h(g)}{q(g) = \theta^T w(g,X) + ln h(g)} to be the \emph{graph potential}; in some cases, separate potentials (\eqn{q_f} and \eqn{q_d}) may be employed for the formation and dissolution of edges.  For brevity, we define the normalizing factor of the equilibrium ERGM distribution by \eqn{Z = \sum_{g' \in \mathcal{G}} q(g')}{Z = sum_{g' in Supp(G)} q(g')}, let \eqn{w_e} be the edge count statistic, and let \eqn{\mathcal{H}(g)}{H(g)} be the Hamming neighborhood of \eqn{g} (i.e., the set of graphs reachable by single edge changes, or \dQuote{toggles}).  All transition rates from graph \eqn{a} to \eqn{b} (denoted \eqn{a \to b}{a -> b}) zero for \eqn{b \not\in \mathcal{H}(a)}{b not in H(a)}.  Processes not otherwise noted were introduced in Butts (2023).

\describe{
  \item{Longitudinal ERGM (\code{LERGM})}{
    \describe{
      \item{Description:}{ Introduced by Koskinen and Snijders (2007), this process was originally conceived of as a continuum analog to the Gibbs sampler, with transition rates that are increasing with differences in graph potential.  Grazioli et al. (2019) subsequently showed that it can also be derived from a physical model with locally Arrhenius-like kinetics.  This process has a maximum change rate (but no minimum), and may thus be plausible in settings for which changes can only be made when (exogenously determined) opportunities arise.
      }
      \item{Event Rate (\eqn{a \to b}{a -> b}):}{
        \eqn{[1+\exp[q(a)-q(b)]]^{-1}}{1/[1 + exp[q(a) - q(b)]]}
      }
      \item{Equilibrium:}{
        \eqn{\exp[q(a)]/Z}{exp[q(a)]/Z}
      }
    }
  }
  \item{Competing Rate Stochastic Actor-Oriented Model (\code{CRSAOM})}{
    \describe{
      \item{Description:}{ Introduced as an EGP by Snijders (2001), this model was originally proposed as a behavioral process, where vertices represent actors controlling their outgoing edges, the rate at which actors make tie changes is a function of the attractiveness of the networks reachable by making such changes, and (given opportunity to act) edge changes are chosen by a multinomial logit with utility function \eqn{q}.  Dynamics in this model are distinctive in being driven solely by the attractiveness of the target state, which can sometimes lead to rapid state switching when multiple high-potential states are Hamming-adjacent.
      }
      \item{Event Rate (\eqn{a \to b}{a -> b}):}{
        \eqn{\exp[q(b)]}{exp[q(b)]}
      }
      \item{Equilibrium:}{
        \eqn{\exp[q(a)]/Z}{exp[q(a)]/Z}
      }
    }
  }
  \item{Change Inhibition (\code{CI})}{
    \describe{
      \item{Description:}{  In the same sense that the LERGM is analogous to a continuum Gibbs sampler, this process is loosely analogous to continuum Metropolis algorithm.  Downhill transitions with respect to the graph potential occur with a rate that is decreasing with the potential difference; uphill transitions, however, occur at a fixed rate (irrespective of the potential difference).  The process thus works by selectively inhibiting downhill moves, rather than by preferentially moving to graphs of highest local potential.  
      }
      \item{Event Rate (\eqn{a \to b}{a -> b}):}{
        \eqn{\min(1,\exp[q(b)-q(a)])}{min(1, exp[q(b) - q(a)])}
      }
      \item{Equilibrium:}{
        \eqn{\exp[q(a)]/Z}{exp[q(a)]/Z}
      }
    }
  }
  \item{Differential Statibility (\code{DS})}{
    \describe{
      \item{Description:}{ Analogous to a \dQuote{win-stay, lose-shift} process, transition targets in this EGP are chosen uniformly at random, with structure arising entirely from transition \emph{times}.  The time to exit a state \eqn{g} is proportional to \eqn{\exp[q(g)]}{exp[q(g)]}, making high-potential states exponentially more persistent than low-potential states.  Note that this process is in a sense the inverse of the \code{CRSAOM}, being dependent only on the potential of the source state (while the \code{CRSAOM} depends only on the potential of the target state).  Since the transitions themselves from a random walk, it should be noted that this process can generate very large numbers of transition events involving low-potential states that, while taking little phenomenological time, nevertheless are expensive to compute.
      }
      \item{Event Rate (\eqn{a \to b}{a -> b}):}{
        \eqn{|\mathcal{H}(a)|^{-1} \exp[-q(a)]}{exp[-q(a)]/|H(a)|}
      }
      \item{Equilibrium:}{
        \eqn{\exp[q(a)]/Z}{exp[q(a)]/Z}
      }
    }
  }
  \item{Constant Dissolution Continuum STERGM (\code{CDCSTERGM})}{
    \describe{
      \item{Description:}{ This process is a continuum version of the discrete time constant dissolution separable temporal ERGM (STERGM) introduced by Carnegie et al. (2015); here, edges are lost randomly at a fixed rate, with a formation potential \eqn{q_f} that governs edge addition.  This is a special case of the continuum STERGMs (below), and is particularly easy to identify from limited information.  
      }
      \item{Event Rate (\eqn{a \to b}{a -> b}):}{
        If \eqn{b} is formed by adding an edge to \eqn{a}, then \eqn{\exp[q_f(b)-q_f(a)]}{exp[q_f(b) - q_f(a)]}; otherwise \eqn{\exp[\theta_d]}{exp[theta_d]}
      }
      \item{Equilibrium:}{
        \eqn{\exp[q_f(a) - \theta_d w_e(a)]/Z}{exp[q_f(a) - theta_d w_e(a)]/Z}
      }
    }
  }
  \item{Constant Formation Continuum STERGM (\code{CFCSTERGM})}{
    \describe{
      \item{Description:}{ This process is analogous to the \code{CDCSTERGM}, except that in this case edge \emph{formation} occurs randomly at a fixed rate, with a dissolution potential \eqn{q_d} governing edge loss.  It is a simple model for settings in which edges arise from essentially idiosyncratic events, with the resulting network structure subsequently stabilizing or destabilizing particular edges.
      }
      \item{Event Rate (\eqn{a \to b}{a -> b}):}{
        If \eqn{b} is formed by adding an edge to \eqn{a}, then \eqn{\exp[\theta_f]}{exp[theta_f]}; otherwise \eqn{\exp[q_d(b)-q_d(a)]}{exp[q_d(b) - q_d(a)]}
      }
      \item{Equilibrium:}{
        \eqn{\exp[q_d(a) + \theta_f w_e(a)]/Z}{exp[q_d(a) + theta_f w_e(a)]/Z}
      }
    }
  }
  \item{Continuum STERGM (\code{CSTERGM})}{
    \describe{
      \item{Description:}{ This process represents a continuum limit of the discrete time separable temporal ERGMs (STERGMs) introduced by Krivitsky and Handcock (2014).  Edge formation is here governed by one potential (\eqn{q_f}), while dissolution is governed by another (\eqn{q_d}), allowing these processes to be governed by different effects.  The resulting equilibrium pmf is based on the sum of both potentials.
      }
      \item{Event Rate (\eqn{a \to b}{a -> b}):}{
        If \eqn{b} is formed by adding an edge to \eqn{a}, then \eqn{\exp[q_f(b)-q_f(a)]}{exp[q_f(b) - q_f(a)]}; otherwise \eqn{\exp[q_d(b)-q_d(a)]}{exp[q_d(b) - q_d(a)]}
      }
      \item{Equilibrium:}{
        \eqn{\exp[q_d(a) + q_f(a)]/Z}{exp[q_d(a) + q_f(a)]/Z}
      }
    }
  }
  \item{Continuum TERGM (\code{CTERGM})}{
    \describe{
      \item{Description:}{ This process is a continuum limit of the discrete time temporal ERGMs (TERGMs) introduced in Robins and Pattison (2001).  The transition rates for this class are particularly natural, with the log rates being equal to the potential differences between states.  Note that the potential of the equilibrium ERGM is scaled by a factor of 2 from the transition potential (as can be obtained from the \code{CSTERGM} by letting \eqn{q_f=q_d}); intuitively, this arises because states of higher potential are both more stable (lower exit rates) and more attractive (higher entrance rates) than states of lower potential.
      }
      \item{Event Rate (\eqn{a \to b}{a -> b}):}{
        \eqn{\exp[q(b)-q(a)]}{exp[q(b) - q(a)]}
      }
      \item{Equilibrium:}{
        \eqn{\exp[2q(a)]/Z}{exp[2 q(a)]/Z}
      }
    }
  }
}

Further details on each process can be found in Butts (2023).  All of the above transition rates are defined up to an arbitrary pacing constant (which is generally specified separately, and taken to be 1 in package tools if not otherwise indicated).  Note that the LERGM and Change Inhibition processes have unit-maximum transition rates, and thus the pacing constant sets the maximum rate of change.

Information on functions for simulation or analysis of EGPs is provided in their respective manual pages.  Information on ERGMs and their specification can be found within the \code{\link[ergm]{ergm}} page in the \code{ergm} library.
}

\seealso{
\code{\link{simEGP}}, \code{\link{EGPHazard}}, \code{\link{EGPRateEst}}, \code{\link[ergm]{ergm}}, \code{\link{durations}}
}

\references{
Butts, Carter T.  (2023).  \dQuote{Continuous Time Graph Processes with Known ERGM Equilibria: Contextual Review, Extensions, and Synthesis.} \emph{Journal of Mathematical Sociology}.  \doi{10.1080/0022250X.2023.2180001}

Carnegie, Nicole B.; Krivitsky, Pavel N.; Hunter, David R.; and Goodreau, Steven M. (2015). \dQuote{An Approximation Method for Improving Dynamic Network Model Fitting.}  \emph{Journal of Computational and Graphical Statistics}, 24(2):502-519.

Grazioli, Gianmarc; Yu, Yue; Unhelkar, Megha H.; Martin, Rachel W.; and Butts, Carter T. (2019). \dQuote{Network-based Classification and Modeling of Amyloid Fibrils.}  \emph{Journal of Physical Chemistry, B}, 123(26):5452-5462.

Koskinen, Johan H. and Snijders, Tom A. (2007). \dQuote{Bayesian Inference for Dynamic Social Network Data.}  \emph{Journal of Statistical Planning and Inference}, 137(12):393--3938. 5th St. Petersburg Workshop on Simulation, Part II.

Krivitsky, Pavel N. and Handcock, Mark S. (2014). \dQuote{A Separable Model for Dynamic Networks.} \emph{Journal of the Royal Statistical Society, Series B,} 76(1):29-46.

Robins, Garry L. and Pattison, Philippa E. (2001). \dQuote{Random Graph Models for Temporal Processes in Social Networks.} \emph{Journal of Mathematical Sociology}, 25:5-41.

Snijders, Tom A. B. (2001). \dQuote{The Statistical Evaluation of Social Network Dynamics.}  \emph{Sociological Methodology,} 31:361-395.
}

\author{
Carter T. Butts \email{buttsc@uci.edu}
}
\keyword{package}
