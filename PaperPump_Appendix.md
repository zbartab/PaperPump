&include "paperfigs/MS_title.md"

---
standalone: true
bibliography: PaperPump.yaml
natbib: true
biblio-title: References
fontsize: 12pt
papersize: a4paper
header-includes:
- \usepackage{double_spaced}
geometry:
- margin=1in
---

The following table shows basic statistics for several attributes of the two real publication networks. Column 'Network' contains the network ID, _n_ is the sample size, _m_ gives the mean, M is the median, 'min' is the minimum, while 'max' is the maximum, 'LQ' and 'UQ' are the lower and upper quartiles, respectively. Terms in italic give the variables: _authors_ is the number of authors a paper has; _papers_ is the number of papers an author has; _weighted papers_ is the number of papers weighted by the reciprocal of their authors; _strengths_ is the strength of a given node (author), i.e. sum of weights the links connected to a given node have; _degrees_ gives the number of neighbours of a node; and _clustering coefficients_ is the local clustering coefficient of a node.

!pagebreak

&include "paperfigs/real_nets_description.txt"

![!FIGURE(A1) The distribution of degrees for two real word collaboration networks (MTMT: www.mtmt.hu, dblp: dblp.org). The lines give the empirical complementary cummulative distribution functions.](paperfigs/real_nets_description.!EXT)

!pagebreak

These tables describe the fit of several theoretical distributions to the degrees of the collaboration networks (Figure 1). The distributions fitted are [@clauset2009]: log-normal, Weibull, exponential, power-law and power-low with exponential cutoff. Column _p_ gives the probability of goodness of fit, _D_~obs~ is the test statistics, _x_~min~ is the minimum value included in the fit, and _n_ is the sample size. Column **Parameters** gives the estimated parameter values for the given distribution.

The MTMT dataset.

&include "paperfigs/MTMT-descr-degrees.txt"

The dblp dataset.

&include "paperfigs/dblp-descr-degrees.txt"

# References
