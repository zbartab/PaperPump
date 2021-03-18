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

&define TCG tightly connected group

# Receiver operating characteristic (ROC) of cut-off rule for identifying cartels

To investigate the robustness of the cut-off rule to identify cartels we generated many random publication networks with cartels added and perform ROC analyses on them. The results indicate that the performance of this simple rule is rather independent of two important parameters of the network generation (Figure A1 and A2).

![!FIGURE(A1) The effect of $\alpha$ parameter of Beta distirbution on the results of ROC analyses of the cut-off algorithm for identifying cartel members. The insert illustrates the probability distributions of the probabilities of sharing a paper's authorship for several $\alpha$s. Dots represent means, while error bars SDs for 25 independent simmulations.](paperfigs/cartel_footprint-alpha.!EXT)

![!FIGURE(A2) The effect of parameter _D_ on the results of ROC analyses of the cut-off algorithm for identifying cartel members. Parameter _D_ symbolise the density of publication networks. Dots represent means, while error bars SDs for 25 independent simmulations.](paperfigs/cartel_footprint-p.!EXT)

!pagebreak

# Description of real publication networks

The following table shows basic statistics for several attributes of the two real publication networks. Column 'Network' contains the network ID, _n_ is the sample size, _m_ gives the mean, M is the median, 'min' is the minimum, while 'max' is the maximum, 'LQ' and 'UQ' are the lower and upper quartiles, respectively. Terms in italic give the variables: _authors_ is the number of authors a paper has; _papers_ is the number of papers an author has; _weighted papers_ is the number of papers weighted by the reciprocal of their authors; _strengths_ is the strength of a given node (author), i.e. sum of weights the links connected to a given node have; _degrees_ gives the number of neighbours of a node; and _clustering coefficients_ is the local clustering coefficient of a node.

!pagebreak

&include "paperfigs/real_nets_description.txt"

![!FIGURE(A3) The distribution of degrees for two real word collaboration networks (MTMT: www.mtmt.hu, dblp: dblp.org). The lines give the empirical complementary cummulative distribution functions.](paperfigs/real_nets_description.!EXT)

!pagebreak

These tables describe the fit of several theoretical distributions to the degrees of the collaboration networks (Figure 1). The distributions fitted are [@clauset2009]: log-normal, Weibull, exponential, power-law and power-low with exponential cutoff. Column _p_ gives the probability of goodness of fit, _D_~obs~ is the test statistics, _x_~min~ is the minimum value included in the fit, and _n_ is the sample size. Column **Parameters** gives the estimated parameter values for the given distribution.

The MTMT dataset.

&include "paperfigs/MTMT-descr-degrees.txt"

The dblp dataset.

&include "paperfigs/dblp-descr-degrees.txt"

# Tightly connected groups in different subject areas

![!FIGURE(A4) Signs of tightly connected groups in real networks by subject areas. The figure shows the empirical complementary cummulative distribution function (eCCDF) of link weights in the two real publication network (dblp and MTMT) and in subnetworks of MTMT deliniated by Elsevier subject areas.](paperfigs/subject_nets_MTMT-descr.!EXT)

The following table shows the size of the subnetworks in MTMT and the proportion of authors participating in tightly connected groups.

&include "paperfigs/subject_nets_MTMT-descr.txt"

# Level of paper sharing

The relation between average group productivity (i.e. the number of distinct papers produced by a group divided by the group size) and the average number of papers authors authored in the group suggests that authors share a constant proportion of their papers with others in the group. This is supported by the following calculations.

![!FIGURE(A5) The relationship between the average group productivity and average number of papers for random groups and !TCGs for the MTMT (left) and the dblp (right) dataset. The lines are fitted by linear models. Note, the lines for the random groups have slopes of 1, while the slopes of lines for the !TCGs are less than 2.](paperfigs/paper_sharing-relations.!EXT)

Consider a group of _a_ authors, where author _i_ produces _p_~i~ papers alone. Author can share _s_~i~ of her papers with others in the group. Under these conditions the group's productivity is given as

$$
\sum_{i=1}^{a} p_i,
$$

while the average group productivity, _G_~P~, is

$$
G_P = \frac{1}{a}\sum_{i=1}^{a} p_i.
$$

The number of papers produced by author _i_ is

$$
p_i + \sum_{i \ne j} s_j = p_i + \sum_{j=1}^{a} s_j - s_i.
$$

From this the average number of papers, _N_~P~, authored by an author in the group follows:

$$
N_P = \frac{1}{a} \sum_{i=1}^{a}{\left(p_i +\sum_{j=1}^{a}s_j - s_i\right)}.
$$

After some simplification we got:

$$
N_P = \frac{1}{a} \left[\sum_{i=1}^{a} p_i + (a-1)\sum_{i=1}^{a} s_i\right].
$$

We consider three cases. First, _s_~i~ = _p_~i~:

$$
N_P = \frac{1}{a}\left[\sum_{i=1}^{a} p_i + (a-1)\sum_{i=1}^{a} p_i\right] = \sum_{i=1}^{a} p_i.
$$

It follows that _N_~P~ = _a_ _G_~P~. In other words, the average number of papers is a linear function of the average group productivity with group size as slope. As $a \ge 2$, the slope should be larger than 2.

For the second case we assume that each individual in the group shares the same number of papers with others, i.e. $s_i = s_j = s,\, \forall i \ne j$. Then

$$
N_P = \frac{1}{a} \left[\sum_{i=1}^{a} p_i + (a-1)a s\right] = (a-1)s + \frac{1}{a}\sum_{i=1}^{a} p_i = (a-1)s + G_P.
$$

Here _N_~P~ is also a linear function of _G_~P~ with slope equals 1.

For the third case we assume that each author shares a constant portion, _c_, of her papers with others, i.e. _s_~i~ = _c_ _p_~i~, _c_ < 1. Under this condition,

$$
N_P = \frac{1}{a} \left[\sum_{i=1}^{a} p_i + (a-1)\sum_{i=1}^{a} c p_i\right].
$$

After some simplification we get:

$$
N_P = \left[1 + c(a-1)\right] G_P.
$$

Let $\rho = 1 + c(a-1)$. Then $\rho$ can be estimated as $\rho = N_P / G_P$. From this the level of sharing, _c_, can be obtained as

$$
c = \frac{\rho -1}{a-1}.
$$

# References
