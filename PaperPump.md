% Why correcting for coauthorships fails?
% Z Barta
% 2019-10-13

---
standalone: true
toc: true
bibliography: /home/apa/notes/articles.bib
natbib: true
biblio-files: /home/apa/notes/articles.bib
biblio-style: /home/apa/lib/texinputs/AmNat
biblio-title: References
fontsize: 12pt
papersize: a4paper
header-includes:
- \usepackage{double_spaced}
geometry:
- margin=1in
---

<!--
# TODOs

- write introduction
-->

&define TODO(tobedone)
&define pubcart publication cartel
&define pub publication
&define pubnet publication network
&define colnet collaboration network
&define pubmat publication matrix
&define pubmats publication matrices
&define colmat collaboration matrix
&define colmats collaboration matrices
&define colgraph collaboration graph
&define susgroup suspicious group

# Introduction

Under current conditions a scientist' reputation is proportional to the number of papers he or she (co)authored. This puts a great pressure on scientists to increase their number of publications, a process nicely summarised by the term "publish or perish". An easy way to increase the number of publications is to form "!pubcarts" where authors cooperate to form a group and reciprocally invite each others to be coauthor on their papers.

# The model

We compare the publication performance of two authors, author _A_ and author _B_. Authors work in separate groups (group _A_ and group _B_, respectively) each of which contains _G_~i~ (_i_ = _A_ or _B_) people (including the focal author). Each author in group _A_ produces _p_~A~ papers in a year by collaborating with _c_~A~ authors from outside of the group. Similarly, authors in group _B_ produce _p_~B~ papers by collaborating with _c_~B~ people outside of the group. The difference between author _A_ and _B_ is that authors in group _A_ work independently of each other, while authors in group _B_ invite all other group members to be a coauthor on their papers independently of their contribution to that paper (Figure 1). In other words, authors in group _B_ form a !pubcart.

![!FIGURE(1) The publication relationships in group _A_ and _B_. Nodes are authors, while edges symbolise shared publications. Groups of four authors are marked by the underlying shapes. In group _A_ authors work with several coauthors from outside of the group but they do not invite group mates to be coauthors on their own papers. Contrarary, authors in group _B_ form a !pubcart i.e. each author invites all other aouthors in the group to be a coauthor (note the connections between group members).](paperfigs/groups.!EXT)

For simplicity, we initially assume that _p_~A~ = _p_~B~ = _p_, _G_~A~ = _G_~B~ = _G_ and _c_~A~ = _c_~B~ = _c_, i.e. author groups are of the same size, authors produce the same number of papers and they have the same number of coauthors from outside of the group. In this case the total numbers of papers produced by the groups are equal (_Gp_ = _G_~A~_p_~A~ and _G_~B~_p_~B~, respectively). The total numbers of papers (co)authored by author _A_ and _B_ are, however, different. Author _A_ authored _n_~A~ = _p_~A~ = _p_ papers. On the other hand, author _B_ (co)authored  _n_~B~ = _p_~B~ + (_G_~B~ - 1) _p_~B~ = _G_~B~ _p_~B~ = _Gp_ papers. In the case of author _B_ the term (_G_~B~ - 1)_p_~B~ represents the papers on which author _B_ was invited. It is easy to see that as far as _G_ > 1, author _B_ will have more paper than author _A_, i.e _n_~B~ > _n_~A~. 

A natural way to correct for this bias is to taking into account the number of authors each paper has and instead of counting the papers themselves sum the inverse of the number of authors [@Vavrycuk18e0195509]:

$$
w = \sum_{i=1}^{n} \frac{1}{1+C}.
$$

Here, 1 in the denominator symbolises the focal author, while _C_ is the number of coauthors. For author _A_, _C_ = _c_~A~ = _c_. On the other hand, for author _B_, _C_ = (_G_~B~ - 1) + _c_~B~ = (_G_ - 1) + _c_. If _c_ = 0, then the division by the number of coauthors works, we regain the number of papers the authors produced without inviting their group members.

For author _A_:
$$
w_{A} = \sum_{i=1}^{n_A} \frac{1}{1} = \sum_{i=1}^{n_A} 1 = n_A = p.
$$

For author _B_:
$$
w_{B} = \sum_{i=1}^{n_B} \frac{1}{1+G-1} = \sum_{i=1}^{Gp} \frac{1}{G} = \frac{Gp}{G} = p.
$$

On the other hand, if the focal authors collaborate with others outside of their groups, as Figure 1 illustrates, the situation changes (Figure 2):

For author _A_:
$$
w_{A} = \sum_{i=1}^{p}\frac{1}{1+c} = \frac{p}{1+c}.
$$

For author _B_:
$$
w_{B} = \sum_{i=1}^{Gp}\frac{1}{G+c} = \frac{Gp}{G+c}.
$$

The weighted number of papers produced by author _B_ relative to author _A_, _w_~B~/_w_~A~, is:

$$
\frac{w_B}{w_A} = \frac{\frac{Gp}{G+c}}{\frac{p}{1+c}} = \frac{Gp}{G+c} \times \frac{1+c}{p} = \frac{G(1+c)}{G+c} = \frac{G+Gc}{G+c}.
$$

The proportion of _w_~B~/_w_~A~ is greater than one if _G_+_Gc_ > _G_+_c_, which is always true if _c_ > 0 (Figure 2). This means that if authors collaborate anyone from outside their groups then authors in group _B_ will always have higher publication performance than authors in group _A_, despite the fact that the two groups have the same productivity.

![!FIGURE(2) Publication performance when authors collaborate with people from outside of their groups. Weighted publication performance of author _A_ and _B_ (a). Weighted publication performance of author _B_ relative to that of author _A_ (b). The weighted publication performance is calculated by taking into account the number of coauthors. During this calculation first authorship can be rewarded by a bonus, _b_. If _b_ = 0, then each coauthors receive the same weight for a given publication. On the other hand, if _b_ > 0, the weight of the first author is higher then that of the coauthors, i.e. the first author of a paper is rewarded. On subpanel (a) _b_ = 0.2, on (b) _b_ is given on the right margin.](paperfigs/weighted_production.!EXT)

To compensate for this productivity bias, author _A_ should produce _w_~B~/_w_~A~ times more papers, _p_~A~ = _p_~B~(_G_ + _Gc_)/(_G_ + _c_). This surplus of papers needed for compensating the productivity bias increases with _c_ and it keeps to _G_.

Authors in group _A_ can also compensate for the productivity bias by decreasing their collaborators from outside of the group. This reduction must be by a factor of _G_: _c_~A~ = _c_~B~ / _G_.

A useful modification to the weighted performance scheme is the so called _first-author-emphasis_ scheme [@Vavrycuk18e0195509]. In this scheme, the first authors receive a bonus, _b_, to recognise their leading role in producing the papers. Under this scheme the weighted publication performance for author _A_, _w_~A~', is:

$$
w_{A}' = \sum_{i=1}^{n_A} \left(b + \frac{1-b}{1+c_A}\right) = \sum_{i=1}^{p}
\left(b + \frac{1-b}{1+c}\right) = \frac{p(1+bc)}{1+c}.
$$

Here, the first author, who is author _A_ for all his papers, get a bonus _b_ for contributing most to the paper, and the rest of the credit, 1-_b_, is divided equally between all authors (including the first one). The weighted publication performance for author _B_ under the first author scheme, _w_~B~', is:

$$
w_{B}' = \sum_{i=1}^{p_B} \left(b + \frac{1-b}{G_B + c_B}\right) + \sum_{i=1}^{(G_B - 1)p_B} \frac{1-b}{G_B + c_B},
$$

where the first term gives the credit for first author papers, while the second one is for the coauthored papers. After simplification, we obtain:

$$
w_{B}' = \frac{p(G+bc)}{G+c}.
$$

By comparing to _w_~B~' to _w_~A~' it is easy to show that author _B_ will always have a higher publication performance than author _A_, i.e. _w_~B~'/_w_~A~' > 1, if _G_ > 1. Further analysis,

$$
\frac{w_{B}'}{w_{A}'} = \frac{p(G+bc)}{G+c} \times \frac{1+c}{p(1+bc)} = \frac{G+c[G+b(1+c)]}{G+c[1+b(G+c)]},
$$

shows that for _w_~B~'/_w_~A~' > 1, the condition _c_ > 0 should also be fulfilled. As numerical computation indicates (Figure 2) the bias is decreased by introducing the first authorship bonus, but it is still significant. @Vavrycuk18e0195509, for instance, recommend a bonus of _b_ = 0.2, but in this case author _B_ sill has around 50% more credit for the same work than author _A_ has. The difference between author _A_ and _B_ decreases as _b_ increases (Figure 2b), but this way coauthorship is worth less and less, undermining the benefits of collaborations.

To summarise, this simple model shows that the formation of !pubcarts can be an advantageous strategy in terms of increasing publication productivity even if one control for the number of coauthors of papers. Note, however, that this model might be overly simplified as all authors have the same primary productivity and we do not investigated how productivity of authors outside of the cartels changes as a consequence of forming cartels. To obtain a more realistic understanding of forming !pubcarts next we develop a simulation of the publication process.

# The simulations

&define MP _M_~P~
&define MC _M_~C~
&define GC _G_~C~

We start simulating the !pub process with constructing a !pubmat, !MP (Figure 3). The _a_~ij~ element of !MP is one if author _j_ is on the author list of paper _i_ and zero otherwise. Therefore !MP can be considered as a co-occurrence matrix, i.e. a bipartite graph, where rows and columns represent the two types of nodes, papers and authors, respectively. To construct !MP we set up _n_~C~ communities which sizes are given by _c_~i~ (_i_ = 1 ... _n_~C~). In other words, community _i_ consists of _c_~i~ authors. The number of papers written by each author in community _i_ is given by _k_~j~. For each community we first construct an empty matrix (all _a_~ij~ = 0) of size _p_~i~ and _c_~i~, where _p_~i~ > max(_k_~j~). Then for each columns we randomly distributed _k_~j~ number of ones over the _p_~i~ empty places. Having constructed !MP we create a weighted collaboration (or co-authorship) matrix, !MC, by projecting !MP to the nodes of authors. The weights of !MC, _w_~ij~, are Jaccard similarity indices  calculated between each pair of authors _i_ and _j_ (_i_ $\ne$ _j_) as

$$
w_{ij} = \frac{ | A_i \cap A_j | }{ | A_i \cup A_j | }.
$$

Here, _A_~i~ is the set of papers to which author _i_ contributed. In other words, the weight between two authors is the proportion of shared papers to the total number of unique papers to which either of authors _i_ or _j_ contributed to. It varies between zero (i.e. no common publication between author _i_ and _j_) and one (i.e. all publications by the two authors are shared). From !MC we construct a !colgraph, !GC.

![!FIGURE(3) The construction of publication network. The top left panel shows the publication matrix. Each row and column of this matrix represents a paper and an author, respectively. Values of 1s indicate that an author is on the author list of a given paper. From the publication matrix one can derive the collaboration matrix (bottom right panel) by calculating the Jaccard simmilarity (top right) for each possible pairs of authors. The bottom left panel shows the resulting collaboration network, a weighted, undirected graph. The red squares exemplifies the calculation of Jaccar simmilarity.](paperfigs/sample_publication_network-01.!EXT)

We simulated the formation of cartels by randomly choosing $|\kappa|$ authors within a community (Figure 4). Let $\kappa$ is the set of cartel members. Then, with probability _p_~c~, we changed each element _a_~ij~ = 0 of !MP to _a_~ij~ = 1 where the following conditions met: _j_ $\in \kappa$ and at least one _a_~ik~ = 1 with _k_ $\in \kappa$ but _k_ $\ne$ _j_. We project the resulting !pubmat, !MP' to !MC' and constructed the corresponding !colgraph, !GC'.

![!FIGURE(4) The formation of cartels. The panels on the left illustrate a publication network without cartel. The panels on the right show how a cartel between authors A1, A4 and A6 can be formed: Author A6 invites authors A1 and A4 to be coauthors on paper p2, while author A1 do the same with authors A4 and A6 on paper p8. The small red rectangles marks the authorships gained this way. The bottom right panel shows the resulting collaboration graph, where the red edges connect cartel members. Note the strong connections between members.](paperfigs/sample_publication_network-02.!EXT)

# References

