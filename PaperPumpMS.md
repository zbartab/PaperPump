---
title: "Publication games: in the web of reciprocity"
author: Z Barta
date: !__DATE__
---
<!--
project start date: 2019-10-13
-->

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

We compare the publication performance of two authors, author _A_ and author _B_. Authors work in separate groups (group _A_ and group _B_, respectively) each of which contains _G_~i~ (_i_ = _A_ or _B_) people (including the focal author). Each author in group _A_ produces _p_~A~ papers in a year by collaborating with _c_~A~ authors from outside of the group, i.e. his/her prinary production is _p_~A~. Similarly, each author in group _B_ primarily produces _p_~B~ papers by collaborating with _c_~B~ people outside of the group. The difference between author _A_ and _B_ is that authors in group _A_ work independently of each other, while authors in group _B_ invite all other group members to be a coauthor on their papers independently of their contribution to that paper (Figure 1). In other words, authors in group _B_ form a !pubcart.

![!FIGURE(1) The publication relationships in group _A_ and _B_. Nodes are authors, while edges symbolise shared publications. Groups of four authors are marked by the underlying shapes. In group _A_ authors work with several coauthors from outside of the group but they do not invite group mates to be coauthors on their own papers. Contrarary, authors in group _B_ form a !pubcart i.e. each author invites all other aouthors in the group to be a coauthor (note the connections between group members).](paperfigs/groups.!EXT)

For simplicity, we assume that _p_~A~ = _p_~B~ = _p_, _G_~A~ = _G_~B~ = _G_ (_G_ > 1) and _c_~A~ = _c_~B~ = _c_, i.e. author groups are of the same size, authors produce the same number of primary papers and they have the same number of coauthors from outside of the group. In this case the total numbers of papers produced by the groups are equal (_Gp_ = _G_~A~_p_~A~ and _G_~B~_p_~B~, respectively). The total numbers of papers (co)authored by author _A_ and _B_ are, however, different. Author _A_ writes _n_~A~ = _p_~A~ = _p_ papers. On the other hand, author _B_ (co)authors  _n_~B~ = _p_~B~ + (_G_~B~ - 1) _p_~B~ = _G_~B~ _p_~B~ = _Gp_ papers. In the case of author _B_ the term (_G_~B~ - 1)_p_~B~ represents the papers on which author _B_ is invited. It is easy to see that as far as _G_ > 1, author _B_ will have more paper than author _A_, i.e _n_~B~ > _n_~A~. 

A natural way to correct for this bias is to taking into account the number of authors each paper has and instead of counting the papers themselves as a measure of productivity sum the inverse of the number of authors [@Vavrycuk18e0195509]:

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

Authors in group _A_ can also compensate for the productivity bias by decreasing the number of their collaborators from outside of the group. This reduction must be by a factor of _G_: _c_~A~ = _c_~B~ / _G_.

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

By comparing _w_~B~' to _w_~A~' it is easy to show that author _B_ will always have a higher publication performance than author _A_, i.e. _w_~B~'/_w_~A~' > 1, if _G_ > 1. Further analysis,

$$
\frac{w_{B}'}{w_{A}'} = \frac{p(G+bc)}{G+c} \times \frac{1+c}{p(1+bc)} = \frac{G+c[G+b(1+c)]}{G+c[1+b(G+c)]},
$$

shows that for _w_~B~'/_w_~A~' > 1, the condition _c_ > 0 should also be fulfilled. As numerical computation indicates (Figure 2) the bias is decreased by introducing the first authorship bonus, but it is still significant. @Vavrycuk18e0195509, for instance, recommend a bonus of _b_ = 0.2, but in this case author _B_ sill has around 50% more credit for the same work than author _A_ has. The difference between author _A_ and _B_ decreases as _b_ increases (Figure 2b), but this way coauthorship is worth less and less, undermining the benefits of collaborations.

To summarise, this simple model shows that the formation of !pubcarts can be an advantageous, but unethical, strategy in terms of increasing publication productivity even if one control for the number of coauthors of papers. Note, however, that this model might be overly simplified as all authors have the same primary productivity and we do not investigated how productivity of authors outside of the cartels changes as a consequence of founding cartels. To obtain a more realistic understanding of !pubcarts next we develop a simulation of the publication process.

# The simulation

&define MP _M_~P~
&define MC _M_~C~
&define GC _G_~C~

We start simulating the !pub process with constructing a !pubmat of papers and authors, !MP (Figure 3). Element _a_~ij~ of !MP is one if author _j_ is on the author list of paper _i_ and zero otherwise. Therefore, !MP can be considered as a matrix representation of a bipartite graph, where rows and columns represent the two types of nodes, papers and authors, respectively. To construct !MP we set up _n_~C~ communities sizes of which are given by _c_~l~ (_l_ = 1 ... _n_~C~), i.e. community _l_ consists of _c_~l~ authors. The number of papers written by each author in community _l_ is given by _k_~j~. For each community we first construct an empty matrix (all _a_~ij~ = 0) of size _p_~l~ and _c_~l~, where _p_~l~ > max(_k_~j~). Then for each columns _j_ we randomly distributed _k_~j~ number of ones over the _p_~l~ empty places. When we have constructed a submatrix for each community, we combine them along their diagonal to form the whole !pubmat by filling the appearing off-diagonal empty cells with zeros. Finaly, we randomly switch 10% of the ones with zeros such that the marginal sums of rows and columns do not change. This degree maintianing rewiring connects the separate communities into a more-or-less connected graph. Having constructed !MP we create a weighted collaboration (or co-authorship) matrix, !MC, by projecting !MP to the nodes of authors. The weights of !MC, _w_~ij~, are Jaccard similarity indices  calculated between each pair of authors _i_ and _j_ (_i_ $\ne$ _j_) as

$$
w_{ij} = \frac{ | A_i \cap A_j | }{ | A_i \cup A_j | }.
$$

Here, _A_~i~ is the set of papers to which author _i_ contributed. In other words, the weight between two authors is the proportion of shared papers to the total number of unique papers to which either of authors _i_ or _j_ contributed to. It varies between zero (i.e. no common publication between author _i_ and _j_) and one (i.e. all publications by the two authors are shared). From !MC we construct a !colgraph, !GC.

![!FIGURE(3) The construction of publication network. The top left panel shows the publication matrix, !MP. Each row and column of this matrix represents a paper and an author, respectively. Values of 1s indicate that an author is on the author list of a given paper, while dots symbolise zeros. From the publication matrix one can derive the collaboration matrix, !MC (bottom right panel) by calculating the Jaccard simmilarity (top right) for each possible pairs of authors. The bottom left panel shows the resulting weighted, undirected collaboration graph, !GC. The red squares exemplifies the calculation of Jaccar simmilarity.](paperfigs/sample_publication_network-01.!EXT)

We simulated the formation of cartels by randomly choosing $|\kappa|$ authors several times within a community (Figure 4). Let $\kappa$ is the set of cartel members. Then, with probability _p_~c~, we changed each element _a_~ij~ = 0 of !MP to _a_~ij~ = 1 where the following conditions met: _j_ $\in \kappa$ and at least one _a_~ik~ = 1 with _k_ $\in \kappa$ but _k_ $\ne$ _j_. We project the resulting !pubmat, !MP' to !MC' and constructed the corresponding !colgraph, !GC'.

![!FIGURE(4) The formation of cartels. The panels on the left illustrate a publication network without cartel. The panels on the right show how a cartel between authors A1, A4 and A6 can be formed: Author A6 invites authors A1 and A4 to be coauthors on paper p2, while author A1 do the same with authors A4 and A6 on paper p8. The small red rectangles mark the authorships gained this way. The bottom right panel shows the resulting collaboration graph, where the red edges connect cartel members. Note (i) the strong connections between members and (ii) adding cartels also changes the connections of non-members.](paperfigs/sample_publication_network-02.!EXT)

By setting all _k_~j~ = _k_ and _p_~l~ >> _k_ we can simulate the case of equal productivity and no collaboration from outside of the group. Here the simulation produces the same results as the model: productivity of cartel members increased but this can be accounted for by using weighted number of publications.

To induce collaboration between authors we next set _p_~l~ < $\sum$_k_~j~ (the authors still have the same productivity prior to cartel formation). Under these conditions, if we consider the number of papers, the productivity of cartel members increases significantly by forming cartel while productivity of non-members does not change (Figure 5). In accordance with the model, the productivity of cartel members increases even if we consider the weighted number of papers. Interestingly, the productivity of many non-members decreases when cartel is formed (Figure 5).

![!FIGURE(5) The effect of cartel formation on the productivity of cartel members and non-members: equal prior productivity of authors. The top panels illustrate the collaboration graph before and after cartel formation. The middle panels show how the number of papers produced by members and non-members changes because of founding cartel. The bottom panels illustrate the same but using the weighted number of papers as a measure of productivity. Collaboration graph formed with _n_~C~ = 1, _c_~l~ = 30, _k_ = 3, _p_~l~ = 60, _p_~c~ = 1 and $\kappa =\{1,2,3,29,30\}$.](paperfigs/simulation_analyses-01.!EXT)

We further generalise the simulation results by setting the prior productivity of authors to different values (Figure 6). Using the number of papers leads to the same conclusions: members' productivity increases after cartel formation, non-members' productivity does not change. On the other hand, using the weighted number of papers reveal an interesting effect: the productivity of members with high prior productivity have their productivity being decreased because of cartel foundation (Figure 6). Similarly to the previous case, productivity of many non-members decreases as a consequence of cartel formation (Figure 6).

![!FIGURE(6) The effect of cartel formation on the productivity of cartel members and non-members: prior productivity of authors differs. The top panels illustrate the collaboration graph before and after cartel formation. The middle panels show how the number of papers produced by members and non-members changes because of founding cartel. The bottom panels illustrate the same but using the weighted number of papers as a measure of productivity. Collaboration graph formed with _n_~C~ = 1, _c_~l~ = 30, _k_~j~ = j, _p_~l~ = 60, _p_~c~ = 1 and $\kappa =\{1,2,3,29,30\}$.](paperfigs/simulation_analyses-02.!EXT)

The simulations strengthens the results of the model above, namely under real world conditions (i.e. productivity of authors differs, out-of-group collaborations) the number of papers as a measure of productivity cannot correct for the formation of cartels; members' productivity always increases by founding cartels. This effect, because the great benefit in terms of productivity it offers, can put a huge pressure on scientists to establish cartels under the current circumstances found in the academic world, where most measures of productivity, e.g. citation counts, _h_-index, is based on the unweighted number of papers an author produces!TODO(we need a citation for this). Unfortunately, switching to the use of weighted number of papers does not solve the problems arising from the possibility of cartel formation, the productivity of members increases in most of the cases when they found cartels. Despite of this result, using weighted number of papers can offer a significant opportunity to reduce unethical publication practice, because cartel formation reduces the productivity of out-of-cartel collaborators, generating a conflict of interest between cartel members and non-members. This might lead to a situation when cartel members become isolated without the advantages of out-of-group collaborations and hence resulting a significant decrease in the benefits of cartels formation. The results of the simulation also suggest that, given that weighted number of papers used to rate authors, cartel members should be of similar prior productivity, because cartel establishment with lowly ranked authors can have a detrimental effect on the productivity of prolific authors. This would also introduce a conflict of interest between authors leading to that founding cartels is only worth for low productivity authors among themselves, because (i) highly ranked authors do not need to manipulate their productivity indices (they are already high) and (ii) they can actually loose on cartel formation, given that productivity is measured by weighted number of papers. As currently cartel formation cannot be penalised in science!TODO(we need a citation for this) we highly recommend to switch to use weighted number of papers as a basis of measuring of productivity, because it introduces conflicts of interest that might drive the publication system towards self-purification.

# The prevalence of cartel formation

Given the significant benefit founding cartels provides, naturally arise the questions of (i) whether cartels can be identified in real !pubnets and (ii) how prevalent they are there. To answer these questions we first perform simulation studies then investigates real !pubnets for the sign of cartel formation.

## Footprints of cartels

We use the procedure outlined above to simulate !pubnets with the following differences: Now we simulate many communities and many more authors. The community sizes, _c_~l~, are drawn from a saturated power-law distribution with $\gamma$ = 2.5 and saturation point of 10 [@AlbertLaszloBarabasi16]. The primary productivity of authors follows a saturated power-law distribution with an exponential cutoff with parameters $\gamma$ = 2.5, _k_~sat~ = 10 and _k_~cut~ = 450 [@AlbertLaszloBarabasi16]. The number of potential papers in community _l_, _p_~l~, is given as _p_~l~ = 6 _c_~l~. After the combination of the simulated communities we apply a degree-preserving randomisation [@AlbertLaszloBarabasi16] to rewire 10% of the links between authors and papers and hence connecting the communities together. We refer these networks as "random" networks hereafter. After generating the random networks we forms cartels in each community. The size of these cartels follows a power-law distribution ($\gamma$ = 2.5) truncated at 20. We added one to these values to avoid cartels of size one. Around 0.1% of the authors in each community take part in the cartels.

The formation of cartels increases 

![!FIGURE(7) The effects of founding cartels. The panels show the distribution of link weihts, node strengthes and the number of authors per papers (from top to bottom). Orange lines show the empirical complementary cummulative distribution function for ten separate simulations (dashed lines) and the combination of these runs (solid lines). Green lines illustrate the same for publication networks where cartels added (for details see main text).](paperfigs/cartel_footprint-01.!EXT)

# References
