&include "paperfigs/MS_title.md"

<!--
project start date: 2019-10-13
-->

---
standalone: true
bibliography: PaperPump.yaml
natbib: true
biblio-title: References
fontsize: 12pt
papersize: a4paper
header-includes:
- \usepackage{double_spaced}
- \usepackage{lineno}
- \linenumbers
geometry:
- margin=1in
---

<!--
# TODOs


-->

<!--

# Definitions

Proba

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
&define TCG tightly connected group
&define AUTHOR(grp) author _grp_~1~
&define AUTHORS(G1, G2) authors _G1_~1~ and _G2_~1~
&define OPN 1/_n_

-->

# Abstract

The present processes of research assessment, i.e. focusing on one or a few, related, scientometrics, foster questionable authorship practices, like gifting authorship to non-contributing people. An especially harmful one of these unethical practices is the formation of publication cartels, where authors offer gift authorship to each other reciprocally. Here, by developing a simple model and a simulation of the publication process I investigate how beneficial cartels can be and what measure can be used to restrict them. My results indicate that publication cartels can significantly boost members' productivity even if paper counts are weighted by the inverse of author number (the !OPN rule). Nevertheless, applying the !OPN rule generates conflicts of interest both among cartel members themselves and between members and non-members which might lead to the self-purification of the academic publishing industry.


# Introduction

Research integrity [ethical behaviour, sound methodology and rigorous peer review; @szomszor2020] provides assurance that scientific activities lead to trustable and replicable results. Research integrity is, however, under threat as a result of how science currently operates. The recent, unprecedented expansion of science, exemplified, for instance, by the  exponentially growing number of scientific articles [@fire2019], gives way to the wide-spread use of scientometry for assessing the productivity and impact of researchers [@aubertbonn2021]. As science is usually funded by public resources, the desire to measure the performance of its actors is well justified. Introducing the assessment of scientists by one or a few metrics, like number of publications or citations, together with the hyper-competitiveness of science had, however, somehow unexpected consequences [@biagioli2019].

As, among others, Charles Goodhart observed, if a metric is used as a target then it becomes a bad metric [@edwards2017; @fire2019]. This happens because people, in response to introduction of a target, alter their own behaviour to affect the metric directly instead to modify the activity the change of which was intended by introducing the metric [@werner2015]. In the recent process of corporatisation of science two such metrics became relevant: the numbers of papers and citations [@grossman2019].

Goodhart's law is well illustrated by the introduction of the number of papers as a measure of productivity in science. Using this measure is based on the assumption that characteristics of scientific papers (like length or number of coauthors) are fixed and hence targeting more papers automatically leads to the generation of more new knowledge. Unfortunately, this was not what had happened, scientists responded in some unexpected, nevertheless clearly rational but sometimes unethical, ways [@fong2017; @gopalakrishna2021]. For instance, they reduced the length of papers [@fire2019], i.e. they are publishing the same amount of knowledge in more papers (salami articles). Furthermore, mangling with authorship appeared where offering authorship to those who did not contributed to the given paper considerably (honorary authorship) can quickly increase their number of publications, again without any increase in knowledge produced [@aubertbonn2021; @biagioli2019; @fong2017; @gopalakrishna2021]. A possible sign of this questionable authorship practice can be the recent raise of number of authors per paper [@fire2019]. One may argue that more authors per paper is the sign of science becoming more interdisciplinary. A recent analysis is, however, unlikely to support this conclusion; the number of coauthors increases with time even after controlling for attributes related to complexity of science [@papatheodorou2008]. Another reason for the increased number of coauthors might be the increased efficiency that can follow from the increased possibility for division of labour facilitated by more authors [@demesnard2017]. In this case, however, it is expected that the number of papers per author also increases, which seems not to be the case [@fire2019].

Questionable authorship practice, on the other hand, appears to be common. Recent surveys suggest that about 30% of authors were involved in these unethical practices [@biagioli2019; @fong2017; @gopalakrishna2021; @grossman2019; @halaweh2020; @marusic2011]. One of these practices is ghost authorship when someone who has significantly contributed to the article is excluded from the author bylist [@jabbehdari2017]. In other forms (honorary authorship) just the opposite happens; those are offered authorship who have not (considerably) contributed to the work published [@fong2017; @gopalakrishna2021]. Several reasons can be behind gifting authorship to someone. Junior authors might include more senior ones because of respect or they are forced to do so [@pan2020]. Senior authors may gift authorship to juniors in order to help them obtain post-doctoral scholarships or tenure [@vonbergen2017].

A very efficient way to increase the number of publications may be to practice honorary authorship reciprocally. The most organised form of this behaviour is founding publication cartels. The cartel is formed by a group of people who agree to mutually invite each others to their own publications as guest authors without any contribution. As in recent assessment practice coauthored papers count as a whole publication to every coauthor on the bylist, publication cartels can significantly boost the productivity of cartel members. This is the phenomenon which is called as 'publication club' by @demesnard2017. As the noun of 'club' involves a positive connotation I prefer to use 'cartel' for this under studied but highly unethical behaviour. Simple argument suggests that sharing the credit of a publication among the coauthors can decrease the incentive of forming cartels [@demesnard2017]. The simplest scenario for sharing is the 1/_n_ rule under which only 1/_n_ part of a publication is attributed to each of the _n_ coauthors of the given paper [@demesnard2017].

In this paper I develop a simple model of publication cartels to understand how effective they are to increase members' productivity and whether it is possible to eliminate them by applying different measures, like the 1/_n_ rule. I then extend my study to situations resembling more to real world conditions by developing a computer simulation of cartels. I use this simulation to investigate how using different metrics of productivity affect authors outside of cartels.


# The model

We compare the publication performance of two authors, !AUTHOR(A) and !AUTHOR(B). Authors work in separate groups (group _A_ and group _B_, respectively) each of which contains _G_~i~ (_i_ = _A_ or _B_) people (including the focal author). Each author in group _A_ produces _p_~A~ papers in a year by collaborating with _c_~A~ authors from outside of the group, i.e. their primary production is _p_~A~. Similarly, each author in group _B_ primarily produces _p_~B~ papers by collaborating with _c_~B~ people outside of the group. The difference between !AUTHORS(A,B) is that authors in group _A_ work independently of each other, while authors in group _B_ invite all other group members to be a coauthor on their papers independently of their contribution to that paper (Fig 1). In other words, authors in group _B_ form a !pubcart.

![!FIGURE(1) The publication relationships in groups _A_ and _B_ of the model. Nodes are authors, while edges symbolise shared publications. Groups of four authors are marked by the underlying shapes. In group _A_ authors work with several coauthors from outside of the group but they do not invite group mates to be coauthors on their own papers. Contrarary, authors in group _B_ form a !pubcart i.e. each author invites all other authors in the group to be a coauthor (note the connections between group members).](paperfigs/groups.!EXT)

For simplicity, we assume that _G_~A~ = _G_~B~ = _G_ (_G_ > 1), _p_~A~ = _p_~B~ = _p_ and _c_~A~ = _c_~B~ = _c_, i.e. author groups are of the same size, authors produce the same number of primary papers and they have the same number of coauthors from outside of the group. In this case the total numbers of papers produced by the groups, the group productivity, are equal (_Gp_ = _G_~A~_p_~A~ and _G_~B~_p_~B~, respectively). The total numbers of papers (co)authored by !AUTHORS(A,B) are, however, different. Author _A_~1~ writes _n_~A~ = _p_~A~ = _p_ papers. On the other hand, !AUTHOR(B) (co)authors  _n_~B~ = _p_~B~ + (_G_~B~ - 1) _p_~B~ = _G_~B~ _p_~B~ = _Gp_ papers. In the case of !AUTHOR(B) the term (_G_~B~ - 1)_p_~B~ represents the papers on which !AUTHOR(B) is invited as honorary author. It is easy to see that as far as _G_ > 1, !AUTHOR(B) will have many more paper than !AUTHOR(A), i.e _n_~B~ > _n_~A~. 

A natural way to correct for this bias is to taking into account the number of authors each paper has and instead of counting the papers themselves as a measure of productivity one sums the inverse of the number of authors [the 1/_n_ rule, @demesnard2017; @vavrycuk2018]:

$$
w = \sum_{i=1}^{n} \frac{1}{1+C}.
$$

Here, number 1 in the denominator symbolises the focal author, while _C_ is the number of coauthors. For !AUTHOR(A), _C_ = _c_~A~ = _c_. On the other hand, for !AUTHOR(B), _C_ = (_G_~B~ - 1) + _c_~B~ = (_G_ - 1) + _c_. If _c_ = 0, then the division by the number of coauthors works, we regain the number of papers the authors produced without inviting their group members.

For !AUTHOR(A):
$$
w_{A} = \sum_{i=1}^{n_A} \frac{1}{1} = \sum_{i=1}^{n_A} 1 = n_A = p.
$$

For !AUTHOR(B):
$$
w_{B} = \sum_{i=1}^{n_B} \frac{1}{1+G-1} = \sum_{i=1}^{Gp} \frac{1}{G} = \frac{Gp}{G} = p.
$$

On the other hand, if the focal authors collaborate with others outside of their groups, as Fig 1 illustrates, the situation changes (Fig 2):

For !AUTHOR(A):
$$
w_{A} = \sum_{i=1}^{p}\frac{1}{1+c} = \frac{p}{1+c}.
$$

For !AUTHOR(B):
$$
w_{B} = \sum_{i=1}^{Gp}\frac{1}{G+c} = \frac{Gp}{G+c}.
$$

The weighted number of papers produced by !AUTHOR(B) relative to !AUTHOR(A), _w_~B~/_w_~A~, is:

$$
\frac{w_B}{w_A} = \frac{\frac{Gp}{G+c}}{\frac{p}{1+c}} = \frac{Gp}{G+c} \times \frac{1+c}{p} = \frac{G(1+c)}{G+c} = \frac{G+Gc}{G+c}.
$$

The proportion of _w_~B~/_w_~A~ is greater than one if _G_+_Gc_ > _G_+_c_, which is always true if _c_ > 0 (as we already assumed _G_ > 1; Fig 2). This means that if authors collaborate anyone from outside of their groups then authors in group _B_ will always have higher publication performance than authors in group _A_, despite the fact that the two groups have the same productivity.

![!FIGURE(2) Publication performance when authors collaborate with people from outside of their groups. Weighted publication performance of !AUTHORS(A,B) (a). Weighted publication performance of !AUTHOR(B) relative to that of !AUTHOR(A) (b). The weighted publication performance is calculated by taking into account the number of coauthors. During this calculation first authorship can be rewarded by a bonus, _b_. If _b_ = 0, then each coauthors receive the same weight for a given publication. On the other hand, if _b_ > 0, the weight of the first author is higher then that of the coauthors, i.e. the first author of a paper is rewarded. On subpanel (a) _b_ = 0.2, on (b) _b_ is given on the right margin.](paperfigs/weighted_production.!EXT)

To compensate for this productivity bias, !AUTHOR(A) should produce _w_~B~/_w_~A~ times more papers, _p_~A~ = _p_~B~(_G_ + _Gc_)/(_G_ + _c_). This surplus of papers needed for compensating the productivity bias increases with _c_ and it keeps to _G_.

Authors in group _A_ can also compensate for the productivity bias by decreasing the number of their collaborators from outside of the group. This reduction must be by a factor of _G_: _c_~A~ = _c_~B~ / _G_.

A useful modification to the 1/_n_ rule is the so called _first-author-emphasis_ scheme [@vavrycuk2018]. In this scheme, the first authors receive a bonus, _b_, to recognise their leading role in producing the papers. Under this scheme the weighted publication performance for !AUTHOR(A), _w_~A~', is:

$$
w_{A}' = \sum_{i=1}^{n_A} \left(b + \frac{1-b}{1+c_A}\right) = \sum_{i=1}^{p}
\left(b + \frac{1-b}{1+c}\right) = \frac{p(1+bc)}{1+c}.
$$

Here, the first author, who is author _A_ for all his papers, get a bonus _b_ for contributing most to the paper, and the rest of the credit, 1-_b_, is divided equally between all authors [including the first author, @vavrycuk2018]. The weighted publication performance for !AUTHOR(B) under the first author scheme, _w_~B~', is:

$$
w_{B}' = \sum_{i=1}^{p_B} \left(b + \frac{1-b}{G_B + c_B}\right) + \sum_{i=1}^{(G_B - 1)p_B} \frac{1-b}{G_B + c_B},
$$

where the first term gives the credit for first author papers, while the second one is for the coauthored papers. After simplification, we obtain:

$$
w_{B}' = \frac{p(G+bc)}{G+c}.
$$

By comparing _w_~B~' to _w_~A~' it is easy to show that !AUTHOR(B) will always have a higher publication performance than !AUTHOR(A), i.e. _w_~B~'/_w_~A~' > 1, if _G_ > 1 and _b_ < 1. Further analysis,

$$
\frac{w_{B}'}{w_{A}'} = \frac{p(G+bc)}{G+c} \times \frac{1+c}{p(1+bc)} = \frac{G+c[G+b(1+c)]}{G+c[1+b(G+c)]},
$$

shows that for _w_~B~'/_w_~A~' > 1, the condition _c_ > 0 should also be fulfilled. As numerical computation indicates (Fig 2) the bias is decreased by introducing the first authorship bonus, but it is still significant. @vavrycuk2018, for instance, recommend a bonus of _b_ = 0.2, but in this case !AUTHOR(B) sill has around 50% more credit for the same work than !AUTHOR(A) has. The difference between !AUTHORS(A,B) decreases as _b_ increases (Fig 2b), but this way coauthorship is worth less and less, undermining the possible benefits of collaborations.

To summarise, this simple model shows that the formation of !pubcarts can be an advantageous, but unethical, strategy to increase publication productivity even if one control for the number of coauthors of papers. Note, however, that this model might be overly simplified as all authors have the same primary productivity and we do not investigated how productivity of authors outside of the cartels changes as a consequence of founding cartels. To obtain a more realistic understanding of !pubcarts next I develop a simulation of the publication process.

# The simulation

&define MP _M_~P~
&define MC _M_~C~
&define GC _G_~C~

We start simulating the !pub process with constructing a !pubmat of papers and authors, !MP (Fig 3). Element _a_~ij~ of !MP is one if author _j_ is on the bylist of paper _i_ and zero otherwise. Therefore, !MP can be considered as a matrix representation of a bipartite graph, where rows and columns represent the two types of nodes, papers and authors, respectively. To construct !MP we consider a community of _c_ authors. The number of papers written by author _j_ in the community is given by _k_~j~. For the community we construct an empty matrix (all _a_~ij~ = 0) of size _p_ and _c_, where _p_ > max(_k_~j~). Then for each column _j_ we randomly distributed _k_~j~ number of ones over the _p_ empty places. Having constructed !MP we create a weighted collaboration (or co-authorship) matrix, !MC, by projecting !MP to the nodes of authors. The weights of !MC, _J_~ij~, are Jaccard similarity indices  calculated between each pair of authors _i_ and _j_ (_i_ $\ne$ _j_) as

$$
J_{ij} = \frac{ | P_i \cap P_j | }{ | P_i \cup P_j | }.
$$

Here, _P_~i~ is the set of papers to which author _i_ contributed. In other words, the weight between two authors is the proportion of shared papers to the total number of unique papers to which either of authors _i_ or _j_ contributed to. It varies between zero (i.e. no common publication between author _i_ and _j_) and one (i.e. all publications by the two authors are shared). Note, Jaccard similarity between authors in group A of the above model is zero, while between authors in group B is one. From !MC we construct a !colgraph, !GC.

![!FIGURE(3) The construction of publication network. The top left panel shows the publication matrix, !MP. Each row and column of this matrix represents a paper and an author, respectively. Values of 1 indicate that an author is on the author list of a given paper, while dots symbolise zeros. From the publication matrix one can derive the collaboration matrix, !MC (bottom right panel) by calculating the Jaccard simmilarity !TODO(update the symbols in the graph for Jaccard index) (top right) for each possible pairs of authors. The bottom left panel shows the resulting weighted, undirected collaboration graph, !GC. The red rectangles exemplifies the calculation of Jaccar simmilarity.](paperfigs/sample_publication_network-01.!EXT)

I simulated the formation of cartels by choosing $|\kappa|$ authors from the community (Fig 4). Let $\kappa$ is the set of cartel members. Then, with probability _p_~c~, I changed each element _a_~ij~ = 0 of !MP to _a_~ij~ = 1 where the following conditions met: _j_ $\in \kappa$ and at least one _a_~ik~ = 1 with _k_ $\in \kappa$ but _k_ $\ne$ _j_. I project the resulting !pubmat, !MP' to !MC' and constructed the corresponding !colgraph, !GC'.

![!FIGURE(4) The formation of cartels. The panels on the left illustrate a publication network without cartel. The panels on the right show how a cartel between authors A1, A4 and A6 can be formed: Author A6 invites authors A1 and A4 to be coauthors on paper p2, while author A1 do the same with authors A4 and A6 on paper p8. The small red rectangles mark the authorships gained this way. The bottom right panel shows the resulting collaboration graph, where the red edges connect cartel members. Note (i) the strong connections between members and (ii) adding cartels also changes the connections of non-members.](paperfigs/sample_publication_network-02.!EXT)

By setting all _k_~j~ = _k_ and _p_ > > _k_ we can simulate the case of equal productivity and no collaboration from outside of the group. Here the simulation produces the same results as the model: productivity of cartel members increased but this can be accounted for by using weighted number of publications.

To induce collaboration between authors I next set _p_ < $\sum$_k_ (the authors still have the same productivity prior to cartel formation). Under these conditions, if we consider the number of papers, the productivity of cartel members increases significantly by forming cartel while productivity of non-members does not change (Fig 5). In accordance with the model, the productivity of cartel members increases even if we consider the weighted number of papers. Interestingly, the productivity of many non-members decreases when cartel is formed (Fig 5).

![!FIGURE(5) The effect of cartel formation on the productivity of cartel members and non-members: equal prior productivity of authors. The top panels illustrate the collaboration graph before and after cartel formation. The middle panels show how the number of papers produced by members and non-members changes because of founding cartel. The bottom panels illustrate the same but using the weighted number of papers as a measure of productivity. Collaboration graph formed with _c_ = 30, _k_ = 3, _p_ = 60, _p_~c~ = 1 and $\kappa =\{1,2,3,29,30\}$.](paperfigs/simulation_analyses-01.!EXT)

I further generalise the simulation results by setting the prior productivity of authors to different values (Fig 6). Using the number of papers as metric leads to the same conclusions: members' productivity increases after cartel formation, non-members' productivity does not change. On the other hand, using the weighted number of papers reveal an interesting effect: the productivity of cartel members with high prior productivity have their productivity being decreased because of cartel foundation (Fig 6). Similarly to the previous case, productivity of many non-members decreases as a consequence of cartel formation (Fig 6).

![!FIGURE(6) The effect of cartel formation on the productivity of cartel members and non-members: prior productivity of authors differs. The top panels illustrate the collaboration graph before and after cartel formation. The middle panels show how the number of papers produced by members and non-members changes because of founding cartel. The bottom panels illustrate the same but using the weighted number of papers as a measure of productivity. Collaboration graph formed with _c_ = 30, _k_~jj~ = j, _p_ = 60, _p_~c~ = 1 and $\kappa =\{1,2,3,29,30\}$.](paperfigs/simulation_analyses-02.!EXT)

# Conclusions

<!--
Both the model and the simulation presented here indicate under a rather wide range of conditions that members' productivity always increases by founding cartels. 

As under the current circumstances found in the academic world most measures of productivity/impact (e.g. citation counts, _h_-index) are based on or correlates strongly with the unweighted number of papers an author produces there is a huge pressure on scientists to establish cartels under the current circumstances found in the academic world, where most measures of productivity, e.g. citation counts, _h_-index, is based on the unweighted number of papers an author produces [@aubertbonn2021; @fong2017; @grossman2019]. Unfortunately, switching to the use of weighted number of papers does not solve the problems arising from the possibility of cartel formation, the productivity of members increases in most of the cases when they found cartels. Despite of this result, using weighted number of papers as metric can offer a significant opportunity to reduce unethical authorship practice, because cartel formation reduces the productivity of out-of-cartel collaborators, generating a conflict of interest between cartel members and non-members [@demesnard2017]. This might lead to a situation when cartel members become isolated without the advantages of out-of-group collaborations and hence resulting a significant decrease in the benefits of cartel formation. The results of the simulation also suggest that, given that weighted number of papers used to rate authors, cartel members should be of similar prior productivity, because cartel establishment with lowly ranked authors can have a detrimental effect on the productivity of prolific authors. This would also introduce a conflict of interest between authors leading to that founding cartels is only worth for low productivity authors among themselves, because (i) highly ranked authors do not need to manipulate their productivity indices (they are already high) and (ii) they can actually loose on cartel formation, given that productivity is measured by weighted number of papers. As currently cartel formation cannot be penalised in science [@fister2016] we highly recommend to switch to use weighted number of papers as a basis of measuring productivity, because it introduces conflicts of interest that might drive the publication system towards self-purification.
-->

Under the current climate of wide spread use of scientometry indices to assess academics publication cartels can provide huge, although unethical benefits. As my results indicate, members of cartels by reciprocally inviting each other as honorary authors can easily boost their own publication productivity, i.e. the number of papers they appear on as (co)author. As many scientometrics currently in use are strongly associated with the number of publications a scholar has produced [@aubertbonn2021; @fong2017; @grossman2019] becoming cartel member can have a very general positive effect on one's academic career.

One may consider that fighting off cartels is not necessary because of "no harm no foul": research integrity may not be inevitably damaged by cartel foundation, cartels can produce high quality research. Nevertheless, cartels do distort the research competition landscape. This might result in that highly competent, talented researchers, who are not members of any cartels, are forced into inferior roles which, in turn, compromises the society's ability to produce more novel and innovative results. Therefore, cartel formation should be restricted.

Fighting against cartels is, however, not trivial. First, identifying cartels, not to mention to prove for a group of researchers that they are cartelling is inherently difficult. Investigating properties of coauthor networks might help. Nevertheless, a possible way to restrict cartels without their identification is to use such scientometrics which penalise cartel formation. An obvious choice can be to weight the number of publications an author has by the inverse of the number of authors on the bylists of these papers, the so-called 1/_n_ rule. As my calculation shows this rule can only be effective if coauthorship only occurs between cartel members. As soon as collaboration is wide spread among both cartelling and non-cartelling authors my results indicate that the !OPN rule breaks down and cartel members still gain undeserved benefits. On the other hand, my computations also show that the !OPN rule can still be useful against publication cartels, because it generates conflicts of interest among the parties. Collaborators of cartel members suffer a loss if the !OPN rule is applied which might force them either to change the unethical behaviour of cartel members or abandon to collaborate with them.

The results of the simulation also suggest that, given that the !OPN rule is used to rate authors, cartel members should be of similar prior productivity, because cartel establishment with lowly ranked authors can have a detrimental effect on the productivity of prolific authors. This would also introduce a conflict of interest between authors leading to that founding cartels is only worth for low productivity authors among themselves, because (i) highly ranked authors do not need to manipulate their productivity indices (they are already high) and (ii) they can actually loose on cartel formation, given that productivity is measured by weighted number of papers, i.e. according to the !OPN rule. 

To summarise, I strongly argue for using the !OPN rule as the basis of scientometry. Unfortunately, its general use is opposed by many parties for many reasons. It still remains to see whether these reasons are valid or not, but my calculations indicate that the application of !OPN rule can generate such processes which may ultimately lead to the self-purification of the academic publication industry. Of course, abandoning the current, metric-only research assessment system can also help.

# Acknowledgments

I thank Miklós Bán, Gábor Lövei, Tibor Magura and Jácint Tökölyi to review a previous version of the manuscript. The research was supported by the Thematic Excellence Programme (TKP2020-IKA-04) of the Ministry for Innovation and Technology in Hungary.

# References
