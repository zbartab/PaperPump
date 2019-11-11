% Why correcting for coauthorships fail?
% Z Barta
% 2019-10-11

---
toc: true
bibliography: /home/apa/notes/articles.bib
natbib: true
biblio-files: /home/apa/notes/articles.bib
biblio-style: /home/apa/lib/texinputs/AmNat
biblio-title: References
fontsize: 12pt
papersize: a4paper
include-before:
- \linenumbers
header-includes:
- \usepackage{lineno}
- \usepackage{double_spaced}
geometry:
- margin=1in
---

# TODOs

- write introduction

# Introduction

# The model

We compare the publication performance of two authors, author _A_ and author _B_. Authors work in separate groups (group _A_ and group _B_, respectively) each of which contains _G_~i~ (_i_ = _A_ or _B_) people (including the focal author). Each author in group _A_ produces _a_~A~ papers in a year by collaborating with _c_~A~ authors from out side of the group. Similarly, authors in group _B_ produce _a_~B~ papers by collaborating with _c_~B~ people outside of the group. The difference between author _A_ and _B_ is that authors in group _A_ work independently of each other, while authors in group _B_ invite all other group members to be a coauthor on their papers independently of their role in producing that paper. Figure 1 illustrates the relationships between authors and coauthors in the groups.

![The publication relationships in group _A_ and _B_. Nodes symbolise authors, while edges show shared publications. Groups of four authors are marked by the underlying shapes. In group _A_ authors work with several coauthors from outside of the group but they do not invite group mates to be coauthors on their own papers. Contrarary, in group _B_ each author invites all other aouthors in the group to be a coauthor (note the connections between group members).](groups.png)

For simplicity, we assume that _a_~A~ = _a_~B~ = _a_, _G_~A~ = _G_~B~ = _G_ and _c_~A~ = _c_~B~, i.e. author groups are of the same size, authors produce the same number of papers and they have the same number of coauthors from outside of the group. In this case the total numbers of papers produced by the groups are equal (_Ga_ = _G_~A~_a_~A~ and _G_~B~_a_~B~, respectively). The total numbers of papers (co)authored by author _A_ and _B_ are, however, different. Author _A_ authored _n_~A~ = _a_~A~ = _a_ papers. On the other hand, author _B_ (co)authored  _n_~B~ = _a_~B~ + (_G_~B~ - 1) _a_~B~) = _G_~B~ _a_~B~ = _Ga_ papers. In the case of author _B_ the term (_G_~B~ - 1)_a_~B~ represents the papers on which author _B_ was invited. It is easy to see that as far as _G_ > 1, author _B_ will have more paper than author _A_, i.e _n_~B~ > _n_~A~. 

A natural way to correct for this bias is to taking into account the number of authors each paper has and instead of counting the papers themselves sum the inverse of the number of authors:

$$
w = \sum_{i=1}^{n} \frac{1}{1+C}.
$$

Here, 1 in the denominator symbolises the focal author, while _C_ is the number of coauthors. For author _A_, _C_ = _c_~A~ = _c_. On the other hand, for author _B_, _C_ = (_G_~B~ - 1) + _c_~B~ = (_G_ - 1) + _c_. If _c_ = 0, then the division by the number of coauthors works, we regain the number of papers the authors produced without inviting their group members.

For author _A_:
$$
w_{A} = \sum_{i=1}^{n_A} \frac{1}{1} = \sum_{i=1}^{n_A} 1 = n_A = a.
$$

For author _B_:
$$
w_{B} = \sum_{i=1}^{n_B} \frac{1}{1+G-1} = \sum_{i=1}^{Ga} \frac{1}{G} = \frac{Ga}{G} = a.
$$

On the other hand, if our focal authors collaborate with others outside of their groups the situation changes:

For author _A_:
$$
w_{A} = \sum_{i=1}^{a}\frac{1}{1+c} = \frac{a}{1+c}.
$$

For author _B_:
$$
w_{B} = \sum_{i=1}^{Ga}\frac{1}{G+c} = \frac{Ga}{G+c}.
$$

The weighted number of papers produced by author _B_ relative to author _A_, _w_~B~/_w_~A~, is:

$$
\frac{w_B}{w_A} = \frac{\frac{Ga}{G+c}}{\frac{a}{1+c}} = \frac{Ga}{G+c} \times \frac{1+c}{a} = \frac{G(1+c)}{G+c} = \frac{G+Gc}{G+c}.
$$

The proportion of _w_~B~/_w_~A~ is greater than one if _G_+_Gc_ > _G_+_c_, which is always true if _c_ > 0. This means that if authors collaborate anyone from outside their groups then authors in group _B_ will always have higher publication performance than authors in group _A_, despite the fact that the two groups have the same productivity. To compensate for this productivity bias, author _A_ should produce _w_~B~/_w_~A~ times more papers, _a_~A~/_a_~B~ = (_G_ + _Gc_)/(_G_ + _c_). This surplus of papers needed for compensating the productivity bias increases with _c_ and it keeps to _G_.

Authors in group _A_ can also compensate for the productivity bias by decreasing their collaborators from outside of the group. This reduction must be by a factor of _G_: _c_~A~ = _c_~B~ / _G_.
