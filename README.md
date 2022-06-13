
# PaperPump

## Introduction

Under recent research assessment protocols questionable authorship practices seem to flower. An especially problematic one these is when a group of authors form a publication cartel where each of them reciprocally gifts honorary authorship to others in the cartel.

In this manuscript I present a simple model to show that members of publication cartels are always better off than researchers not cartelling even if the 1/n rule (publications are weighted by the inverse of their number of authors) is used to assess their productivity. With a more realistic simulation, on the other hand, I found that the use of the 1/n rule can generate conflicts of interest among cartel members and between members and non-members which might lead to the self-purification of the academic publishing industry.

This repository contains everything to reproduce my paper.

## Make the paper

The code below produces the manuscript in pdf format. You should have `R` and `julia` installed. It is only tested on Linux.

```shell
make init
make
```

To reproduce my PLOS ONE submission run the following:

```shell
make plosone
```
