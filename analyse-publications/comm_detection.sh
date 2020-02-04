#! /bin/sh

# run community detection algorithm on collaboration matrix

# It is based on the following paper:
# Fast unfolding of communities in large networks,
# Vincent D Blondel, Jean-Loup Guillaume, Renaud Lambiotte, Etienne Lefebvre,
# Journal of Statistical Mechanics: Theory and Experiment 2008 (10),
# P10008 (12pp)
# doi: 10.1088/1742-5468/2008/10/P10008.
# ArXiv: http://arxiv.org/abs/0803.0476

# the code was downloaded from:
# https://sourceforge.net/p/louvain/
# on 2020-02-02

MATF=$1

TXTF=`echo $MATF|sed 's/\.mat$/.txt/'`
BINF=`echo $MATF|sed 's/\.mat$/.bin/'`
WEIGHTF=`echo $MATF|sed 's/\.mat$/.weights/'`
TREEF=`echo $MATF|sed 's/\.mat$/.tree/'`

sed -n '/.*,.*,.*/p' $MATF > $TXTF
sed -i 1d $TXTF
sed -i 's/,/ /g' $TXTF

convert-louvain -i $TXTF -o $BINF -w $WEIGHTF

louvain $BINF -l -1 -q id_qual -w $WEIGHTF > $TREEF

hierarchy-louvain $TREEF
