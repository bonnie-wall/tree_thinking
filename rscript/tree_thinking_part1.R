# This tutorial is adapted from the following:
# http://lukejharmon.github.io/ilhabela/instruction/2015/07/02/introduction-phylogenies-in-R/
# by Luke Harmon
# http://schmitzlab.info/phylo.html
# by Lars Schmitz

## Part I – Importing and visualising trees in R

# Install the following packages
# to run, remove the comment character '#'
# install.packages("ape")
# install.packages("geiger")
# install.packages("nlme")
# install.packages("phytools")

# Load necessary libraries
library(ape)
library(geiger)
library(nlme)
library(phytools)

# Let’s begin by simulating a tree. There are many options for doing this.
# We will read an example vertebrate tree from a Newick text string
# Also called "bracket form" or “parenthetic format”
tt = "(((((((cow, pig),whale),(bat,(lemur,human))),(robin,iguana)),coelacanth),gold_fish),shark);"


# We’ll use the read.tree() function of the ape package.
# The tree is stored as an object called “vert.tree”.
vert.tree <- read.tree(text=tt)

# plot this phylogeny in three different styles
plot(vert.tree) #"phylogram" is default 
plot(vert.tree,type="cladogram")
plot(unroot(vert.tree),type="unrooted")
plot(vert.tree,type="fan")

# Before you continue, let’s pause for a second to explain Newick trees
tt2 <- "(((((((cow:0.74, pig:0.74):0.76,whale:1.5):0.58,(bat:1.38,(lemur:0.46,human:0.46):0.92):0.70):0.75,(robin:2.09,iguana:2.09):0.74):0.63,coelacanth:3.46):0.56,gold_fish:4.02):0.98,shark:5):0.58;"
tree <- read.tree(text=tt2)
plot(tree)
edgelabels(tree$edge.length,adj=c(0.5, -0.25),cex=0.7) # cex adjusts text size

tree$tip.label
# The vector tip.label contains the labels for all the tips in the tree.
# The order of tip.label is the order of the tips numbered 1 through n in edge.
tree$Nnode
# The integer Nnode contains the number of internal nodes in the tree,
# including the root of the tree if the tree is rooted.

# If it seems difficult to imagine that this object could contain all the
# information in the tree, here is a plotted tree with the edge numbers overlain:
plot(vert.tree,label.offset=0.3)
nodelabels() #add node numbers
tiplabels() #add tip numbers
# We can see that all of the relationship information among the taxa
# in the tree is containing in the starting & ending nodes for each edge.
# Edges that share a common starting node number are descended from an
# immediate common ancestor, etc.