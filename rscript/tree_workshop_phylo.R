# This tutorial is adapted from the following:
# http://lukejharmon.github.io/ilhabela/instruction/2015/07/02/introduction-phylogenies-in-R/
# by Luke Harmon
# http://schmitzlab.info/phylo.html
# by Lars Schmitz

# First, we will need a few libraries installed
# Install the following packages
# to run, remove the comment character '#'
# the shortcut for this is under "Code>Comment/Uncomment Lines"
# install.packages("ape")
# install.packages("geiger")
# install.packages("nlme")
# install.packages("phytools")

# Now we load all of the necessary libraries
library(ape)
library(geiger)
library(phytools)

# Let’s begin by simulating a tree. There are many options for doing this.
# We will read an example vertebrate tree from a Newick text string
tt="(((((((cow, pig),whale),(bat,(lemur,human))),(robin,iguana)),coelacanth),gold_fish),shark);"

# We’ll use the read.tree() function of the ape package.
# The tree is stored as an object called “vert.tree”.
vert.tree<-read.tree(text=tt)

# plot this phylogeny in three different styles
plot(vert.tree) #"phylogram" is the default 
plot(vert.tree,type="cladogram")
plot(unroot(vert.tree),type="unrooted")
plot(vert.tree,type="fan")


# So, now we wonder how the phylogenetic information is encoded.
# First let’s find out what the class of the object “vert.tree” is.
class(vert.tree)

# By typing ‘vert.tree’ in the command line one can retrieve some basic information about the object.
vert.tree

# But how is this information organized within the phylo object?
# We can find out with the str() function, which displays the structure of an R object.
str(vert.tree)


# Before you continue, let’s pause for a second to explain
# Newick trees (also called "bracket form" or “parenthetic format”)
tree<-read.tree(text = "((((A,B), C), (D,E)),F);")
plot(tree,type="cladogram",edge.width=2)
tree$edge # We see a matrix of 10 rows and 2 columns. 
# This matrix represents unique combinations of node- and
# tip numbers, defining each branch segment of the tree.
# By convention, the tips of the tree are numbered 1 through n
# for n tips; and the nodes are numbered n.
# 1 through n + m for m nodes.
# m = n - 1 for a fully bifurcating tree.
# This is just to keep track of which nodes are internal and which are leaves.

tree$tip.label
# The vector tip.label contains the labels for all the tips in the tree.
# The order of tip.label is the order of the tips numbered 1 through n in edge.
tree$Nnode
# The integer Nnode contains the number of internal nodes in the tree,
# including the root of the tree if the tree is rooted.

# If it seems difficult to imagine that this object could contain all the
# information in the tree, here is a plotted tree with the edge numbers overlain:
plot(tree,edge.width=2,label.offset=0.1,type="cladogram")
nodelabels() #add node numbers
tiplabels() #add tip numbers
# We can see that all of the relationship information among the taxa
# in the tree is containing in the starting & ending nodes for each edge.
# Edges that share a common starting node number are descended from an
# immediate common ancestor, etc.

# Swapping sisterclades, extractig and identifying clades/tips, dropping tips
plot(tree, label.offset=0.2)  #same as before
nodelabels() #add node numbers
tiplabels() #add tip numbers

# How about we swap clades (D, E) and (A, B, C)?
# Their most recent common ancestor is found at node 8.
rot.phy <- rotate(tree, node=8)
# And now let’s see what happenend:
plot(rot.phy, label.offset=0.2)
nodelabels()          
tiplabels()  
# Move back and fourth between plots

# Ok, now extract the clade descended from node 9
node9 <-extract.clade(tree,9)
plot(node9, label.offset=0.2) 
# What are the tip labels?

# It will also be very helpful to select all tips in a given clade.
# This is implemented in the geiger package; the tips() function
# finds all descendants of a node.
cladeAB <- tips(rot.phy, node=10) # node 10 defines the clade composed of (A, B)
cladeAB

# Another helpful command allows for tree pruning,
# i.e. cutting of tips or entire clades.
# For example, one can delete the entire cladeABC group:
plot(rot.phy, label.offset=0.2)   #still the same as before
pruned.phy <- drop.tip(rot.phy, cladeAB)
plot(pruned.phy, label.offset=0.2)
nodelabels()
tiplabels() #did it work?

# Or we can drop tips (one or multiple) randomly.
# To prune tips, say m=2 random tips, enter:
m <- 2
pruned.phy2 <- drop.tip(rot.phy, sample(rot.phy$tip.label)[1:m]) 
# m=1 drops a single tip, of course!
plot(pruned.phy2, label.offset=0.2)   

# It may also be useful to select all branches in a specific group of tips.
# This is implemented in the ape package;
# the which.edge() function finds all edges in a specified group.
# For example, let’s identify all branches of the cladeABC group.
cladeABC <- tips(rot.phy, node=9) # node 9 defines the clade composed of (A, B, C)
cladeABCbranches <- which.edge(rot.phy, cladeABC)
cladeABCbranches #this should be a numerical vector containing 6, 7, 8, 9
# And as we can see, rows 6-9 of the edge.lengthedge.length matrix
# represent the expected branches.

# Let’s first plot the tree, and then look at the edgeedge matrix for cross-checking.
plot(rot.phy, label.offset=0.2) #we are sticking to the same plot settings
nodelabels()
tiplabels()
rot.phy$edge ## returns matrix

# Here would be a good opportunity to show you how to assign different branch colors.
# For example, how can one emphasize the branches of the clade formed by A, B, and C?
# We first create a color vector, only consisting of grey colors.
# Then we’ll assign black to all branches of clade ABC.
clcolr <- rep("darkgrey", dim(rot.phy$edge)[1]) #defining a color vector; look up rep() and dim() functions!
clcolr[cladeABCbranches] <- "red" #specifying color of branches in clade (A, B, C)
plot(rot.phy, edge.color=clcolr, edge.width=2) #also making the branches a little thicker

