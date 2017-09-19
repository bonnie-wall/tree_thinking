# This tutorial is adapted from the following:
# http://schmitzlab.info/asr.html
# by Lars Schmitz

# Load required packages
library(ape)
library(geiger)
library(phytools)

# After loading libraries, let’s get our data.
# getting tree
haemulidTreess <- read.nexus(url("http://schmitzlab.info/Trees4dryad.nex"))
haemulidTrees
# getting data
haemulidData <- read.csv(url("http://schmitzlab.info/haemulids.csv"), header=TRUE)
head(haemulidData)

# Let's begin by plotting the first tree in the list
plot(haemulidTrees[[1]], cex=0.4)
add.scale.bar()
# OK, that can and should look much cleaner...

# Let's ladderize all trees in the distribution with lapply() and plot the
# first tree again
haemulidTreesLadderized <- lapply(haemulidTrees, ladderize)
plot(haemulidTreesLadderized[[1]], cex=0.4) #Plot the first ladderized tree
add.scale.bar()

# Looks better, doesn't it?

# However, does the order of tips in the plot match the order of tiplabels
# in the tree object?

# Please compare
haemulidTreesLadderized[[1]]$tip.label

# As you will see, by ladderizing the tree, the plotting order and the
# tiplabel order are no longer aligned. That’s a problem when plotting
# phylogenetic and phenotypic data side by side!
  
# There is an easy fix for that, though. First we will compare the tree
# to the data (a necessary step we learned about previously),
# then save the ladderized tree and finally re-load it.

# Grab the tree from above
base.tree <- haemulidTreesLadderized[[1]]

# assign row names to the dataframe
rownames(haemulidData) <- haemulidData$Taxon

# Compare tree with data
compare <- treedata(base.tree, haemulidData, sort=TRUE) 
# R will inform you of any inconsistencies and how they were solved!

# Extract vetted data and tree
haemulid.data <- as.data.frame(compare$data)
haemulid.data[,2:4] <- lapply(haemulid.data[,2:4], as.character)
haemulid.data[,2:4] <- lapply(haemulid.data[,2:4], as.numeric)

# Geiger has a function to check that the names match
# between the tree and the data frame.
name.check(base.tree, haemulid.data)

tree <- compare$phy

# Save the new tree and bring it back
write.tree(tree, "haemulid.tre")
tree <- read.tree("haemulid.tre")

# Retrieve info on habitat and adding names
z <- as.factor(haemulid.data$Habitat)
names(z) <- tree$tip.label

# Create a vector to be filled with colors for plotting purposes,
# matching the diel activity pattern
mycol<-character(length(z))
mycol[haemulid.data$Habitat=="reef"]<-"blue"
mycol[haemulid.data$Habitat=="non-reef"]<-"red"

# make very narrow margins
par(mar=c(1,0,0,1))

# And finally we can plot with tiplabels :)
plot(tree, label.offset=1.5, cex=0.5)
tiplabels(pch=21, col="black", adj=1, bg=mycol, cex=1.3)

# Now we are ready to carry out ancestral state reconstructions
# of the habitat of grunts. Did grunts originate on reefs or away from reefs?

# calculating marginal likelihoods for ancestral states with an all-rates-equal model (ER)
fit.ER <- ace(z, tree, model="ER", type="discrete")
fit.ER

fit.ER$lik.anc

# Does that make sense to you?
# Well, something doesn’t seem quite right.
# All likelihood estimates are the same, which suggests that the
# likelihood optimization got stuck along the way. Let’s try to troubleshoot.
# I suspect the branch length may be an issue.

# divide all branches by 100
scaled.tree <- tree
scaled.tree$edge.length <- scaled.tree$edge.length/100

# try again...
# calculating marginal likelihoods for ancestral states with
# an all-rates-equal model (ER)
fit.ER <- ace(z, scaled.tree, model="ER", type="discrete")
fit.ER

fit.ER$lik.anc

# OK that looks more realistic. So, what can you say about the origin of grunts?
# On reefs or not? The likelihoods for either are pretty much the same so we
# cannot tell with any amount of certainty. Some other marginal likelihoods
# seem to much more strongly facor one habitat over the other though,
# especially for shallower nodes. Let’s take a look at that.

# make very narrow margins
par(mar=c(1,0,0,1))

# plot tree
plot(tree, label.offset=1.5, cex=0.5); add.scale.bar()

# add tiplabels
tiplabels(pch=21, col="black", adj=1, bg=mycol, cex=1.3)

# add nodelabels (pie charts representing likelihoods)
nodelabels(node=1:tree$Nnode+Ntip(tree), pie=fit.ER$lik.anc, cex=0.35, piecol=c("red", "blue"))

# Please describe the evolutionary history of habitat occupation
# in grunts in one paragraph.

# Next, please explore the ace-function in more detail.
# What rate transition models other than “ER” are available?
# Do you get different results when choosing other models?

# I would like to briefly introduce an alternate approach for
# reconstructing discrete character states — stochastic character mapping,
# oiginally introduced by Huelsenbeck et al. (2003).
# Stochastic character mapping is great as it provides many samples of more
# complete histories of discrete traits over all branches of the tree,
# not just the nodes. Evolutionary transitions will also occur along the
# branches, and sometimes there may even be more than one change on the
# same branch. One often performs stochastic character mapping over entire
# tree sets, creating thousands of samples, but let’s stick to a single tree.

sim.hab <- make.simmap(scaled.tree, z, model="ER")
cols <- setNames(c("blue","red"),c("reef", "non-reef"))
plotSimmap(sim.hab, cols, fsize=0.7, ftype="i")
# That’s just one possible history — how about we repeat this 25 times?
sim.hab25 <- make.simmap(scaled.tree, z, model="ER", nsim=25)
par(mfrow=c(5,5))
for (i in 1:25){
  plotSimmap(sim.hab25[[i]], cols, ftype="off")
}

# It’s pretty opbvious that the reconstructions of habitat evolution
# are quite different from one another, which strongly suggests we
# should repeat this process many times and summarize results over
# results obtained from each simulation.

# make 100 simulations of trait histories
sim.hab100 <- make.simmap(scaled.tree, z, model="ER", nsim=100)

# summarize results
sim.summary <- describe.simmap(sim.hab100, plot=FALSE)
sim.summary

# proceed to plotting
par(mar=c(1,0,0,1))
plot(tree, label.offset=1.5, cex=0.5); add.scale.bar()

# add tiplabels
tiplabels(pch=21, col="black", adj=1, bg=mycol, cex=1.3)

# add nodelabels (pie charts representing likelihoods)
nodelabels(node=1:tree$Nnode+Ntip(tree), pie=sim.summary$ace, cex=0.35, piecol=c("red", "blue"))

# Is this result different from the maximum likelihood reconstruction?
# Also, what information is provided in the summary of the simmap trees?

