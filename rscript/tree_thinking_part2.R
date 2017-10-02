# This tutorial is adapted from the following:
# http://lukejharmon.github.io/ilhabela/
# by Luke Harmon
# http://schmitzlab.info/phylo.html
# by Lars Schmitz

# Install libraries
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

## Part II - Reconstructing ancestral character states

# Load and exlpore data
# First, read in the tree and trait data
gruntTree<-read.tree(url("https://bonnie-wall.github.io/tree_thinking/data/grunts.phy"))
gruntData<-read.csv(url("https://bonnie-wall.github.io/tree_thinking/data/grunts-edit.csv"), row.names=1)
# Let’s see what this tree looks like
par(mar=c(0,0,0,0)) # sets margins
plot(gruntTree, cex=0.7) # cex adjusts text size

# Lets look at our data
head(gruntData)
# Geiger has a function to check that the names match between the tree and the data frame
name.check(gruntTree, gruntData)

# Retrieve info on habitat and add names
habitat<-gruntData[,2]
names(habitat)<-rownames(gruntData)

### Visualize the data
# assign colors to habitat names
cols<-setNames(palette()[1:length(unique(habitat))],sort(unique(habitat)))
# plot tree
plotTree(gruntTree, type="fan", fsize=0.7,ftype="i"); add.scale.bar()
# add tip labels
tiplabels(pie=to.matrix(habitat,sort(unique(habitat))),piecol=cols,cex=0.3)
# add legend
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(gruntTree)),fsize=0.8)

### Acestral state reconstructions of the habitat of grunts
# Did grunts originate on reefs or away from reefs?

# Next, let's fit a single-rate model & reconstruct ancestral states at internal nodes in the tree.
# Calculate marginal likelihoods for ancestral states with an all-rates-equal model (ER)
fitER<-ace(habitat,gruntTree,model="ER",type="discrete")
fitER

round(fitER$lik.anc,3)
# The element lik.anc gives us the marginal ancestral states,
# also known as the 'empirical Bayesian posterior probabilities'

# It is fairly straightforward to overlay these posterior probabilities on the tree:
# plot tree
plotTree(gruntTree, type="fan", fsize=0.7,ftype="i"); add.scale.bar()
# add tip labels
tiplabels(pie=to.matrix(habitat,sort(unique(habitat))),piecol=cols,cex=0.3)
# add legend
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(gruntTree)),fsize=0.8)
# add nodelabels (pie charts representing likelihoods)
nodelabels(node=1:gruntTree$Nnode+Ntip(gruntTree),
           pie=fitER$lik.anc,piecol=cols,cex=0.5)
# add nodelabels (pie charts representing likelihoods)
nodelabels(node=1:gruntTree$Nnode+Ntip(gruntTree),
           pie=fitER$lik.anc,piecol=cols,cex=0.5)

### Stochastic character mapping - MCMC approach
# An alternative approach to the one outline above is to use an MCMC approach
# to sample character histories from their posterior probability distribution.
# This is called stochastic character mapping (Huelsenbeck et al. 2003).
# The model is the same but in this case we get a sample of unambiguous histories
# for our discrete character's evolution on the tree - rather than a probability
# distribution for the character at nodes.

# For instance, given the data simulated above - we can generate the stochastic character map as follows:
# simulate single stochastic character map using empirical Bayes method
mtree<-make.simmap(gruntTree,habitat,model="ER")

plotSimmap(mtree,cols,type="fan",fsize=0.8,ftype="i")

add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(gruntTree)),fsize=0.8)

## Part III - Testing for correlated evolution of continuous traits

# Is there a correlation between 'standard length' and 'raker length'?

par(mar=c(5.1,4.1,4.1,2.1)) # sets margins
plot(Raker_Length~Standard_length, data=gruntData)
# It looks like there might be...

### Phylogenetically independent contrasts (PICs)
# We can do this analysis easily with PICs

# Extract columns
Standard_length <- gruntData[, "Standard_length"]
Raker_Length <- gruntData[, "Raker_Length"]

# Give them names
names(Standard_length) <- names(Raker_Length) <- rownames(gruntData)

# Calculate phylogenetically independent contrasts (PICs)
t1Pic <- pic(Standard_length, gruntTree)
t2Pic <- pic(Raker_Length, gruntTree)

# Make a model
picModel <- lm(t1Pic ~ t2Pic - 1)

# Is it significant?
summary(picModel)

# plot results
plot(t1Pic ~ t2Pic)
abline(a = 0, b = coef(picModel))

### Phylogenetic generalized least squares (PGLS)
## This whole procedure can be carried out more simply and with more flexibility using PGLS

#### Brownian Correlation Strucure
pglsResult<-gls(Standard_length~Raker_Length,
                cor=corBrownian(1, phy=gruntTree),
                data=gruntData,
                method="ML")

#### Martins Correlation Structure
# We can also assume that the error structure follows an OU model rather than Brownian motion
pglsResultOU<-gls(Standard_length~Raker_Length,
                  cor=corMartins(1, phy=gruntTree, fixed=T),
                  data=gruntData,
                  method="ML")

#### Pagel's “lambda” Correlation Structure
pglsResultPagel<-gls(Standard_length~Raker_Length,
                     cor=corPagel(1, phy=gruntTree, fixed=T),
                     data=gruntData,
                     method="ML")

anova(pglsResult, pglsResultOU, pglsResultPagel)

# the main one is the lowest AIC, choose that

anova(pglsResultPagel)

# Is there is a relationship?

