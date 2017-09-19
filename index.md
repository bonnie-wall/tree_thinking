# Workshop 5: BIOL20001 Evolution: Making Sense of Life
## Tree thinking: A tutorial in R
### Learning objectives
-	Introduce you to using trees to address evolutionary questions and to the R statistical environment
-	Develop basic skills in tree visualisation
-	Mapping character states and reconstructing ancestral characters
-	Exploring convergent evolution
-	Testing for correlated evolution


Download the following R scripts:

A brief introduction to phylogenetic trees in R:

[tree_workshop_phylo.R](rscript/tree_workshop_phylo.R)

Reconstructing ancestral states of discrete characters:

[tree_workshop_asr.R](/rscript/tree_workshop_asr.R)



### A brief introduction to phylogenetic trees in R

Install the necessary libraries below

```
## to run, remove the comment character '#'
#install.packages("ape")
#install.packages("geiger")
#install.packages("nlme")
#install.packages("phytools")
```
Now we load all of the necessary libraries
```
library(ape)
library(geiger)
library(phytools)
```
Let’s begin by simulating a tree. There are many options for doing this. 
We will read an example vertebrate tree from a Newick text string
```
tt="(((((((cow, pig),whale),(bat,(lemur,human))),(robin,iguana)),coelacanth),gold_fish),shark);"
```

We’ll use the ```read.tree()``` function of the ape package. 
The tree is stored as an object called ```vert.tree```. 

```
vert.tree<-read.tree(text=tt)

#plot this phylogeny in three different styles
plot(vert.tree) #"phylogram" is the default 
plot(vert.tree,type="cladogram")
plot(unroot(vert.tree),type="unrooted")
plot(vert.tree,type="fan")
```

So, now we wonder how the phylogenetic information is encoded.
First let’s find out what the class of the object “vert.tree” is.
```
class(vert.tree)
```

By typing ‘vert.tree’ in the command line one can retrieve some 
basic information about the object.
```
vert.tree
```

But how is this information organized within the phylo object? 
We can find out with the str() function, which displays the structure of an R object.
```
str(vert.tree)
```

Before you continue, let’s pause for a second to explain 
Newick trees (also called "bracket form" or “parenthetic format”)
```
tree<-read.tree(text = "((((A,B), C), (D,E)),F);")
plot(tree,type="cladogram",edge.width=2)
tree$edge # We see a matrix of 10 rows and 2 columns. 
```
This matrix represents unique combinations of node- and tip numbers, defining each branch segment of the tree.By convention, the tips of the tree are numbered 1 through n for n tips; and the nodes are numbered n. 1 through n + m for m nodes. m = n - 1 for a fully bifurcating tree. This is just to keep track of which nodes are internal and which are leaves.

```
tree$tip.label
```
The vector tip.label contains the labels for all the tips in the tree.
The order of tip.label is the order of the tips numbered 1 through n in edge.
```
tree$Nnode
```
The integer Nnode contains the number of internal nodes in the tree, including the root of the tree if the tree is rooted.

If it seems difficult to imagine that this object could contain all the information in the tree, here is a plotted tree with the edge numbers overlain:
```
plot(tree,edge.width=2,label.offset=0.1,type="cladogram")
nodelabels() #add node numbers
tiplabels() #add tip numbers
```
We can see that all of the relationship information among the taxa in the tree is containing in the starting & ending nodes for each edge. Edges that share a common starting node number are descended from an immediate common ancestor, etc.

Swapping sisterclades, extractig and identifying clades/tips, dropping tips
```
plot(tree, label.offset=0.2)  #same as before
nodelabels() #add node numbers
tiplabels() #add tip numbers
```

How about we swap clades (D, E) and (A, B, C)? Their most recent common ancestor is found at node 8.
```
rot.phy <- rotate(tree, node=8)
```
And now let’s see what happenend:
```
plot(rot.phy, label.offset=0.2)
nodelabels()          
tiplabels()  
```
Move back and fourth between plots

Ok, now extract the clade descended from node 9
```
node9 <-extract.clade(tree,9)
plot(node9, label.offset=0.2) 
```
What are the tip labels?

It will also be very helpful to select all tips in a given clade. This is implemented in the geiger package; the tips() function finds all descendants of a node.
```
cladeAB <- tips(rot.phy, node=10) # node 10 defines the clade composed of (A, B)
cladeAB
```

Another helpful command allows for tree pruning, i.e. cutting of tips or entire clades. For example, one can delete the entire cladeABC group:
```
plot(rot.phy, label.offset=0.2)   #still the same as before
pruned.phy <- drop.tip(rot.phy, cladeAB)
plot(pruned.phy, label.offset=0.2)
nodelabels()
tiplabels() #did it work?
```

Or we can drop tips (one or multiple) randomly. To prune tips, say m=2 random tips, enter:
```
m <- 2
pruned.phy2 <- drop.tip(rot.phy, sample(rot.phy$tip.label)[1:m]) 
# m=1 drops a single tip, of course!
plot(pruned.phy2, label.offset=0.2) 
```

It may also be useful to select all branches in a specific group of tips. This is implemented in the ape package; the which.edge() function finds all edges in a specified group. For example, let’s identify all branches of the cladeABC group.
```
cladeABC <- tips(rot.phy, node=9) # node 9 defines the clade composed of (A, B, C)
cladeABCbranches <- which.edge(rot.phy, cladeABC)
cladeABCbranches #this should be a numerical vector containing 6, 7, 8, 9
```
And as we can see, rows 6-9 of the edge.lengthedge.length matrix represent the expected branches. 

Let’s first plot the tree, and then look at the edgeedge matrix for cross-checking.
```
plot(rot.phy, label.offset=0.2) #we are sticking to the same plot settings
nodelabels()
tiplabels()
rot.phy$edge ## returns matrix
```

Here would be a good opportunity to show you how to assign different branch colors. 
For example, how can one emphasize the branches of the clade formed by A, B, and C? 
We first create a color vector, only consisting of grey colors. 
Then we’ll assign black to all branches of clade ABC.
```
clcolr <- rep("darkgrey", dim(rot.phy$edge)[1]) #defining a color vector; look up rep() and dim() functions!
clcolr[cladeABCbranches] <- "red" #specifying color of branches in clade (A, B, C)
plot(rot.phy, edge.color=clcolr, edge.width=2) #also making the branches a little thicker
```

### Reconstructing ancestral states of discrete characters

After loading libraries, let’s get our data.
getting tree
```
haemulidTreess <- read.nexus(url("http://schmitzlab.info/Trees4dryad.nex"))
haemulidTrees
```
getting data
```
haemulidData <- read.csv(url("http://schmitzlab.info/haemulids.csv"), header=TRUE)
head(haemulidData)
```

Let's begin by plotting the first tree in the list
```
plot(haemulidTrees[[1]], cex=0.4)
add.scale.bar()
```
OK, that can and should look much cleaner...

Let's ladderize all trees in the distribution with lapply() and plot the first tree again
```
haemulidTreesLadderized <- lapply(haemulidTrees, ladderize)
plot(haemulidTreesLadderized[[1]], cex=0.4) #Plot the first ladderized tree
add.scale.bar()
```
Looks better, doesn't it?

However, does the order of tips in the plot match the order of tiplabels in the tree object?
```
# Please compare
haemulidTreesLadderized[[1]]$tip.label
```

As you will see, by ladderizing the tree, the plotting order and the tiplabel order are no longer aligned. That’s a problem when plotting phylogenetic and phenotypic data side by side!

There is an easy fix for that, though. First we will compare the tree to the data, then save the ladderized tree and finally re-load it.

```
#grab the tree from above
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
```

Now we are ready to carry out ancestral state reconstructions of the habitat of grunts. Did grunts originate on reefs or away from reefs?

Calculating marginal likelihoods for ancestral states with an all-rates-equal model (ER)
```
fit.ER <- ace(z, tree, model="ER", type="discrete")
fit.ER

fit.ER$lik.anc
```

Does that make sense to you? Well, something doesn’t seem quite right. All likelihood estimates are the same, which suggests that the likelihood optimization got stuck along the way. Let’s try to troubleshoot. I suspect the branch length may be an issue.

```
# divide all branches by 100
scaled.tree <- tree
scaled.tree$edge.length <- scaled.tree$edge.length/100
```

try again...
Calculating marginal likelihoods for ancestral states with an all-rates-equal model (ER)
```
fit.ER <- ace(z, scaled.tree, model="ER", type="discrete")
fit.ER

fit.ER$lik.anc
```

OK that looks more realistic. So, what can you say about the origin of grunts? On reefs or not? The likelihoods for either are pretty much the same so we cannot tell with any amount of certainty. Some other marginal likelihoods seem to much more strongly facor one habitat over the other though, especially for shallower nodes. Let’s take a look at that.

```
# make very narrow margins
par(mar=c(1,0,0,1))

# plot tree
plot(tree, label.offset=1.5, cex=0.5); add.scale.bar()

# add tiplabels
tiplabels(pch=21, col="black", adj=1, bg=mycol, cex=1.3)

# add nodelabels (pie charts representing likelihoods)
nodelabels(node=1:tree$Nnode+Ntip(tree), pie=fit.ER$lik.anc, cex=0.35, piecol=c("red", "blue"))
```

Please describe the evolutionary history of habitat occupation in grunts in one paragraph.

Next, please explore the ace-function in more detail. What rate transition models other than “ER” are available? Do you get different results when choosing other models?

I would like to briefly introduce an alternate approach for reconstructing discrete character states — stochastic character mapping, oiginally introduced by Huelsenbeck et al. (2003). Stochastic character mapping is great as it provides many samples of more complete histories of discrete traits over all branches of the tree, not just the nodes. Evolutionary transitions will also occur along the branches, and sometimes there may even be more than one change on the same branch. One often performs stochastic character mapping over entire tree sets, creating thousands of samples, but let’s stick to a single tree.

```
sim.hab <- make.simmap(scaled.tree, z, model="ER")
cols <- setNames(c("blue","red"),c("reef", "non-reef"))
plotSimmap(sim.hab, cols, fsize=0.7, ftype="i")
# That’s just one possible history — how about we repeat this 25 times?
sim.hab25 <- make.simmap(scaled.tree, z, model="ER", nsim=25)
par(mfrow=c(5,5))
for (i in 1:25){
  plotSimmap(sim.hab25[[i]], cols, ftype="off")
}
```

It’s pretty opbvious that the reconstructions of habitat evolution are quite different from one another, which strongly suggests we should repeat this process many times and summarize results over results obtained from each simulation.

```
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
```

Is this result different from the maximum likelihood reconstruction?
Also, what information is provided in the summary of the simmap trees?
