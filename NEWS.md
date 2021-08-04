# RRphylo 2.5.8
What's new in version 2.5.8

New function added!

The auxilliary function node.paths collates nodes along individual lineages from the youngest (i.e. furthest from the tree root) to the oldest.

search.trend received some makeover. It now returns phenotypic values, absolute rate values, and rescaled rate values, all collated into a 'trend.data' object of class 'RRphyloList'. It stores a plot of both rescaled absolute rates and unscaled absolute rates versus time according to the argument filename.

fix.poly includes two new arguments: at resolving polytomies, tol lets specify the tolerance to consider a branch length significantly greater than zero, and random indicates whether polytomies are to be resolved randomly.

# RRphylo 2.5.7
What's new in version 2.5.7

overfitRR is implemented to test for the robustness of PGLS_fossil results. The output of overfitRR is now stored as a handy kind of list of class 'RRphyloList'.

# RRphylo 2.5.4
What's new in version 2.5.4

New function added!

The function rate.map selects RW(PC) axes linked to the highest (and lowest) evolutionary rate values and maps where and how the phenotype changed the most between any pair of taxa.

# RRphylo 2.5.0
What's new in version 2.5.0

New functions added!

The function tree.merger allows to merge phylogenetic information from different sources to create a single supertree.
The function fix.poly either collapses individual clades under a polytomy or resolves polytomous clades to non-zero length branches, dichotomous clades.
Check the "Tree-Manipulation"" vignette for details.

# RRphylo 2.4.11
What's new in version 2.4.11

In search.trend it is possible to factor out the effect of a predictor (x1.residuals argument) on the response variable in the evaluation of trends in phenotypes.

The vignette now includes the illustration of fuction cutPhylo application.

# RRphylo 2.4.9
What's new in version 2.4.9

New functions added!

The function random.evolvability.test assesses the importance of phylogenetic 
structuring (signal) on Respondability Evolvability, and Flexibility, by means of randomization.

The function conv.map selects RW(PC) axes which best account for convergence and maps convergent areas on the corresponding 3D surfaces.

RRphylo is now able to deal with categorical variables and incorporate more than one predictor (x1 argument) at the same time.

RRphylo now deals consistently with both newick and nexus tree format files.

# RRphylo 2.4.6
What's new in version 2.4.6

html-vignettes are now available for RRphylo main functions!

# RRphylo 2.4.3
What's new in version 2.4.2

New function added!

The function cutPhylo cuts all the branches of the phylogeny which are younger than a specific age or node.

# RRphylo 2.4.0
What's new in version 2.4.0

New function added!

The function phyloclust tests for phylogenetic clustering for the distribution of discrete "states" among tips and resolves clustering by removing tips randomly.

New versions of the functions search.conv, overfitRR and swapONE. search.conv now accounts for phylogenetic clustering; overfitRR now works with multivariate data and tests the robustness of search.conv results; swapONE allows the user to keep specific clades monophyletic upon indication.

# RRphylo 2.3.0
What's new in version 2.3.0

New function added!

The function scaleTree rescales a phylogenetic tree according to nodes and tips calibration ages.
The function PGLS_fossil has been implemented to work with multivariate data and to use different types of correlation structure.

# RRphylo 2.2.0
What's new in version 2.2.0

The multiple version of RRphylo has been implemented!

# RRphylo 2.1.0
What's new in version 2.1.0

New functions added!

The function StableTraitsR runs the StableTraits software (Elliot and Mooers 2014) from within the R environment.
The function overfitRR provide a test of robustness to sampling effects and phylogenetic uncertainty for both search.trend and search.shift functions.

# RRphylo 2.0.6
What's new in version 2.0.6

Improvements in search.shift and search.trend outputs and plots.

# RRphylo 2.0.0
What’s new in version 2.0.0

The package namesake function RRphylo now allows to specify values at internal nodes as derived from the fossil record.

# RRphylo 1.6.0

What’s new in version 1.6.0

A new function and new versions of "search.conv" and "distNodes" added! 

Now "search.conv" further scans the morphological convergence between species evolving under specific states 
The new version of "distNodes"  now computes distances between pairs of nodes, pairs of tips, or between nodes and tips. 
The new "PGLS_fossil" function performs pgls for non-ultrametric trees with lambda correlation structure.


# RRphylo 1.5.0

What’s new in version 1.5.0

We added a new tool. The "search.conv" function tests for morphological convergence between pairs of distant nodes.


# RRphylo 1.4.0

What’s new in version 1.4.0

New function added! The "distNodes" function returns the distance between pairs of nodes. The distance is meant as both patristic distance and the number of nodes intervening between the pair.


# RRphylo 1.3.0

What’s new in version 1.3.0

Faster and more efficient versions of setBM and search.trend.

# RRphylo 1.2.0

What’s new in version 1.2.0

RRphylo function now can be informed about the phenotype at the root. 
Users might have a fossil specimen as a reference.Otherwise, if no phenotypic 
value for the root is specified, RRphylo reads the phenotypes of one tenth of 
the tip values (those closest to the root), and calculates a weighted phenotypic 
average depending on both the tip values and the distance from the root (closest 
tips count more). Either the user-specified value or the weighted average 
phenotype will be used to constrain the root state.

# RRphylo 1.1.0

What’s new in version 1.1.0

1. In the new version of RRphylo the regularization factor lambda is fitted via maximum
likelihood within limits (from 0 to 10). This makes the function appreciably faster
than before
2. The function setBM still produced desired phenotypic trends (in either mean or
variance). However, in the new version the user is left to specify the intensity of the
simulated pattern. For the “trend” type, a scalar es defines the exponential
relationship between the trait variance and age (distance from the root). Negative es
values simulate an exponential decrease in variance. The opposite is true for es&gt;0
For the “drift” type, the scalar ds define the chance per unit time in phenotypic
means. With negative ds a pattern of decreasing mean is simulated, the other way
around for ds &gt; 0.
3. The brand new function search.trend scans for temporal patterns of trait change on
the tree (and for part of it). Provided with a RRphylo-fitted object, search.trend
computes the regression of evolutionary rates versus age, and phenotype versus
age, and compares regression slopes to randomly-generated slopes and against the
null that slope = 0. Subtrees can further be compared to each other by means of
standard major axis regression.
