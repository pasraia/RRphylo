# RRphylo 1.1.0

What’s new in version 1.1.0

1. In the new version of RRphylo the regularization factor lambda is fitted via maximum
likelihood within limits (from 0 to 10). This makes the function appreciably faster
than before
2. The function setBM still produced desired phenotypic trends (in either mean or
variance). However, in the new version the user il left to specify the intensity of the
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
