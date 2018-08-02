---
title: "UMAP for t-SNE"
output:
  html_document:
    theme: cosmo
    toc: true
    toc_float:
      collapsed: false
---

This is part of the documentation for [UWOT](https://github.com/jlmelville/uwot).

The [UMAP paper](https://arxiv.org/abs/1802.03426) does not go into much 
implementation detail. If you are coming to UMAP from t-SNE and don't know
much about topology or fuzzy sets (and I certainly don't), you may find yourself
hankering for some insight into how UMAP achieves its results and what
connection there is with t-SNE and related methods.

Here are some details I have taken
from scouring the [Python source code](https://github.com/lmcinnes/umap) and 
from asking UMAP creator [Leland McInnes](https://github.com/lmcinnes). In what 
follows, I assume that you are already familiar with how t-SNE works.

Broadly, the UMAP implementation uses a similar approach to
[LargeVis](https://arxiv.org/abs/1602.00370). LargeVis in turn uses
concepts from [t-SNE](https://lvdmaaten.github.io/tsne/). t-SNE itself is a
modification of the original 
[Stochastic Neighbor Embedding (SNE)](https://papers.nips.cc/paper/2276-stochastic-neighbor-embedding)
method.

Different papers use different symbols and nomenclature, so to make sure we're
all on the same page, I will restate the relevant definitions from SNE and
t-SNE before we take a look at LargeVis and UMAP.

## (Asymmetric) SNE

Given $N$ observations of some high dimensional data, for any pair,
$\mathbf{x_{i}}$ and $\mathbf{x_{j}}$, SNE defines the similarity (aka an
affinity or weight) between them, $v_{j|i}$, using a Gaussian kernel function:

$$
v_{j|i} = \exp(-\beta_{i} r_{ij}^2)
$$

where $r_{ij}$ is the distance between $\mathbf{x_{i}}$ and $\mathbf{x_{j}}$ and
$\beta_{i}$ must be determined by some method (we'll get back to that). The
notation of $v_{j|i}$ rather than $v_{ij}$, is to indicate that this quantity is
not symmetric, i.e. $v_{j|i} \neq v_{i|j}$. I've borrowed this notation from the
conditional versus joint probability definitions used in symmetric SNE (see
below) but we'll also need it for quantities other than probabilities. The
$r_{ij}$ notation indicates that the distances are symmetric and this convention
will be used for other symmetric values.

The weights are normalized to form $N$ probability distributions:

$$
p_{j|i} = \frac{v_{j|i}}{\sum_{k}^{N} v_{k|i}} 
$$
$\beta_{i}$ is chosen by finding that value that results in the probability
distribution having a specific perplexity. The perplexity has to be chosen by
the user, but is interpreted as being a continuous version of the number of
nearest neighbors, and generally is chosen to take values between 5 and 50.

$p_{j|i}$ is a conditional probability, and is interpreted as meaning "the
probability that you would pick item $j$ as being similar to item $i$, given
that you've already picked $i$".

In the output space of the embedded coordinates, the similarity between the
points $\mathbf{y_i}$ and $\mathbf{y_j}$ is also defined as a Gaussian:

$$
w_{ij} = \exp(-d_{ij}^2)
$$

where $d_{ij}$ is the Euclidean distance between $\mathbf{y_i}$ and
$\mathbf{y_j}$. There is no $\beta$ in this weight definition so these weights 
are symmetric. The output probabilities, $q_{j|i}$ are calculated from $w_{ij}$
in the same way that we go from $v_{j|i}$ to $p_{j|i}$, again creating $N$
probability distributions. Due to normalizing by rows, the $q_{j|i}$ are 
asymmetric despite the symmetric weights they are generated from.

The SNE cost function is the sum of the Kullback-Leibler divergences of the $N$
distributions:

$$
C_{SNE} = \sum_{i}^{N} \sum_{j}^{N} p_{j|i} \log \frac{p_{j|i}}{q_{j|i}} 
$$
In all of the above (and in what follows), weights and probabilities when $i =
j$ are not defined. I don't want to clutter the notation further, so assume they
are excluded from any sums.

## Symmetric SNE

In Symmetric SNE, the input probability matrix is symmetrized by averaging 
$p_{j|i}$ and $p_{i|j}$ and then re-normalized over all pairs of points, 
to create a single (joint) probability distribution, $p_{ij}$:

$$
p_{ij} = \frac{p_{j|i} + p_{i|j}}{2N}
$$

The output probabilities, $q_{ij}$ are now defined by normalizing the output
weights over all pairs, again creating a single probability distribution:

$$
q_{ij} = \frac{w_{ij}}{\sum_{k}^N \sum_{l}^N w_{kl}} 
$$

The cost function for SSNE is then:

$$
C_{SSNE} = \sum_{i}^{N} \sum_{j}^{N} p_{ij} \log \frac{p_{ij}}{q_{ij}} 
$$

## t-SNE

For the purposes of this discussion, t-SNE only differs from symmetric SNE by
its weight function:

$$
w_{ij} = \frac{1}{1 + d_{ij}^2}
$$

## SNE Optimization

Optimization in t-SNE proceeds by:

* Calibrate the input probabilities, $p_{ij}$ according to the desired
perplexity (this only needs to be done once). 
* Iteratively: 
* Calculate all pairwise distances, $d_{ij}$ from $\mathbf{y_{i}}$ and 
$\mathbf{y_{j}}$ 
* Calculate the weights, $w_{ij}$ 
* Calculate the output probabilities $q_{ij}$
* Use gradient descent to update the $\mathbf{y_{i}}$ and $\mathbf{y_{j}}$

This is fundamentally $O(N^2)$ because of the need to calculate pairwise
distances: the t-SNE gradient requires the $q_{ij}$ to be calculated and the
normalization step that converts $w_{ij}$ to $q_{ij}$ requires 
all of $w_{ij}$ to be calculated so you also need all the distances.

Approaches like [Barnes-Hut t-SNE](https://arxiv.org/abs/1301.3342) and others 
(e.g. [Flt-SNE](https://arxiv.org/abs/1712.09005)) attempt to improve on this by
taking advantage of the t-SNE gradient:

* The attractive part of the gradient depends on $p_{ij}$, which is constant and
only large for neighbors that are close in the input space. Therefore it's only
necessary to calculate the attractive gradient for nearest neighbors of 
$\mathbf{x_i}$. In Barnes-Hut t-SNE, the number of nearest neighbors used is
three times whatever the perplexity is. For larger datasets, a perplexity of
50 is common, so usually you are looking for the 150-nearest neighbors of each point.
* The repulsive part of the gradient is dependent on $q_{ij}$ which changes
with each iteration, so the improvements here focus on grouping together points
which are distant in the output space and treating them as a single point for
the purposes of the gradient calculation.

As you will be able to tell from perusing the publications linked to above,
these approaches are increasing in sophistication and complexity.

## LargeVis

LargeVis takes a different approach: it re-uses a lot of the same definitions as
t-SNE, but makes sufficient modifications so that it's possible to use
stochastic gradient descent.

Also, rather than talk about probabilities, LargeVis uses the language of graph
theory. Each observation in our dataset is now considered to be a vertex or node
and the similarity between them is the weight of the edge between the two
vertices. Conceptually we're still talking about elements in a matrix, but I
will start slipping into the language of "edges" and "vertices".

The key change is the cost function, which is now a maximum likelihood function:

$$
L_{LV} = \sum_{ \left(i, j\right) \in E} p_{ij} \log w_{ij} 
+\gamma \sum_{\left(i, j\right) \in \bar{E}} \log \left(1 - w_{ij} \right)
$$
$p_{ij}$ and $w_{ij}$ is the same as in t-SNE (the authors try some alternative
$w_{ij}$ definitions, but they aren't as effective).

The new concepts here are $\gamma$ and $E$. $\gamma$ is a user-defined positive
scalar to weight repulsive versus attractive forces. Its default in the 
[reference implementation](https://github.com/lferry007/LargeVis) is 7.

$E$ is the set of edges with a non-zero weight. This is the graph theory way to
talk about nearest neighbors in the input space. Just as with Barnes-Hut t-SNE,
we find a set of nearest neighbors for each point $\mathbf{x_i}$ and only 
define input weights and probabilities for pairs of points which are nearest
neighbors. As with the official Barnes-Hut t-SNE implementation, the LargeVis 
reference implementation uses a default perplexity of 50, and the default number of 
nearest neighbors is 3 times the perplexity.

This cost function therefore consists of two disjoint contributions: nearest
neighbors in the input space contribute to the attractive part of the cost
function (the first part). Everything else contributes to the second, repulsive 
part.

The key advantage of this cost function over the KL divergence is that it
doesn't contain $q_{ij}$. With no output normalization, we don't need to
calculate all the output pairwise distances. So this cost function is amenable
to stochastic gradient descent techniques.

## The LargeVis sampling strategy

To calculate a stochastic gradient, LargeVis does the following:

* Samples an edge in $E$, i.e. chooses $i$ and $j$ such that $p_{ij} \neq 0$. 
This is called a "positive edge" in the LargeVis paper. $i$ and $j$ are used to 
calculate the attractive part of the gradient.
* Samples one of the $N$ vertices, let's call it $k$. As datasets grow larger,
the probability that $k$ is in $E$ grows smaller, but that's not actually
checked for (only that $k \neq i$). This is used to calculate the repulsive
part of the gradient. These are called the "negative samples" in the LargeVis
paper.
* Repeat the negative sample step a number of times. The default in LargeVis 
is to sample 5 negatives for each positive edge.

The coordinates of $i$, $j$ and the various $k$ are then updated according to
the gradients. This concludes one iteration of the SGD.

The attractive and repulsive gradients for LargeVis are respectively:

$$
\frac{\partial L_{LV}}{\partial \mathbf{y_i}}^+ =
\frac{-2}{1 + d_{ij}^2}p_{ij} \left(\mathbf{y_i - y_j}\right) \\
\frac{\partial L_{LV}}{\partial \mathbf{y_i}}^- =
\frac{2\gamma}{\left(0.1 + d_{ij}^2\right)\left(1 + d_{ij}^2\right)} \left(\mathbf{y_i - y_j}\right)
$$
The value of 0.1 that appears in the repulsive gradient is there to prevent
division by zero.

Sampling of edges and vertices is not uniform. For the attractive gradient, the
authors note that the factor of $p_{ij}$ that appears means that the magnitude
of the gradient can differ hugely between samples to the extent that choosing an
appropriate learning rate can be difficult. Instead they sample the edges
proportionally to $p_{ij}$ and then for the gradient calculation, treat each
edge as if the weights were all equal. The attractive gradient used as part of
LargeVis SGD is therefore:

$$
\frac{\partial L_{LV}}{\partial \mathbf{y_i}}^+ =
\frac{-2}{1 + d_{ij}^2} \left(\mathbf{y_i - y_j}\right)
$$

As $p_{ij}$ doesn't appear in the repulsive part of the gradient, so it would
seem that uniform sampling would work for the negative sampling. However,
vertices are sampled using a "noisy" distribution proportional to their degree ^
0.75, where the degree of the vertex is the sum of the weights of the edges
incident to them. There doesn't seem to be a theoretical reason to use the
degree ^ 0.75. It's based on results from the field of word embeddings: the
LargeVis authors reference a
[skip-gram](http://papers.nips.cc/paper/5021-distributed-representations-of-words-andphrases)
paper, but the same power also shows up in 
[GloVE](https://nlp.stanford.edu/projects/glove/). In both cases it is justified
purely empirically. The `uwot` version of LargeVis (`lvish`) samples the 
negative edges uniformly, and it doesn't seem to cause any problems.

## UMAP (at last)

The UMAP cost function is the cross-entropy of two fuzzy sets, which
can be represented as symmetric weight matrices:

$$
C_{UMAP} = 
\sum_{ij} \left[ v_{ij} \log \left( \frac{v_{ij}}{w_{ij}} \right) + 
(1 - v_{ij}) \log \left( \frac{1 - v_{ij}}{1 - w_{ij}} \right) \right]
$$
$v_{ij}$ are symmetrized input affinities, and are not probabilities. The graph
interpretation of them as weights of edges in a graph still applies, though. 
These are arrived at differently to t-SNE and LargeVis. The unsymmetrized UMAP 
input weights are given by:

$$
v_{j|i} = \exp \left[ -\left( r_{ij} - \rho_{i} \right) / \sigma_{i} \right]
$$

where $r_{ij}$ are the input distances, $\rho_{i}$ is the distance to the
nearest neighbor (ignoring zero distances where neighbors are duplicates) and
$\sigma_{i}$ is analogous to $\beta_{i}$ in the perplexity calibration used in
SNE. In this case, $\sigma_{i}$ is determined such that 
$\sum_{j} v_{j|i} = \log_{2} k$ where $k$ is the number of nearest neighbors.

These weights are symmetrized by a slightly different method to SNE: 

$$
v_{ij} = \left(v_{j|i} + v_{i|j}\right) - v_{j|i}v_{i|j}
$$
or as a matrix operation:

$$
V_{symm} = V + V^{T} - V \circ V^{T}
$$

where $T$ indicates the transpose and $\circ$ is the Hadamard (i.e. entry-wise) 
product. This effectively carries out a fuzzy set union.

The output weights are given by:

$$
w_{ij} = 1 / \left(1 + ad_{ij}^{2b}\right)
$$

where $a$ and $b$ are determined by a non-linear least squares fit based on the
`min_dist` and `spread` parameters that control the tightness of the squashing
function. By setting $a = 1$ and $b = 1$ you get t-SNE weighting back.
The current UMAP defaults result in $a = 1.929$ and $b = 0.7915$.

The attractive and repulsive UMAP gradient expressions are, respectively:

$$
\frac{\partial C_{UMAP}}{\partial \mathbf{y_i}}^+ = 
\frac{-2abd_{ij}^{2\left(b - 1\right)}}{1 + d_{ij}^2}  v_{ij} \left(\mathbf{y_i - y_j}\right) \\
\frac{\partial C_{UMAP}}{\partial \mathbf{y_i}}^- = 
\frac{b}{\left(0.001 + d_{ij}^2\right)\left(1 + d_{ij}^2\right)}\left(1 - v_{ij}\right)\left(\mathbf{y_i - y_j}\right)
$$
While more complex-looking than the LargeVis gradient, there are obvious
similarities. The 0.001 term in the denominator of the repulsive gradient plays
the same role as the 0.1 in the LargeVis gradient (preventing division by zero).

UMAP uses the same sampling strategy as LargeVis, where sampling of positive
edges is proportional to the weight of the edge (in this case $v_{ij}$), and
then the value of the gradient is calculated by assuming that $v_{ij} = 1$
for all edges. So for SGD purposes, the attractive gradient for UMAP is:

$$
\frac{\partial C_{UMAP}}{\partial \mathbf{y_i}}^+ = 
\frac{-2abd_{ij}^{2\left(b - 1\right)}}{1 + d_{ij}^2}\left(\mathbf{y_i - y_j}\right)
$$

The repulsive part of the gradient contains a $1 - v_{ij}$ term, but because
$v_{ij} = 0$ for most pairs of edges, that term effectively disappears, leaving:

$$
\frac{\partial C_{UMAP}}{\partial \mathbf{y_i}}^- = 
\frac{b}{\left(0.001 + d_{ij}^2\right)\left(1 + d_{ij}^2\right)}
  \left(\mathbf{y_i - y_j}\right)
$$

Unlike LargeVis, negative sampling in UMAP uses a uniform distribution.

It's worth considering that, although LargeVis uses $p_{ij}$ and UMAP uses
$v_{ij}$ in their cost functions, the difference isn't that important, because
sampling proportionally to $p_{ij}$ is exactly the same as sampling
proportionally to $v_{ij}$. In fact, if you look at the LargeVis reference
implementation (or `lvish` in `uwot`), the input affinities are symmetrized, but
not divided by $N$. Nonetheless, because the affinities are only used to
construct the sampling probabilities, the presence of the $\gamma$ parameter in
the repulsive part of the gradient means that you are effectively using
$p_{ij}$ in the LargeVis cost function that is being optimized, not $v_{ij}$.

## Minor UMAP Variations

There are some extra parameters in UMAP that make some minor changes if modified
from their default-values:

* `local_connectivity`
affects $\rho_{i}$ by using the distance to the `local_connectivity`th
non-zero near neighbor (or by interpolating between distances if
`local_connectivity` is non-integral).

* `set_op_mix_ratio`
This changes the form of the symmetrization from fuzzy set union to fuzzy set
intersection which is just $v_{j|i}v_{i|j}$, and can also blend between the two.

* `gamma`
This works exactly like the value in LargeVis, up-weighting the repulsive 
contribution of the gradient.

*2 August 2018*: The follow parameter no longer appears in the reference UMAP implementation:

* `bandwidth` affects $v_{j|i}$ by multiplying the value of $\sigma_{i}$:

$$
v_{j|i} = \exp \left[ -\left( r_{ij} - \rho_{i} \right) / \beta \sigma_{i} \right]
$$

where `bandwidth` is represented as $\beta$. The value of $\sigma_{i}$
is determined using calculations of $v_{j|i}$ without $\beta$, before recalculating
the $v_{j|i}$ using `bandwidth`.

I'm not sure how useful any of these are when changed from the defaults.

My thanks to [Dmitry Kobak](https://github.com/dkobak) for some very helpful 
discussions and typo-spotting.
