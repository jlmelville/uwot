# Mixed Data Types

The default approach of UMAP is that all your data is numeric and will
be treated as one block using the Euclidean distance metric. To use a
different metric, set the `metric` parameter, e.g. `metric = "cosine"`.

Treating the data as one block may not always be appropriate. `uwot` now
supports a highly experimental approach to mixed data types. It is not
based on any deep understanding of topology and sets, so consider it
subject to change, breakage or completely disappearing.

To use different metrics for different parts of a data frame, pass a
list to the `metric` parameter. The name of each item is the metric to
use and the value is a vector containing the names of the columns (or
their integer id, but I strongly recommend names) to apply that metric
to, e.g.:

``` r
metric = list("euclidean" = c("A1", "A2"), "cosine" = c("B1", "B2", "B3"))
```

this will treat columns `A1` and `A2` as one block of data, and generate
neighbor data using the Euclidean distance, while a different set of
neighbors will be generated with columns `B1`, `B2` and `B3`, using the
cosine distance. This will create two different simplicial sets. The
final set used for optimization is the intersection of these two sets.
This is exactly the same process that is used when carrying out
supervised UMAP (except the contribution is always equal between the two
sets and can’t be controlled by the user).

You can repeat the same metric multiple times. For example, to treat the
petal and sepal data separately in the `iris` dataset, but to use
Euclidean distances for both, use:

``` r
metric = list("euclidean" = c("Petal.Width", "Petal.Length"),
              "euclidean" = c("Sepal.Width", "Sepal.Length"))
```

## Indexing

As the `iris` example shows, using column names can be very verbose.
Integer indexing is supported, so the equivalent of the above using
integer indexing into the columns of `iris` is:

``` r
metric = list("euclidean" = 3:4, "euclidean" = 1:2)
```

but internally, `uwot` strips out the non-numeric columns from the data,
and if you use Z-scaling (i.e. specify `scale = "Z"`), zero variance
columns will also be removed. This is very likely to change the index of
the columns. If you really want to use numeric column indexes, I
strongly advise not using the `scale` argument and re-arranging your
data frame if necessary so that all non-numeric columns come after the
numeric columns.

## Categorical columns

supervised UMAP allows for a factor column to be used. You may now also
specify factor columns in the `X` data. Use the special `metric` name
`"categorical"`. For example, to use the `Species` factor in standard
UMAP for `iris` along with the usual four numeric columns, use:

``` r
metric = list("euclidean" = 1:4, "categorical" = "Species")
```

Factor columns are treated differently from numeric columns:

- They are always treated separately, one column at a time. If you have
  two factor columns, `cat1`, and `cat2`, and you would like them
  included in UMAP, you should write:

``` r
metric = list("categorical" = "cat1", "categorical" = "cat2", ...)
```

As a convenience, you can also write:

``` r
metric = list("categorical" = c("cat1", "cat2"), ...)
```

but that doesn’t combine `cat1` and `cat2` into one block, just saves
some typing.

- Because of the way categorical data is intersected into a simplicial
  set, you cannot have an X `metric` that specifies only `categorical`
  entries. You *must* specify at least one of the standard Annoy metrics
  for numeric data. For `iris`, the following is an error:

``` r
# wrong and bad
metric = list("categorical" = "Species")
```

Specifying some numeric columns is required:

``` r
# OK
metric = list("categorical" = "Species", "euclidean" = 1:4)
```

- Factor columns not explicitly included in the `metric` are still
  removed as usual.
- Categorical data does not appear in the model returned when
  `ret_model = TRUE` and so does not affect the project of data used in
  `umap_transform`. You can still use the UMAP model to project new
  data, but factor columns in the new data are ignored (effectively
  working like supervised UMAP).

## Overriding global options

Some global parameters can be overridden for a specific data block by
providing a list as the value for the metric, containing the vector of
columns as the only unnamed element, and then the over-riding keyword
arguments. An example:

``` r
  umap(
    X,
    pca = 40,
    pca_center = TRUE,
    metric = list(
      euclidean = 1:200,
      euclidean = list(201:300, pca = NULL),
      manhattan = list(300:500, pca_center = FALSE)
    )
  )
```

In this case, the first `euclidean` block with be reduced to 40
dimensions by PCA with centering applied. The second `euclidean` block
will not have PCA applied to it. The `manhattan` block will have PCA
applied to it, but no centering is carried out.

Currently, only `pca` and `pca_center` are supported for overriding by
this method, because this feature exists only to allow for the case
where you have mixed real-valued and binary data, and you want to carry
out PCA on both. It’s typical to carry out centering on real-value data
before PCA, but *not* to do so with binary data.

## `y` data

The handling of `y` data has been extended to allow for data frames, and
`target_metric` works like `metric`: multiple numeric blocks with
different metrics can be specified, and categorical data can be
specified with `categorical`. However, unlike `X`, the default behavior
for `y` is to include all factor columns. Any numeric data found will be
treated as one block, so if you have multiple numeric columns that you
want treated separately, you should specify each column separately:

``` r
target_metric = list("euclidean" = 1, "euclidean" = 2, ...)
```

I suspect that the vast majority of `y` data is one column, so the
default behavior will be fine most of the time.
