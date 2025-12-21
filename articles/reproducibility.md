# Reproducibility

`uwot` relies on the underlying compiler and C++ standard library on
your machine and this can result in differences in output even with the
same input data, arguments, packages and R version. If you require
reproducibility between machines, it is strongly suggested that you
stick with the same OS and compiler version on all of them (e.g. a fixed
LTS of a Linux distro and gcc version). Otherwise, the following can
help:

- Use the `tumap` method instead of `umap`. This avoid the use of
  `std::pow` in gradient calculations. This also has the advantage of
  being faster to optimize. However, this gives larger clusters in the
  output, and you don’t have the ability to control that with `a` and
  `b` (or `spread` and `min_dist`) parameters.
- For `umap`, it’s better to provide `a` and `b` directly with a fixed
  precision rather than allowing them to be calculated via the `spread`
  and `min_dist` parameters. For default UMAP, use
  `a = 1.8956, b = 0.8006`.
- Use `approx_pow = TRUE`, which avoids the use of the `std::pow`
  function.
- Use `init = "spca"` rather than `init = "spectral"` (although the
  latter is the default and preferred method for UMAP initialization).
- If `n_sgd_threads` is set larger than `1`, then even if you use
  `set.seed`, results of the embeddings are not repeatable, This is
  because there is no locking carried out on the underlying coordinate
  matrix, and work is partitioned by edge not vertex and a given vertex
  may be processed by different threads. The order in which reads and
  writes occur is of course at the whim of the thread scheduler. This is
  the same behavior as LargeVis.
- Use `rng_type = "deterministic`, which will make vertex sampling
  during the optimization deterministic. Note that this will *not*
  affect the use of a random number generator in other parts of the
  algorithm, such as approximate nearest neighbor search and
  initialization. This may give slightly less accurate results due to
  the lack of random sampling but the trade-off may be worth it (and
  it’s also a bit faster).
- For random number generation, you can provide a `seed` parameter. This
  doesn’t do anything other than call `set.seed` for you inside the
  routine, but you may find it convenient to fix the seed in the call to
  `umap`.
- If you use `nn_method = "hnsw"` or `nn_method = "nndescent"`, then by
  default these use multiple threads for index building. This can lead
  to non-determinism in the nearest neighbor calculation. Set
  `n_build_threads = 1` to avoid this. Note that the default nearest
  neighbor method using Annoy is unaffected by this parameter, because
  it always uses one thread to build the index.

In summary, your chances of reproducibility are increased by using:

``` r
mnist_umap <- umap(mnist, a = 1.8956, b = 0.8006, approx_pow = TRUE, init = "spca", batch = TRUE, rng_type = "deterministic", seed = 42)
# or
mnist_tumap <- tumap(mnist, init = "spca", batch = TRUE, rng_type = "deterministic", seed = 42)
```

If you use `nn_method = "hnsw"` (or `"nndescent"`) then also add
`n_build_threads = 1` to the above calls.
