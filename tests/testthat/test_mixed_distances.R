library(uwot)
context("mixed distance calculations")

set.seed(1337)
res <- umap(iris10, n_neighbors = 4, n_epochs = 2, init = "spca",
            verbose = FALSE, n_threads = 0)
expect_ok_matrix(res)

set.seed(1337)
resmli <- umap(iris10, n_neighbors = 4, n_epochs = 2, init = "spca",
               metric = list("euclidean" = 1:4),
               verbose = FALSE, n_threads = 0)
expect_equal(resmli, res)

set.seed(1337)
resmls <- umap(iris10, n_neighbors = 4, n_epochs = 2, init = "spca",
               metric = list("euclidean" = c("Sepal.Length", "Sepal.Width",
                                             "Petal.Length", "Petal.Width")),
               verbose = FALSE, n_threads = 0)
expect_equal(resmls, res)

set.seed(1337)
jiris10 <- jitter(iris10)
metric2 <- list("euclidean" = c(1, 2),
               "euclidean" = c("Petal.Length", "Petal.Width"))
reseuc2 <- umap(jiris10, n_neighbors = 4, n_epochs = 2, init = "spca",
               metric = metric2,
               verbose = FALSE, n_threads = 0,
               ret_nn = TRUE, ret_model = TRUE)
expect_ok_matrix(reseuc2$embedding)
expect_equal(reseuc2$metric, metric2)
expect_is(reseuc2$nn, "list")
expect_equal(names(reseuc2$nn), c("euclidean", "euclidean"))
expect_is_nn(reseuc2$nn[[1]], 10, 4)
expect_is_nn(reseuc2$nn[[2]], 10, 4)
expect_ok_matrix(umap_transform(jiris10, reseuc2))

i10factor <- factor(c(rep("foo", 3), rep("bar", 3), rep("baz", 4)))

res_y2 <- umap(iris10[, -1], y = cbind(i10factor, iris$Sepal.Length[1:10]),
               n_neighbors = 4, n_epochs = 2, init = "spca",
               verbose = FALSE, n_threads = 0)
expect_ok_matrix(res_y2)
