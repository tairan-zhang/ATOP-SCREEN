design_env <- new.env(parent = globalenv())
design_env$geomean_custom_for_gsea <- function(x) exp(mean(log(x[x > 0])))
design_env$`%||%` <- function(x,y) if (is.null(x)) y else x
for (file in c("preprocessing.R","validation.R")) {
    sys.source(testthat::test_path("..","..","R","atop",file), design_env)
}
counts <- data.frame(guide=letters[1:6], gene=rep(c("A","B"),each=3),
                     end1=c(9,12,1,40,6,30), end2=c(11,40,4,30,9,8),
                     base1=c(10,6,20,15,2,10), base2=c(8,7,15,30,8,11))
score_design <- function(mode, ends, bases) {
    p <- design_env$atop_validate_request(counts,"guide","gene",NULL,ends,bases,0,1,1,"cpp",FALSE,NULL,NULL,0,mode)
    p$ncs_delta0 <- .1; p$ncs_delta1 <- .3; p$ncs_delta2 <- .2
    prepared <- design_env$atop_prepare_input(counts,p,function(...) NULL)
    design_env$atop_score_guides(prepared,p,function(...) NULL)
}
test_that("shared baseline retains endpoint contrasts without duplicating baseline normalization", {
    result <- score_design("auto", c("end1","end2"), "base1")
    n <- result$normalized
    expect_identical(names(n),c("guide","gene","end1","end2","base1"))
    x <- log2((n$end1+1)/(n$base1+1))
    y <- log2((n$end2+1)/(n$base1+1))
    expect_equal(result$ranked$deltaLFC_Rep1,x)
    expect_equal(result$ranked$deltaLFC_Rep2,y)
    expect_equal(result$ranked$diff_score,(x+y)/2)
    expect_true("NCS" %in% names(result$ranked))
})
test_that("unpaired LFC uses normalized arithmetic group means and has no replicate NCS", {
    result <- score_design("auto","end1",c("base1","base2"))
    n <- result$normalized
    expect_equal(result$ranked$diff_score,log2((n$end1+1)/((n$base1+n$base2)/2+1)))
    expect_false("NCS" %in% names(result$ranked))
    groups <- score_design("group_means",c("end1","end2"),c("base1","base2"))
    n <- groups$normalized
    expect_equal(groups$ranked$diff_score,log2((rowMeans(n[c("end1","end2")])+1)/(rowMeans(n[c("base1","base2")])+1)))
    expect_error(score_design("paired","end1",c("base1","base2")),"equally sized")
})
