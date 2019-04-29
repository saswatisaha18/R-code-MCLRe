if (getRversion() < "2.14") library(multicore) else library(parallel)
# just like `sapply'
mcsapply <- function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE)
{
    FUN <- match.fun(FUN)
    answer <- mclapply(X = X, FUN = FUN, ...)
    if (USE.NAMES && is.character(X) && is.null(names(answer)))
        names(answer) <- X
    if (!identical(simplify, FALSE) && length(answer))
        simplify2array(answer, higher = (simplify == "array"))
    else answer
}
# just like `replicate'
mcreplicate <- function(n, expr, simplify = "array")
    mcsapply(integer(n), eval.parent(substitute(function(...) expr)),
             simplify = simplify)
# expr should evaluate to TRUE or FALSE;
# return proportion of TRUE's in n.rep \times n.sim replications.
mcprop <- function(expr, n.rep = 10, n.sim = 1e5)
{
    foo <- eval.parent(substitute(function(...) expr))
    mean(mcsapply(integer(n.rep), function(...)
                  mean(sapply(integer(n.sim), foo))))
}
