## require(Matrix)

## Internal helper functions

############################## Vector Operations ##############################
#' Vector Euclidean norm
vecnorm <- function(x) sqrt(sum(x^2))

#' Normalize a vector in Euclidean norm
normalize <- function(x) x / sqrt(sum(x^2))

#' Distance vector from a point to a set of points
#' @param x a set of points, either a length-n vector or an n-by-p matrix
#' @param x0 a point, either a scalar or a length-p vector
distvec <- function(x, x0) {
    if (isMatrix(x)) {
        dx <- x - rep(x0, each = nrow(x))
        dist <- apply(dx, 1, vecnorm)
    } else {
        dist <- abs(x - x0)
    }
    return(dist)
}

#' Geometric mean
#' @param x A numeric vector of non-negative numbers.
#' @note `prod()` can under or overflow very quickly.
geometric_mean <- function(x, NA.rm = FALSE){
    if(NA.rm) x <- x[!is.na(x)]
    if(any(x < 0)) return(NaN)
    if(any(x == 0)) return(0)
    exp(mean(log(x)))
}


############################## Matrix Operations ##############################
#' Matrix of ones
#' @return A `dpoMatrix` of order n, all entries are 1.
J <- function(n) tcrossprod(Matrix::Matrix(rep(1, n)))

#' Compute the Hadamard product of XtX and kronecker(S, 1_k 1_k^T),
#' without computing the Kronecker product.
#' @note Much slower than XtX * kronecker(DKDinv, Matrix(1, k, k)).
blockHadamard <- function(XtX, S) {
    M <- Matrix::Matrix(NA_real_, k * l, k * l)
    for(i in seq(l)) {
        for(j in seq(l)) {
            rows <- (i-1) * k + seq(k)
            cols <- (j-1) * k + seq(k)
            M[rows, cols] <- XtX[rows, cols] * S[i,j]
        }
    }
    return(M)
}

## #' Whether a `dgCMatrix` is diagonal.
## isDiagonal <- function(M) {
##     ## Row indices of nonzero elements are sequential.
##     identical(M@i, seq(0, nrow(M) - 1)) &
##         ## Each column has one nonzero element, and its a square matrix.
##         identical(M@p, seq(0, nrow(M)))
## }

#' Whether a `dgCMatrix` is a permutation.
isPermutation <- function(M) {
    ## Nonzero elements are all ones.
    all(M@x == 1) &
        ## Each column has one nonzero element, and its a square matrix.
        identical(M@p, seq(0, nrow(M))) &
        ## Row indices of nonzero elements form a permutation
        setequal(M@i, seq(0, nrow(M) - 1))
}

#' Indicator function of a matrix object
isMatrix <- function(x) {
    is.matrix(x) | is(x, "Matrix")
}

#' Collect a list of matrices into an array.
listMatrix2Array <- function(listMr) {
    nr <- nrow(listMr[[1]])
    nc <- ncol(listMr[[1]])
    vapply(listMr, function(x) as.matrix(x), matrix(NA_real_, nrow = nr, ncol = nc))
}

#' Dummy block product: Boxdot = [r_ij Box_ij]_{i,j = 1, ..., l}
#' @note For the particular uses in this package,
#' `Box` and `R` are assumed to be symmetric matrices.
dummyBlockProd <- function(Box, R, ks) {
    k <- nrow(Box)
    l <- nrow(R)
    stopifnot(sum(ks) == k)
    ## stopifnot(all.equal(diag(R), rep(1, l)))
    id0 <- c(0, cumsum(ks)[-l])
    ## lsID <- purrr::map2(id0, ks, ~seq(.x + 1, .x + .y))
    lsID <- lapply(seq(l), function(i) seq(id0[i] + 1, id0[i] + ks[i]))
    ## NOTE: Can be made more efficient for symmetric matrices.
    Boxdot <- matrix(NA_real_, k, k)
    DTid <- data.table::CJ(i = seq(l), j = seq(l))
    setBlock <- function(i, j) {
        Boxdot[lsID[[i]], lsID[[j]]] <<- Box[lsID[[i]], lsID[[j]]] * R[i,j]
        return(NULL)
    }
    ## purrr::pwalk(DTid, setBlock)
    ## DTid[, lapply(seq(l^2), function(x) setBlock(i[x], j[x]))]
    with(DTid, lapply(seq(l^2), function(x) setBlock(i[x], j[x])))
    return(as(Boxdot, "dsyMatrix"))
}

############################## Optimization ##############################
#' Readable tolerance level for `optim()`.
#' @param x Precision for optimization.
ToleranceLevel <- function(x) list(factr = 10^(15-x))

#' Safely set negative (eigen-)values to zero.
#' @param eigs A vector of real eigenvalues in decreasing order.
safeNonNeg <- function(eigs, factor = 3) {
    if (any(eigs / eigs[1] < - factor * .Machine$double.eps)) {
        message("Large negative eigenvalues")
        stop()
    }
    pmax(eigs, 0)
}
