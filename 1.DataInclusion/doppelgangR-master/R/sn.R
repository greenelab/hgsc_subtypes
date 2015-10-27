dst <-  function (x, location=0, scale=1, shape=0, df=Inf, dp=NULL, log=FALSE)
{
    if(!is.null(dp)) {
        if(!missing(shape))
            stop("You cannot set both component parameters and dp")
        location <- dp[1]
        scale <- dp[2]
        shape <- dp[3]
        df <- dp[4]
    }
    if(df > 1e6 ) return(dsn(x,location,scale,shape, log=log))
    z   <- (x - location)/scale
    pdf <- dt(z, df=df, log=log)
    cdf <- pt(shape*z*sqrt((df+1)/(z^2+df)), df=df+1, log.p=log)
    if(log)
        logb(2) + pdf + cdf -logb(scale)
    else
        2 * pdf * cdf / scale
}


rst <- function (n=1, location = 0, scale = 1, shape = 0, df=Inf, dp=NULL)
{
    if(!is.null(dp)) {
        if(!missing(shape))
            stop("You cannot set both component parameters and dp")
        location <- dp[1]
        scale <- dp[2]
        shape <- dp[3]
        df <- dp[4]
    }
    z <- rsn(n,location=0,scale,shape)
    if(df==Inf) return(z+location)
    v <- rchisq(n,df)/df
    y <- z/sqrt(v)+location
    attr(y,"parameters")<- c(location,scale,shape,df)
    return(y)
}

pst <- function (x, location=0, scale=1, shape=0, df=Inf, dp=NULL, ...)
{
    if(!is.null(dp)) {
        if(!missing(shape))
            stop("You cannot set both component parameters and dp")
        location <- dp[1]
        scale <- dp[2]
        shape <- dp[3]
        df <- dp[4]
    }
    fp <- function(v, shape, df, t.value)
        psn(sqrt(v) * t.value, 0, 1, shape) * dchisq(v * df, df = df) * df
    if (df > 1e6) # (== Inf)
        p <- psn(x, location, scale, shape)
    else
        {
            if(df <= 0) stop("df must be non-negative")
            z <- (x-location)/scale
            p <- numeric(length(z))
            for (i in 1:length(z)){
                p[i]  <-
                    if(round(df)==df)
                        pmst(z[i], 0, matrix(1,1,1), shape, df, ...)
                    else{
                        if(abs(z[i]) == Inf)   (1+sign(z[i]))/2
                        else{
                            if(z[i] < 0)
                                integrate(dst, -Inf, z[i], shape = shape, df = df, ...)$value
                            else
                                integrate(fp, 0, Inf, shape = shape, df = df, t.value = z[i], ...)$value
                        }}
            }
            pmax(0,pmin(1,p))
        }
}

qst <- function (p, location = 0, scale = 1, shape = 0, df=Inf,
                 tol = 1e-06, dp = NULL, ...)
{
    if(!is.null(dp)) {
        if(!missing(shape))
            stop("You cannot set both component parameters and dp")
        location <- dp[1]
        scale <- dp[2]
        shape <- dp[3]
        df <- dp[4]
    }
    if (df > 1e4) # (== Inf)
        return(qsn(p, location, scale, shape))
    max.q <- sqrt(qf(p, 1, df))
    min.q <- -sqrt(qf(1 - p, 1, df))
    if (shape == Inf)
        return(location + scale * max.q)
    if (shape == -Inf)
        return(location + scale * min.q)
    na <- is.na(p) | (p < 0) | (p > 1)
    zero <- (p == 0)
    one <- (p == 1)
    p <- replace(p, (na | zero | one), 0.5)
    cum <- st.cumulants(0, 1, shape, max(df,5), n=4)
    g1 <- cum[3]/cum[2]^(3/2)
    g2 <- cum[4]/cum[2]^2
    x <- qnorm(p)
    x <- (x + (x^2 - 1) * g1/6 + x * (x^2 - 3) * g2/24 -
          x * (2 *  x^2 - 5) * g1^2/36)
    x <- cum[1] + sqrt(cum[2]) * x
    max.err <- 1
    while (max.err > tol) {
        x1 <- x - (pst(x, 0, 1, shape, df, ...) - p)/dst(x, 0, 1, shape, df)
        x1 <- pmin(x1, max.q)
        x1 <- pmax(x1, min.q)
        max.err <- max(abs(x1 - x)/(1 + abs(x)))
        x <- x1
    }
    x <- replace(x, na, NA)
    x <- replace(x, zero, -Inf)
    x <- replace(x, one, Inf)
    return(as.numeric(location + scale * x))
}

st.mle <- function(X, y, freq,  start, fixed.df=NA, trace=FALSE,
                   algorithm = c("nlminb","Nelder-Mead", "BFGS", "CG", "SANN"),
                   control=list())
{
    y.name  <- deparse(substitute(y))
    y <- data.matrix(y)
    if(missing(X)) X<- matrix(1, nrow=length(y), ncol=1)
    dimnames(y)[[2]] <- list(y.name)
    if(missing(start)){
        cp0 <- sn.mle(X=X, y=y, plot.it=FALSE, trace=trace)$cp
        m <- length(cp0)-2
        cp0[m+2] <- cp0[m+2]*0.9
        mle0 <- cp.to.dp(cp0)
        start <- list(beta=mle0[1:m], Omega=matrix(mle0[m+1]^2,1,1),
                      alpha=mle0[m+2], df=10)
    }
    else {
        m <- length(start)-3
        if(m<1) stop("bad start vector")
        start<-  list(beta=start[1:m], Omega=matrix(start[m+1]^2,1,1),
                      alpha=start[m+2], df=start[m+3])
    }
    fit <- mst.mle(X, y, freq, start=start, fixed.df=fixed.df, trace=trace,
                   algorithm=algorithm, control=control)
    mle <- list()
    mle$call<- match.call()
    dp <- fit$dp
    se <- fit$se
    p  <- length(dp$beta)
    dp.names <- c(if(p==1) "location" else dimnames(dp$beta)[[1]],
                  "scale","shape","df")
    mle$dp  <- c(dp$beta, sqrt(as.vector(dp$Omega)), dp$alpha, dp$df)
    names(mle$dp) <- dp.names
    mle$se <- if(all(is.na(se))) NA else
    c(se$beta, mle$dp[p + 1] * se$internal[p + 1],
      se$alpha, dp$df * se$internal[p + 3])
    mle$logL <- fit$logL
    mle$algorithm <- fit$algorithm
    mle
}


mst.mle <- function (X, y, freq, start, fixed.df = NA, trace = FALSE,
                     algorithm = c("nlminb", "Nelder-Mead", "BFGS", "CG", "SANN"),
                     control = list())
{
    algorithm <- match.arg(algorithm)
    y.name <- deparse(substitute(y))
    y.names <- dimnames(y)[[2]]
    y <- data.matrix(y)
    X <- if (missing(X)) matrix(rep(1, nrow(y)), ncol = 1)
    else data.matrix(X)
    if (missing(freq)) freq <- rep(1, nrow(y))
    x.names <- dimnames(X)[[2]]
    d <- ncol(y)
    n <- sum(freq)
    m <- ncol(X)
    if (missing(start)) {
        qrX <- qr(X)
        beta <- as.matrix(qr.coef(qrX, y))
        Omega <- matrix(var(qr.resid(qrX, y)), d, d)
        omega <- sqrt(diag(Omega))
        alpha <- rep(0, d)
        df <- ifelse(is.na(fixed.df), 10, fixed.df)
        if (trace) {
            cat("mst.mle: dp=", "\n")
            print(c(beta, Omega, alpha))
            cat("df:", df, "\n")
        }
    }
    else {
        if (!is.na(fixed.df))
            start$df <- fixed.df
        if (all(names(start) == c("beta", "Omega", "alpha", "df"))) {
            beta <- start$beta
            Omega <- start$Omega
            alpha <- start$alpha
            df <- start$df
        }
        else stop("start parameter is not in the form that I expected")
    }
    eta <- alpha/sqrt(diag(Omega))
    Oinv <- solvePD(Omega)
    upper <- chol(Oinv)
    D <- diag(upper)
    A <- upper/D
    D <- D^2
    if (d > 1)
        param <- c(beta, -log(D)/2, A[!lower.tri(A, diag = TRUE)], eta)
    else
        param <- c(beta, -log(D)/2, eta)
    if (is.na(fixed.df))
        param <- c(param, log(df))
    if(algorithm == "nlminb"){
        opt <- nlminb(param, objective = mst.dev, gradient = mst.dev.grad,
                      control = control,  X = X, y = y, freq = freq,
                      trace = trace, fixed.df = fixed.df)
        info <- num.deriv2(opt$par, FUN="mst.dev.grad", X=X, y=y,
                           freq=freq, fixed.df = fixed.df)/2
        opt$value <-  opt$objective
    }
    else{
        opt <- optim(param, fn = mst.dev, gr = mst.dev.grad,
                     method = algorithm, control = control, hessian = TRUE,
                     X = X, y = y, freq = freq, trace = trace, fixed.df = fixed.df)
        info <- opt$hessian/2
    }
    dev   <- opt$value
    param <- opt$par
    opt$name <- algorithm
    if (trace) {
        cat("Message from optimization routine:", opt$message, "\n")
        cat("deviance:", dev, "\n")
    }
    beta <- matrix(param[1:(m * d)], m, d)
    D <- exp(-2 * param[(m * d + 1):(m * d + d)])
    A <- diag(d)
    i0 <- m*d+d*(d+1)/2
    if(d>1)  A[!lower.tri(A,diag=TRUE)] <- param[(m*d+d+1):i0]
    eta <- param[(i0 + 1):(i0 + d)]
    if (is.na(fixed.df))
        df <- exp(param[i0 + d + 1])
    else df <- fixed.df
    Oinv <- t(A) %*% diag(D,d,d) %*% A
    Omega <- solvePD(Oinv)
    omega <- sqrt(diag(Omega))
    alpha <- eta * omega
    dimnames(beta) <- list(x.names, y.names)
    dimnames(Omega) <- list(y.names, y.names)
    if (length(y.names) > 0) names(alpha) <- y.names
    if (all(is.finite(info))) {
        qr.info <- qr(info)
        info.ok <- (qr.info$rank == length(param))
    }
    else info.ok <- FALSE
    if (info.ok) {
        se2 <- diag(solve(qr.info))
        if (min(se2) < 0)
            se <- NA
        else {
            se <- sqrt(se2)
            se.beta <- matrix(se[1:(m * d)], m, d)
            se.alpha <- se[(i0 + 1):(i0 + d)] * omega
            dimnames(se.beta)[2] <- list(y.names)
            dimnames(se.beta)[1] <- list(x.names)
            names(se.alpha) <- y.names
            se.df <- df * se[i0 + d + 1]
            se <- list(beta = se.beta, alpha = se.alpha, df = se.df,
                       internal = se, info = info)
        }
    }
    else se <- NA
    dp <- list(beta = beta, Omega = Omega, alpha = alpha, df = df)
    list(call = match.call(), logL = -dev/2, deviance = dev,
         dp = dp, se = se, algorithm = opt)
}
