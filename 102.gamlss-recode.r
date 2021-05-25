##
## This script contains re-writes of gamlss functions
##
## NOTE: These could be considered as bugs in gamlss (or at least possible improvements in the package codebase
##

##
## We must load gamlss-package, since we directly copy some functions/objects from it
library("gamlss") ## need gamlss.family objects in all scripts



##
## bfp() := the in-built gamlss function cannot handle NAs
##
bfpNA <- function (x, powers = c(1, 2), shift = 0, scale = 1) ## have change gamlss:bfp() defaults, now fix shift=0 and scale=1 [DATA DEPENDENT SCALING IS BAD!]
{
    ## fp.scale <- function(x) {
    ##     if( all(is.na(x)) ) {stop("All NAs makes no sense")}
    ##     if (min(x,na.rm=TRUE) <= 0) {
    ##         xx <- na.omit(x)
    ##         z <- sort(xx)[-1] - sort(xx)[-length(xx)]
    ##         shift <- min(z[z > 0]) - min(xx)
    ##     }
    ##     else shift <- 0
    ##     range <- max(x,na.rm=TRUE) - min(x,na.rm=TRUE)
    ##     scale <- 10^(sign(log10(range)) * trunc(abs(log10(range))))
    ##     list(shift = shift, scale = scale)
    ## }
    nobs <- length(x)
    npoly <- length(powers)
    X <- matrix(0, nrow = nobs, ncol = npoly)
    if (is.null(scale) | is.null(shift)) {
        stop("WARNING: Using automatic scale/shift will invalidate future refitting")
        out <- fp.scale(x)
        shift <- out$shift
        scale <- out$scale
    }
    ## x <- x + shift ## ASSUME variable is validly scaled and shifted!
    ## x <- x/scale
    x1 <- ifelse(powers[1] != rep(0, nobs), x^powers[1], log(x))
    X[, 1] <- x1
    if (npoly >= 2) {
        for (i in 2:npoly) {
            if (powers[i] == powers[(i - 1)]) 
                x2 <- log(x) * x1
            else x2 <- ifelse(powers[i] != rep(0, nobs), x^powers[i], 
                log(x))
            X[, i] <- x2
            x1 <- x2
        }
    }
    X
}

##
## bfp() := the in-built gamlss function cannot handle NAs
##          we replace the in-built with a warning, to make sure we do not accidentally use the 'broken' version
bfp <- function( ... ) {stop("Default bfp() function cannot handle NAs. We have masked with this fatal error. Use bfpNA() instead. ")}



##
## GG() := the GG-family includes variance() and mean() functions that are prone to buffer overflow/underflow (too big/too small numbers)
##         we replace the GG family with GGalt
##

dGGalt <- dGG
pGGalt <- pGG
qGGalt <- qGG
rGGalt <- rGG

##
## The following is a direct copy of the GG functions, except where alterations are indicated
GGalt <- function (mu.link = "log", sigma.link = "log", nu.link = "identity") 
{
    mstats <- checklink("mu.link", "GGalt", substitute(mu.link), c("1/mu^2", "log", "identity"))
    dstats <- checklink("sigma.link", "GGalt", substitute(sigma.link), c("inverse", "log", "identity"))
    vstats <- checklink("nu.link", "GGalt", substitute(nu.link), c("1/nu^2", "log", "identity"))

    ## Replacement mean and variance functions
    GGalt.mean <- function (mu, sigma, nu) {
        TOP <- log(mu) + lgamma(1/(sigma^2 * nu^2) + 1/nu)
        ##BOTTOM <-  log((1/(sigma^2 * nu^2))^(1/nu)) + lgamma(1/(sigma^2 * nu^2)) [fixed in v10]
        BOTTOM <-  (-1)*((1/nu)*log((sigma^2 * nu^2))) + lgamma(1/(sigma^2 * nu^2))
        ifelse(nu > 0 | (nu < 0 & sigma^2 * abs(nu) < 1),
               exp(TOP-BOTTOM),
               Inf)
    }
    GGalt.variance <- function (mu, sigma, nu) {
        ##AA <- log(mu^2) - log( (1/(sigma^2 * nu^2))^(2/nu) ) - 2*lgamma(1/(sigma^2 * nu^2)) [fixed in v10]
        AA <- 2*log(mu) - ((2/nu)*(-1)*log( (sigma^2 * nu^2) )) - 2*lgamma(1/(sigma^2 * nu^2))
        ww <- lgamma(1/(sigma^2 * nu^2) + 2/nu) + lgamma(1/(sigma^2 * nu^2))
        uu <- 2*lgamma(1/(sigma^2 * nu^2) + 1/nu)
        BB <- ww + log( (1 - exp( uu - ww )) )
        YES <- AA + BB

        ifelse(nu > 0 | (nu < 0 & sigma^2 * abs(nu) < 0.5),
        ifelse(is.nan(BB),NA,exp( YES )),
        Inf)
    }

    
    structure(list(family = c("GGalt", "generalised Gamma Lopatatsidis-Green (altered)"),  ##EDIT
                   parameters = list(mu = TRUE, sigma = TRUE, nu = TRUE), 
                   nopar = 3,
                   type = "Continuous",
                   ##
                   mu.link = as.character(substitute(mu.link)), 
                   sigma.link = as.character(substitute(sigma.link)),
                   nu.link = as.character(substitute(nu.link)), 
                   ##
                   mu.linkfun = mstats$linkfun,
                   sigma.linkfun = dstats$linkfun, 
                   nu.linkfun = vstats$linkfun,
                   ##
                   mu.linkinv = mstats$linkinv, 
                   sigma.linkinv = dstats$linkinv,
                   nu.linkinv = vstats$linkinv,
                   ##
                   mu.dr = mstats$mu.eta,
                   sigma.dr = dstats$mu.eta,
                   nu.dr = vstats$mu.eta,
                   ##
                   dldm = function(y, mu, sigma, nu) {
                       z <- (y/mu)^nu
                       theta <- 1/(sigma^2 * abs(nu)^2)
                       dldm <- ifelse(abs(nu) > 1e-06, (z - 1) * theta * nu/mu, (1/(mu * (sigma^2)) * (log(y) - log(mu))))
                       dldm },
                   d2ldm2 = function(mu, sigma, nu) {
                       d2ldm2 <- ifelse(abs(nu) > 1e-06, -1/((mu^2) * (sigma^2)), -(1/(mu^2 * sigma^2)))
                       d2ldm2 },
                   dldd = function(y, mu, sigma, nu) {
                       z <- (y/mu)^nu
                       theta <- 1/(sigma^2 * abs(nu)^2)
                       dldd <- ifelse(abs(nu) > 1e-06, -2 * theta * (log(theta) + 1 + log(z) - z - digamma(theta))/sigma, -(1/sigma) + (1/sigma^3) * (log(y) - log(mu))^2)
                       dldd },
                   d2ldd2 = function(y, mu, sigma, nu) {
                       theta <- 1/(sigma^2 * abs(nu)^2)
                       d2ldd2 <- ifelse(abs(nu) > 1e-06, 4 * (theta/(sigma^2)) * (1 - theta * trigamma(theta)), -2/sigma^2)
                       d2ldd2 },
                   dldv = function(y, mu, sigma, nu) {
                       z <- (y/mu)^nu
                       theta <- 1/(sigma^2 * abs(nu)^2)
                       dldv <- (1/nu) * (1 + 2 * theta * (digamma(theta) + z - log(theta) - 1 - ((z + 1)/2) * log(z)))
                       dldv },
                   d2ldv2 = function(y, mu, sigma, nu) {
                       theta <- 1/(sigma^2 * abs(nu)^2)
                       d2ldv2 <- -(theta/nu^2) * (trigamma(theta) * (1 + 4 * theta) - (4 + 3/theta) - log(theta) * (2/theta - log(theta)) + digamma(theta) * (digamma(theta) + (2/theta) - 2 * log(theta)))
                       d2ldv2 },
                   d2ldmdd = function(y) rep(0, length(y)),
                   d2ldmdv = function(y, mu, sigma, nu) {
                       theta <- 1/(sigma^2 * abs(nu)^2)
                       ddd <- (theta/mu) * (digamma(theta) + (1/theta) - log(theta))
                       ddd },
                   d2ldddv = function(y, mu, sigma, nu) {
                       theta <- 1/(sigma^2 * abs(nu)^2)
                       d2ldddv <- -2 * sign(nu) * theta^(3/2) * (2 * theta * trigamma(theta) - (1/theta) - 2)
                       d2ldddv },
                   G.dev.incr = function(y, mu, sigma, nu, ...) -2 * dGGalt(y, mu = mu, sigma = sigma, nu = nu, log = TRUE), ## EDIT
                   ##
                   rqres = expression(rqres(pfun = "pGGalt", type = "Continuous", y = y, mu = mu, sigma = sigma, nu = nu)), ## EDIT
                   ##
                   mu.initial = expression(mu <- (y + mean(y))/2),
                   sigma.initial = expression(sigma <- rep(1, length(y))),
                   nu.initial = expression(nu <- rep(1, length(y))), mu.valid = function(mu) all(mu > 0),
                   ##
                   sigma.valid = function(sigma) all(sigma > 0),
                   nu.valid = function(nu) TRUE, 
                   y.valid = function(y) all(y > 0),
                   ##
                   ## ***** replaced elements of list below with improved mean and variance functions
                   ##
                   mean = GGalt.mean, ## EDIT
                   variance = GGalt.variance ## EDIT
                   ##
                   ##
                   ##
                   ),
              class = c("gamlss.family", "family")
              )
}

