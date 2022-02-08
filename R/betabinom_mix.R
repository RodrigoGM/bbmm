## Author: Ronglai Shen <shenr@mskcc.org>
betabinom.mix <- function(x, n) {
    a <- c(0.01, 0.99, 0.99)
    b <- c(0.99, 0.99, 0.01)
    p <- c(1/3, 1/3, 1/3)
    
    iter <- 1; dif=1
    while(dif > 1e-3 & iter < 20){	
        cat(iter, '\n')
        a.old <- a; b.old <- b
        dg0 <- dbetabinom.ab(x, size = n, shape1 = a[1], shape2 = b[1])
        dg1 <- dbetabinom.ab(x, size = n, shape1 = a[2], shape2 = b[2])
        dg2 <- dbetabinom.ab(x, size = n, shape1 = a[3], shape2 = b[3])
        dg <- cbind(dg0, dg1, dg2)
        mdg <- apply(dg, 1, which.max)
        
        ## E-step posterior probability for each genotype
        denom <- dg0*p[1] + dg1*p[2] + dg2*p[3]
        postg0 <- dg0*p[1]/denom
        postg1 <- dg1*p[2]/denom
        postg2 <- dg2*p[3]/denom
        
        ## mixing proportion
        p <- c(mean(postg0), mean(postg1), mean(postg2))
        
        ## fit = vglm(cbind(x0, n0-x0) ~ 1, betabinomial, trace = TRUE)
        ## coef(fit)
        
        ## M-step
        bb_mle <- function(x, n) {
            ll <- function(alpha, beta) {
                -sum(dbetabinom.ab(x, n, alpha, beta, log = TRUE))
            }
            m <- stats4::mle(ll, start = list(alpha = 3, beta = 10), method = "L-BFGS-B", 
                             lower = c(0.001, .001))
            ab <- stats4::coef(m)
            list(a = ab[1], b = ab[2], number = length(x))
        }
        
        bb1 <- bb_mle(x[mdg == 1], n[mdg == 1])
        a[1] <- bb1$a; b[1] <- bb1$b
        
        bb2 <- bb_mle(x[mdg == 2], n[mdg == 2])
        a[2] <- bb2$a; b[2] <- bb2$b
        
        bb3 <- bb_mle(x[mdg == 3], n[mdg == 3])
        a[3] <- bb3$a; b[3] <- bb3$b
        
        dif <- max(abs(a-a.old), abs(b-b.old))
        iter <- iter+1
    }
    
    out <- list(genotype = mdg, pAA = postg0, pAB = postg1, pBB = postg2, shape1 = a, shape2 = b, mixp = p)
    return(out)
}
