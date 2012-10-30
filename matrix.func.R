RemoveSize <- function (cov.matrix){
  # Removes first principal component effect in cov.matrix.
  #
  # Args:
  #   cov.matrix: A simetric covariance matrix
  # Return:
  #   cov.matrix with size removed
  cov.matrix.svd  <-  svd(cov.matrix)
  size.eigen.vector <- cov.matrix.svd$u[, 1]
  size.eigen.value <- cov.matrix.svd$d[1]
  size.factor <- size.eigen.vector * sqrt(size.eigen.value)
  cov.matrix.size.removed <- cov.matrix - size.factor %*% t(size.factor)
  return (cov.matrix.size.removed)
}

TestModularity <- function (cor.matrix, modularity.hipot) {
  # Tests modularity hipotesis using cor.matrix matrix and trait groupings
  #
  # Args:
  #   cor.matrix: Correlation matrix
  #   modularity.hipot: Matrix of hipotesis.
  #                     Each line represents a trait and each column a module.
  #                     if modularity.hipot[i,j] == 1, trait i is in module j.
  # Return:
  #   list with mantel correlation of cor.matrix with binary hipotesis matrices
  no.hip <- dim (modularity.hipot) [2]
  traits <- dim (modularity.hipot) [1]
  m.hip.array <- array (0, c(traits, traits, no.hip + 1))
  for (N in 1:no.hip){
    for (L in 1:traits){
      for (M in 1:traits){
        m.hip.array[L,M,N] <- ifelse (modularity.hipot[L,N] & modularity.hipot[M,N], 1, 0)
      }
    }
  }
  m.hip.array[,,no.hip+1] <- as.integer (as.logical (apply (m.hip.array, c(1,2), sum)))
  no.hip <- no.hip + 1
  output <- array (0, c(5, no.hip))
  for (N in 1:no.hip){
    tmp <- MantelCor (cor.matrix, m.hip.array[,,N], mod = TRUE)
    output[,N] <- tmp
  }
  dimnames(output) <- list (names (tmp), c (colnames (modularity.hipot),"Full Integration"))
  return (output)
}

MultiRsMantel <- function (matrix.list, MatrixCompFunc = RandomSkewers,
                           repeat.vector = NULL, iterations = 1000){
  # Performs multiple comparisons between a set of covariance of correlation matrices.
  #
  # Args:
  #  matrix.list: a list of covariance or correlation matrices
  #  MatrixCompFunc: function to use for comparison
  #  repeat.vector: vector of matrix repeatabities
  #  iterations: number of RandomSkewers of matrix permutations passed to MatrixCompFunc
  #
  # Return:
  #  a list with two matrices containing $\Gamma$-values or average random
  #  skewers correlation and probabilities according to permutation test.
  #  if repeat.vector was also passed, values above the diagonal on the correlation matrix
  #  will contain corrected correlation values.
  matrix.n <- length (matrix.list)
  matrix.names <- names (matrix.list)
  probabilities <- array (0, c(matrix.n, matrix.n))
  correlations <- probabilities
  for (i in 1:(matrix.n - 1)) {
    for (j in (i+1):matrix.n) {
      cat (i, ' ', j, '\n')
      comparing.now <- MatrixCompFunc (matrix.list [[i]],
                                       matrix.list [[j]],
                                       iterations)
      correlations [i, j] <- comparing.now [1]
      probabilities [i, j] <- comparing.now [2]
      if (!is.null (repeat.vector))
        correlations [j, i] <- correlations [i, j] / sqrt (repeat.vector [i] * repeat.vector [j])
    }
  }
  if (!is.null (repeat.vector)) {
    repeat.mat <- sqrt (repeat.vector %*% t (repeat.vector))
    corrected <- correlations / repeat.mat
    correlations <- t (corrected) + correlations
    diag (correlations) <- repeat.vector
  }
  rownames (correlations) <- matrix.names
  colnames (correlations) <- matrix.names
  dimnames (probabilities) <- dimnames (correlations)
  output <- list ('correlations' = correlations, 'probabilities' = probabilities)
  return (output)
}

measure.rep <- function (data1, data2, taille = 30) ## talvez eu devesse ajustar cada modelo
  {
    ind <- dim (data1) [1]
    traits <- dim (data1) [2]
    sel <- sample (1:ind, taille)
    subset1 <- data1 [sel,]
    subset2 <- data2 [sel,]
    fac <- c (rownames (subset1), rownames (subset2))
    fac <- factor (fac, levels = unique (fac))
    maindata <- rbind (subset1, subset2)
    reps <- c ()
    for (N in 1:traits)
      {
        msq <- anova (lm (maindata[,N]~fac))[,3]
        reps[N] <- msq[1] / sum(msq)
      }
    names (reps) <-  dimnames (data1) [[2]]
    output <- list ("Individuals" = levels(fac),"Repetabilities" = reps)
    return (output)
  }

adjust.sex.age <- function (data, sex, age, show.lm = FALSE)
  {
    ind <- dim(data)[1];traits <- dim(data)[2]
    arr <- array(0,c(ind,traits))
    for (N in 1:traits)
      arr[,N] <- data[,N]
    if (length(unique(sex)) == 2)
      {
        sex <- factor(sex, levels = unique(sex))
        put <- lm (arr ~ age * sex)
      }
    else
      put <- lm (arr ~ age)
    out <- residuals (put)
    dimnames (out) <- dimnames (data)
    if (show.lm == TRUE)
      return (list("Residuals" = out, "Models" = put))
    else
      return (out)
  }

bootstrap.rep <- function (data, nb = 100)
  {
    ind <-  dim (data) [1]
    or.vcv <- var (data)
    v.rep <- c()
    for (N in 1:nb)
      {
        strap <- sample (1:ind, ind, TRUE)
        strap.vcv <- var (data[strap,])
        v.rep [N] <- RandomSkewers (or.vcv, strap.vcv, 1000) [1]
      }
    out <- mean (v.rep)
    return (out)
  }

bootstrap.rep.G <- function (data, sex, age, ind, nb = 1000, corr = FALSE)
  {
    IND <- dim (data) [1]
    or.res <- adjust.sex.age (data, sex, age)
    or.vcv <- var (or.res)
    v.rep <- c()
    if (corr)
      {
        or.cor <- cor (or.res)
        c.rep <- c()
      }
    for (N in 1:nb)
      {
        strap <- sample (1:IND, ind, TRUE)
        strap.res <- adjust.sex.age (data[strap,],sex[strap],age[strap])
        strap.vcv <- var (strap.res)
        v.rep [N] <- RandomSkewers (or.vcv, strap.vcv, 1) [1]
        if (corr)
          {
            strap.cor <- cor (strap.res)
            c.rep [N] <- MantelCor (or.cor, strap.cor, 1) [1]
          }
      }
    out <- mean (v.rep)
    if (corr)
      {
        put <- mean (c.rep)
        output <- c(out,put)
        names (output) <- c("VCV", "Corr")
        return (output)
      }
    else
      return (out)
  }

rmvNorm2 <- function (n, theta = rep(0, nrow(sigma)),
                      sigma = diag(length(theta)),
                      method = c("eigen", "svd", "chol")) {
  # Calculates random deviates from a Normal multivariate distribution
  #
  # Args:
  #   n: number os deviates
  #   theta: vetor of means
  #   sigma: covariance matrix
  #   method: generation method
  #
  # Return:
  #   Vector of deviates
  if (length(theta) != nrow(sigma)) {
    stop("theta and sigma have non-conforming size")
  }
  sigma.aux <- sigma
  dimnames(sigma.aux) <- NULL
  if (!isTRUE(all.equal(sigma.aux, t(sigma.aux)))) {
    warning("sigma is numerically not symmetric")
  }
  method <- match.arg(method) #  what is this?
  if (method == "eigen") {
    sigma.eigen <- eigen(sigma, symmetric = TRUE)
    if (!all(sigma.eigen$values >=
             -sqrt(.Machine$double.eps) * abs(sigma.eigen$values[1]))) {
      warning("sigma is numerically not positive definite")
    }
    random.deviates <- sigma.eigen$vectors %*%
                       diag(sqrt(sigma.eigen$values),
                            length(sigma.eigen$values)) %*%
                       t(sigma.eigen$vectors)
  }
  else if (method == "svd") {
    sigma.svd <- svd(sigma)
    if (!all(sigma.svd$d >= -sqrt(.Machine$double.eps) * abs(sigma.svd$d[1]))) {
      warning("sigma is numerically not positive definite")
    }
    random.deviates <- t(sigma.svd$v %*% (t(sigma.svd$u) * sqrt(sigma.svd$d)))
  }
  else if (method == "chol") {
    sigma.chol <- chol(sigma)
    random.deviates <- matrix(rnorm(n * ncol(sigma)), nrow = n) %*% sigma.chol
  }
  random.deviates <- sweep(random.deviates, 2, theta, "+")
  return(random.deviates)
}

monte.carlo.rep <- function (matrix, ind, nit = 100)
  {
    if (sum(diag(matrix)) == dim (matrix) [1])
      {
        Func <- MantelCor
        Type <- cor
      }
    else
      {
        Func <- RandomSkewers
        Type <- var
      }
    R <- c()
    for (N in 1:nit)
      {
        rand.samp <- rmvNorm2 (ind, rep(0, times = dim (matrix)[1]),
                              sigma = matrix, method = "chol")
        rand.matrix <- Type (rand.samp)
        R[N] <- Func (matrix, rand.matrix, 1000)[1]
      }
    return (mean(R))
  }

Norm = function(x){return(sqrt(sum(x*x)))}
Normalize = function(x){return(x/Norm(x))}

RandomSkewers <- function (cov.matrix.1, cov.matrix.2, nsk = 10000){
  # Calculates covariance matrix correlation via random skewers
  # Args:
  #     cov.matrix.(1,2): Two covariance matrices to be compared
  #     nsk: Number of generated random skewers
  # Result:
  #     List with mean value of correlation, p value and standard deviation
  traits <- dim (cov.matrix.1) [1]
  base.vector <- rnorm(traits)
  random.vectors <- array (rnorm (nsk * traits, mean = 0, sd = 1), c(traits, nsk))
  random.vectors <- apply (random.vectors, 2, Normlize)
  dist <- abs (base.vector %*% random.vectors)
  dz1 <- apply (cov.matrix.1 %*% random.vectors, 2, Normlize)
  dz2 <- apply (cov.matrix.2 %*% random.vectors, 2, Normlize)
  real <- abs (apply (dz1 * dz2, 2, sum))
  ac <- mean (real)
  stdev <- sd (real)
  prob <- sum (real < dist) / nsk
  output <- c(ac, prob, stdev)
  names(output) <- c("AC","P","SD")
  return(output)
}

MantelCor <- function (corr1, corr2, nit = 1000, mod = FALSE)
  {
    fix <- Normlize (corr1 [lower.tri (corr1)])
    mob <- Normlize (corr2 [lower.tri (corr2)])
    real <- cor (cbind (fix,mob))[1,2]
    dist <- c()
    for (N in 1:nit)
      {
        shuffle <- sample (1:dim(corr1)[1])
        mobb <- Normlize (corr2 [shuffle, shuffle] [lower.tri (corr1)])
        dist[N] <- cor (cbind (fix,mobb))[1,2]
      }
    prob <- sum (dist > as.vector(real)) / nit
    if (mod == TRUE)
      {
        avgM <- mean (corr1 [lower.tri(corr1)] [mob != 0])
        avgm <- mean (corr1 [lower.tri(corr1)] [mob == 0])
        avg.ratio <- avgM / avgm
        output <- c(real,prob,avgM,avgm,avg.ratio)
        names(output) <- c("R²","Probability","AVG+","AVG-","AVG Ratio")
      }
    else
      {
        output <- c(real,prob)
        names(output) <- c("R²","Probability")
      }
    return (output)
  }

papyrus <- function (vcv1, vcv2, nsk = 1000)
{
  size <- dim (vcv1)[1]
  r2s <- array (0, c(size,nsk))
  beta <- apply (array (rnorm (size*nsk, mean = 0, sd = 1),c(size,nsk)),2, Normlize)
  for (I in 1:nsk)
    {
      beta.matrix <- diag (beta[,I])
      dz1 <- apply (vcv1 %*% beta.matrix, 1, Normlize)
      dz2 <- apply (vcv2 %*% beta.matrix, 1, Normlize)
      r2s[,I] <- colSums (dz1 * dz2)
    }
  # results
  mean.r2 <- apply (r2s, 1, mean)
  sample.conf <- function (x, lower = TRUE)
    {
      ox <- x[order(x)]
      lox <- length (ox)
      if (lower)
        crit <- round (0.025 * lox)
      else
        crit <- round (0.975 * lox)
      return (ox[crit])
    }
  low.r2 <- apply (r2s, 1, sample.conf, lower = TRUE)
  up.r2 <- apply (r2s, 1, sample.conf, lower = FALSE)
  sd.r2 <- apply (r2s,1,sd)
  cmean.r2 <- scale (mean.r2, scale = FALSE)
  csd.r2 <- scale (sd.r2, scale = FALSE)
  cent <- cbind (cmean.r2,csd.r2)
  pca.cent <- princomp (cent, cor = TRUE, scores = TRUE)
  pc1 <- pca.cent$scores[,1]
  if (pca.cent$loadings[1,1] < 0)
    pc1 <- - pc1
  pc1 <- pc1 / sd (pc1)
  pc1.quant <- quantile (pc1,
                         probs = c (1,5,10,20,25,30,40,50
                           ,60,70,75,80,90,95,99)/100, names = FALSE)
  pc1.int <- - 1.96 / sqrt (length (pc1))
  pc1.sig <- ifelse (pc1 < pc1.int, 1, 0)
  model <- list ("quantiles" = pc1.quant,
                 "interval" = pc1.int,
                 "code" = pc1.sig)
  output <- cbind (mean.r2, low.r2, up.r2, sd.r2, cmean.r2, csd.r2)
  colnames (output) <- c("ARC","IC-","IC+","SD","CMEAN","CSD")
  rownames (output) <- rownames (vcv1)
  return (list ("out" = output,
                "pc1" = pc1,
                "model" = model,
                "cormat" = cor (t(r2s))))
}

alpha.rep <- function (cor, tam)
  {
    vec <- cor[lower.tri(cor)]
    mvec <- mean(vec)
    varerro <- (1 - (mvec^2))/(tam-2)
    vec2 <- vec^2
    Ex2 <- mean (vec2)
    varvec <- Ex2 - (mean(vec)^2)
    return((varvec - varerro)/varvec)
  }

krz.comp <- function (vcv1, vcv2, d = 19)
  {
    func <- function (x) return (eigen(x)$vectors[,1:d])
    A <- func (vcv1)
    B <- func (vcv2)
    S <- t(A) %*% B %*% t(B) %*% A
    SL <- sum (eigen(S)$values) / d
    return (SL)
  }

mod.pap <- function (vcv, modhip, nit = 1000, neuro.total = FALSE)
  {
    no.hip <- dim (modhip) [2]
    traits <- dim (modhip) [1]
    if (neuro.total)
      add <- 2
    else
      add <- 0
    m.hip.array <- array (0, c(traits, traits, no.hip+add))
    for (i in 1:no.hip)
      {
        m.hip.array[,,i] <- modhip[,i] %*% t(modhip[,i])
        diag (m.hip.array[,,i]) <- 1
        m.hip.array[,,i] <- m.hip.array[,,i] * sqrt (outer (diag(vcv), diag(vcv)))
      }
    if (neuro.total)
      {
        m.hip.array[,,no.hip+1] <-
          apply(!(!(apply (m.hip.array[,,7:8],c(1,2),sum))),2,as.double) ### neuroface
        m.hip.array[,,no.hip+2] <-
          apply(!(!(apply (m.hip.array[,,1:6],c(1,2),sum))),2,as.double) ### total
      }
    output <- list()
    for (N in 1:(no.hip+add))
      {
        tmp <- papyrus (vcv, m.hip.array[,,N], nsk = nit)
        output[[N]] <- tmp
      }
    names(output) <- colnames (modhip)
    if (neuro.total)
      {
        names (output)[no.hip+(1:2)] <- c("neuroface","total")
      }
    return (output)
  }

plot.pap.mod <- function (output, modhip, modwho = "", now = "", plotmat = FALSE)
  {
    if (plotmat == TRUE)
      {
###     layout
        layout (array (c(rep (1, times = 16), 4,
                         rep (3, times = 24), rep (2, times = 41)),
                       c(41,2)))
        par (mar = c(0.2, 5.1, 4.1, 2.1))
      }
    else
      {
        layout (array (c(1,1,2,2),c(2,2)))
        par (mar = c(4.0, 4.0, 4.9, 0.4))
      }
    mean.r2 <- output$out[,1]
    low.r2 <- output$out[,2]
    up.r2 <- output$out[,3]
    c.mean.r2 <- output$out[,5]
    c.sd.r2 <- output$out[,6]
    if (is.null (rownames (output$out)))
      dists <- 1:length (mean.r2)
    else
      dists <- rownames (output$out)
    mod.cex <- ifelse (modhip == 1, 2, 1)
### plot scores
###                b    l    t    r
    plot (mean.r2, type = "p", lty = 2, pch = 17, cex = mod.cex,
          ylab = "", xlab = "", xaxt = "n", ylim = c(-1,1))
    for (i in 1:length (mean.r2))
      {
        abline (v = i, lty = 1, col = rgb (0.8,0.8,0.8))
      }
    arrows (x0 = 1:length (mean.r2),
            y0 = low.r2, y1 = up.r2,
            angle = 90, length = 0.05, code = 3)
    abline (h = mean (mean.r2), lty = 3)
    axis (3, 1:length (mean.r2), dists, las = 2, cex.axis = 0.9)
    mtext (side = 2, at = mean (mean.r2), text = round (mean (mean.r2), 2), las = 2)
### plot av sd
###                b    l    t    r
    par (mar = c(4.0, 0.0, 4.9, 4.6))
    pc.pch <- output$model$code + 17
    plot (c.sd.r2 ~ c.mean.r2, pch = pc.pch, cex = mod.cex, xlab = "", ylab = "",
          yaxt = "n", main = now)
    abline (v = 0, lty = 2)
    abline (h = 0, lty = 2)
    text (c.mean.r2, c.sd.r2, labels = dists, pos = 4, cex = 0.9)
    text (c.mean.r2, c.sd.r2, labels = modwho, pos = 1, cex = 0.9)
    axis (4, las = 2)
    if (plotmat == TRUE)
      {
### matriz
###                                                     b    l    t    r
        par (xaxt = "n", yaxt = "n", cex = 1, mar = c(3.3, 3.4, 0.2, 1.4))
        color2D.matplot (output$cormat,c(-1,1),c(-1,1),c(-1,1), xlab = "", ylab = "")
        par (las = 2, xaxt = "s", yaxt = "s", cex = 0.5)
        axis (side = 1, 1:length (dists), dists)
        axis (side = 2, (1:length (dists))[length (dists):1], dists)
### barra
        .seq <- array (seq (from = (-1), to = 1, by = 0.1),c (1,21))
###                                                     b    l    t    r
        par (xaxt = "n", yaxt = "n", cex = 1, mar = c(0.2, 3.4, 0.2, 1.4))
        color2D.matplot (.seq, c(-1,1),c(-1,1),c(-1,1), xlab = "", ylab = "")
        text (x = 0.5:20.5, y = 0.5, seq (from = (-1), to = 1, by = 0.1), cex = 0.4, col = rgb (21:0,21:0,21:0, maxColorValue = 21))
        par (xaxt = "s", yaxt = "s")
        return (0)
      }
  }

plot.pap.2 <- function (output, modhip, modwho = "", now = "")
  {
    par (mar = c(4.0, 4.0, 4.9, 0.4))
    mean.r2 <- output$out[,1]
    low.r2 <- output$out[,2]
    up.r2 <- output$out[,3]
    c.mean.r2 <- output$out[,5]
    c.sd.r2 <- output$out[,6]
    if (is.null (rownames (output$out)))
      dists <- 1:length (mean.r2)
    else
      dists <- rownames (output$out)
    mod.cex <- ifelse (modhip == 1, 2, 1)
### plot scores
###                b    l    t    r
    plot (mean.r2, type = "p", lty = 2, pch = 20, cex = mod.cex,
          ylab = "", xlab = "", xaxt = "n", ylim = c(-1,1))
    for (i in 1:length (mean.r2))
      {
        abline (v = i, lty = 1, col = rgb (0.8,0.8,0.8))
      }
    arrows (x0 = 1:length (mean.r2),
            y0 = low.r2, y1 = up.r2,
            angle = 90, length = 0.05, code = 3)
    abline (h = mean (mean.r2), lty = 3)
    axis (3, 1:length (mean.r2), dists, las = 2, cex.axis = 0.9)
    mtext (side = 2, at = 0.8, text = round (mean (mean.r2), 2), las = 2)
### plot av sd
###                b    l    t    r
    par (mar = c(4.0, 0.0, 4.9, 4.6))
    pc.pch <- output$model$code + 20
    plot (c.sd.r2 ~ c.mean.r2, pch = pc.pch, cex = mod.cex, xlab = "", ylab = "",
          yaxt = "n", main = now)
    abline (v = 0, lty = 2)
    abline (h = 0, lty = 2)
    text (c.mean.r2, c.sd.r2, labels = dists, pos = 4, cex = 0.9)
    text (c.mean.r2, c.sd.r2, labels = modwho, pos = 3, cex = 0.9)
    axis (4, las = 2)
  }

rep.nova = function (ID, data.matrix)
  {
    ### Lessels, C. M., & Boag, P. T. (1987).
    ### Unrepeatable repeatabilities: a common mistake.
    ### The Auk, 2(January), 116–121.
    chars = ncol (data.matrix)
    model.gen = function (vec) return (lm (vec ~ ID))
    models.list = apply (data.matrix, 2, model.gen)
    models.list = lapply (models.list, anova)
    rep.itself = function (summ)
      {
        msq = summ$'Mean Sq' ## 1 entre, 2 dentro
        s2a = (msq[1] - msq[2])/2
        out = s2a / (s2a + msq[2])
        return (out)
      }
    out = sapply (models.list, rep.itself)
    names (out) = colnames (data.matrix)
    return (out)
  }
