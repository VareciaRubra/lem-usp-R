Norm <- function(x){return(sqrt(sum(x*x)))}
Normalize <- function(x){return(x/Norm(x))}

RemoveSize <- function (cov.matrix)
  # Removes first principal component effect in cov.matrix.
  #
  # Args:
  #   cov.matrix: A simetric covariance matrix
  # Return:
  #   cov.matrix with size removed
{
  cov.matrix.svd  <-  svd(cov.matrix)
  size.eigen.vector <- cov.matrix.svd$u[, 1]
  size.eigen.value <- cov.matrix.svd$d[1]
  size.factor <- size.eigen.vector * sqrt(size.eigen.value)
  cov.matrix.size.removed <- cov.matrix - size.factor %*% t(size.factor)
  return (cov.matrix.size.removed)
}

MonteCarloR2 <- function (corr.matrix, sample.size, iterations = 1000)
  # Computes a distribution of magnitudes of integration (r2)
  # for a given correlation matrix.
  #
  # Args:
  #   corr.matrix: a square symmetric correlation matrix
  #   sample.size: number of individuals to sample
  #   iterations: number of populations sampled
  # Return:
  #   a vector with iterations + 1 entries; the first entry corresponds
  #   to the actual r2 value calculated from the original matrix. The
  #   remaining entries are values calculated from the resampling procedure.
{
  R2 <- function (matrix)
    return (mean (matrix [lower.tri (matrix)]^2))
  n.traits <- dim (corr.matrix) [1]
  populations <- list ()
  for (i in 1:iterations)
    populations [[i]] <- rmvnorm2 (sample.size, sigma = corr.matrix, method = 'chol')
  it.matrices <- lapply (populations, cor)
  it.r2 <- sapply (it.matrices, R2)
  r2 <- c (R2 (corr.matrix), it.r2)
  return (r2)
}

TestModularity <- function (cor.matrix, modularity.hipot)
  # Tests modularity hipotesis using cor.matrix matrix and trait groupings
  #
  # Args:
  #   cor.matrix: Correlation matrix
  #   modularity.hipot: Matrix of hipotesis.
  #                     Each line represents a trait and each column a module.
  #                     if modularity.hipot[i,j] == 1, trait i is in module j.
  # Return:
  #   list with mantel correlation of cor.matrix with binary hipotesis matrices
{
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

MultiRsMantel <- function (matrix.list, MatrixCompFunc = RandomSkewers, repeat.vector = NULL, iterations = 1000)
  # Performs multiple comparisons between a set of covariance of correlation matrices.
  #
  # Args:
  #  matrix.list: a list of covariance or correlation matrices
  #  MatrixCompFunc: function to use for comparison
  #  repeat.vector: vector of matrix repeatabilities
  #  iterations: number of RandomSkewers or matrix permutations passed to MatrixCompFunc
  #
  # Return:
  #  a list with two matrices containing $\Gamma$-values or average random
  #  skewers correlation and probabilities according to permutation test.
  #  if repeat.vector was also passed, values below the diagonal on the correlation matrix
  #  will contain corrected correlation values.
{
  n.matrix <- length (matrix.list)
  matrix.names <- names (matrix.list)
  probabilities <- array (0, c(n.matrix, n.matrix))
  correlations <- probabilities
  for (i in 1:(n.matrix - 1)) {
    for (j in (i+1):n.matrix) {
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
    diag (correlations) <- repeat.vector
  }
  rownames (correlations) <- matrix.names
  colnames (correlations) <- matrix.names
  dimnames (probabilities) <- dimnames (correlations)
  output <- list ('correlations' = correlations, 'probabilities' = probabilities)
  return (output)
}

rmvNorm2 <- function (n, theta = rep(0, nrow(sigma)),
                      sigma = diag(length(theta)),
                      method = c("eigen", "svd", "chol"))
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
{
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

RandomSkewers <- function (cov.matrix.1, cov.matrix.2, nsk = 10000)
  # Calculates covariance matrix correlation via random skewers
  # Args:
  #     cov.matrix.(1,2): Two covariance matrices to be compared
  #     nsk: Number of generated random skewers
  # Return:
  #     List with mean value of correlation, p value and standard deviation
{
  traits <- dim (cov.matrix.1) [1]
  base.vector <- Normalize(rnorm(traits))
  random.vectors <- array (rnorm (nsk * traits, mean = 0, sd = 1), c(traits, nsk))
  random.vectors <- apply (random.vectors, 2, Normalize)
  dist <- base.vector %*% random.vectors
  dz1 <- apply (cov.matrix.1 %*% random.vectors, 2, Normalize)
  dz2 <- apply (cov.matrix.2 %*% random.vectors, 2, Normalize)
  real <- apply (dz1 * dz2, 2, sum)
  ac <- mean (real)
  stdev <- sd (real)
  prob <- sum (ac < dist) / nsk
  output <- c(ac, prob, stdev)
  names(output) <- c("AC","P","SD")
  return(output)
}

MantelCor <- function (cor.matrix.1, cor.matrix.2, nit = 1000, mod = FALSE)
  # Calculates matrix correlation with confidence intervals using mantel permutations
  #
  # Args:
  #     cor.matrix.(1,2): correlation matrices being compared
  #     nit: number of permutations
  #     mod: for when testing binary modularity hipotesis
  # Return:
  #     matrix pearson correelation and significance.
  #     if mod==TRUE also returns average within, between and average ratio correlations
{
  fixed.matrix <- cor.matrix.1 [lower.tri (cor.matrix.1)]
  shuffled.matrix <- cor.matrix.2 [lower.tri (cor.matrix.2)]
  correlation <- cor (fixed.matrix,shuffled.matrix)
  shuffled.correlation <- c()
  for (N in 1:nit){
    shuffle <- sample (1:dim(cor.matrix.1)[1])
    shuffled.matrix <- cor.matrix.2 [shuffle, shuffle] [lower.tri (cor.matrix.1)]
    shuffled.correlation[N] <- cor (fixed.matrix,shuffled.matrix)
  }
  prob <- sum (shuffled.correlation > as.vector(correlation)) / nit
  if (mod == TRUE){
    avg.plus <- mean (cor.matrix.1 [lower.tri(cor.matrix.1)] [mob != 0])
    avg.minus <- mean (cor.matrix.1 [lower.tri(cor.matrix.1)] [mob == 0])
    avg.ratio <- avg.plus / avg.minus
    output <- c(correlation,prob,avg.plus,avg.minus,avg.ratio)
    names(output) <- c("R²","Probability","AVG+","AVG-","AVG Ratio")
  }
  else{
    output <- c(correlation,prob)
    names(output) <- c("R²","Probability")
  }
  return (output)
}

SRD <- function (cov.matrix.1, cov.matrix.2, nsk = 1000)
  # Calculates the selection response decomposition comparison between covariance matrices
  #
  # Args:
  #     cov.matrix.(1,2): covariance matrices being compared
  #     nsk: number of RandomSkewers random vectors
  # Return:
  #     SRD scores for each trait and significance using mean and sd of SRD scores
{
  size <- dim (cov.matrix.1)[1]
  r2s <- array (0, c(size,nsk))
  beta <- apply (array (rnorm (size*nsk, mean = 0, sd = 1),c(size,nsk)),2, Normalize)
  for (I in 1:nsk){
    beta.matrix <- diag (beta[,I])
    dz1 <- apply (cov.matrix.1 %*% beta.matrix, 1, Normalize)
    dz2 <- apply (cov.matrix.2 %*% beta.matrix, 1, Normalize)
    r2s[,I] <- colSums (dz1 * dz2)
  }
  # results
  mean.r2 <- apply (r2s, 1, mean)
  sample.conf <- function (x, lower = TRUE){
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
                         probs = c (1,5,10,20,25,30,40,50,60,70,75,80,90,95,99)/100,
                         names = FALSE)
  pc1.int <- - 1.96 / sqrt (length (pc1))
  pc1.sig <- ifelse (pc1 < pc1.int, 1, 0)
  model <- list ("quantiles" = pc1.quant,
                 "interval"  = pc1.int,
                 "code"      = pc1.sig)
  output <- cbind (mean.r2, low.r2, up.r2, sd.r2, cmean.r2, csd.r2)
  colnames (output) <- c("ARC","IC-","IC+","SD","CMEAN","CSD")
  rownames (output) <- rownames (cov.matrix.1)
  return (list ("out"    = output,
                "pc1"    = pc1,
                "model"  = model,
                "cormat" = cor (t(r2s))))
}

PlotSRD <- function (output, matrix.label = "")
  # Plots the output of the SRD function in standard format
  #
  # Args:
  #     output: the output from the SRD funtion
  #     matrix.label: string with the names of the matrices that were compared in the SRD function
  # Return:
  #     pretty plot
{
  layout (array (c(1,1,2,2),c(2,2)))
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
  ### plot scores
  ###                b    l    t    r
  plot (mean.r2, type = "p", lty = 2, pch = 17,
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
  plot (c.sd.r2 ~ c.mean.r2, pch = pc.pch, xlab = "", ylab = "",
        yaxt = "n", main = matrix.label)
  abline (v = 0, lty = 2)
  abline (h = 0, lty = 2)
  text (c.mean.r2, c.sd.r2, labels = dists, pos = 4, cex = 0.9)
  axis (4, las = 2)
}

KzrCor <- function (cov.matrix.1, cov.matrix.2, ret.dim = NULL)
  # Calculates the Kzranowski correlation between matrices
  #
  # Args:
  #     cov.matrix.(1,2): covariance being compared
  #     ret.dim: number of retained dimensions in the comparison,
  #              default for nxn matrix is n/2-1
  # Return:
  #     Kzranowski correlation
{
  if (is.null(ret.dim))
    ret.dim = dim(cov.matrix.1)[1]/2 - 1
  EigenVectors <- function (x) return (eigen(x)$vectors[,1:ret.dim])
  A <- EigenVectors (cov.matrix.1)
  B <- EigenVectors (cov.matrix.2)
  S <- t(A) %*% B %*% t(B) %*% A
  SL <- sum (eigen(S)$values) / ret.dim
  return (SL)
}

CalcRepeatability <- function (ID, ind.data)
  # Calculates Repeatabilities acording to:
  #    Lessels, C. M., & Boag, P. T. (1987).
  #    Unrepeatable repeatabilities: a common mistake.
  #    The Auk, 2(January), 116–121.
  # Args:
  #     ID: indentity of individuals
  #     ind.data: individual measurments
  # Return:
  #     vector of repeatabilities
{
  models.list = apply (ind.data, 2, function (vec){return (lm (vec ~ ID))})
  models.list = lapply (models.list, anova)
  rep.itself = function (lm.model){
    msq = lm.model$'Mean Sq' ## 1 entre, 2 dentro
    s2a = (msq[1] - msq[2])/2
    out = s2a / (s2a + msq[2])
    return (out)
  }
  out = sapply (models.list, rep.itself)
  names (out) = colnames (ind.data)
  return (out)
}

AlphaRep <- function (cor.matrix, tam)
  # Calculates the matrix repeatability using the equation in Cheverud 1996
  # Quantitative genetic analysis of cranial morphology in the cotton-top
  # (Saguinus oedipus) and saddle-back (S. fuscicollis) tamarins. Journal of Evolutionary Biology 9, 5-42.
  #
  # Args:
  #     cor.matrix: correlation matrix
  #     tam: sample size
  # Return:
  #     matrix repeatability
{
  vec <- cor.matrix[lower.tri(cor.matrix)]
  var.erro <- (1 - mean(vec)^2)/(tam-2)
  var.vec <- var(vec)
  return((var.vec - var.erro)/var.vec)
}

BootstrapRep <- function (ind.data, nb = 100)
  # Calculates the repeatability of the covariance matrix of the suplied data
  # via bootstrap ressampling
  #
  # Args:
  #     ind.data: original individual data
  #     nb = number of resamples
  # Return:
  #     returns the mean repeatability
{
  n.ind <-  dim (ind.data) [1]
  original.cov.matrix <- var (ind.data)
  v.rep <- c()
  for (N in 1:nb){
    sampled.data <- sample (1:n.ind, n.ind, TRUE)
    sampled.data.cov.matrix <- var (ind.data[sampled.data,])
    v.rep [N] <- RandomSkewers (original.cov.matrix, sampled.data.cov.matrix, 1000) [1]
  }
  out <- mean (v.rep)
  return (out)
}

MonteCarloRep <- function (x.matrix, ind, nit = 100)
  # Calculates x.matrix repeatability using parametric sampling
  #
  # Args:
  #     x.matrix: covariance or correlation matrix. 
  #               if x.matrix is a correlation matrix will use MantelCor, 
  #               else, will use RandomSkewers
  #     ind: number of indivuals on each sample
  #     nit: number of samples
  # Return:
  #     mean correlation of sample covariance matrices with original input x.matrix
{
  if (sum(diag(x.matrix)) == dim (x.matrix) [1]){
    Func <- MantelCor
    Type <- cor
  }
  else{
    Func <- RandomSkewers
    Type <- var
  }
  R <- c()
  for (N in 1:nit){
    rand.samp <- rmvNorm2 (ind, rep(0, times = dim (x.matrix)[1]),
                           sigma = x.matrix, method = "chol")
    rand.matrix <- Type (rand.samp)
    R[N] <- Func (x.matrix, rand.matrix, 1000)[1]
  }
  return (mean(R))
}

ExtendMatrix <- function(cov.matrix, cutoff = NULL)
  # Calculates the noise controled covariance matrix using the extension method
  #
  # Args:
  #     cov.matrix: covariance matrixe being extended.
  #                 must be larger then 10x10
  #     cutoff: number of retained eigenvalues
  #             if is not supplied will be calculated using the gradient variance method
  # Return:
  #     returns the extended covariance matrix
{
  p = dim(cov.matrix)[1]
  if(dim(p<10))
    stop("matrix is too small")
  eigen.cov.matrix = eigen(cov.matrix)
  eVal = eigen.cov.matrix$values
  eVec = eigen.cov.matrix$vectors
  if(is.null(cutoff)){
    grad = array(dim=c(p-2))
    tr.cov.matrix = sum(eVal)
    for (i in 1:(p-2))
      grad[i] = abs(eVal[i]/tr.cov.matrix - 2*(eVal[i+1]/tr.cov.matrix) + eVal[i+2]/tr.cov.matrix)
    var.grad = array(dim=c(p-6))
    for(i in 1:(p-6)){
      var.grad[i] = var(grad[i:(i+4)])
    }
    length(var.grad[var.grad<1e-4])
    x11()
    plot(4:(p-3),var.grad)
    cutoff = floor(locator(1)$x)
  }
  eVal[eVal < eVal[cutoff]] = eVal[cutoff]
  extended.cov.matrix = eVec%*%diag(eVal)%*%t(eVec)
  colnames(extended.cov.matrix) = colnames(cov.matrix)
  rownames(extended.cov.matrix) = rownames(cov.matrix)
  return(extended.cov.matrix)
}

### Hansen, T. F., & Houle, D. (2008).
### Measuring and comparing evolvability and constraint in multivariate characters.
### Journal of Evolutionary Biology, 21(5), 1201–1219.
### doi:10.1111/j.1420-9101.2008.01573.x
Respondability          <- function (beta, cov.matrix) return (Norm (cov.matrix %*% beta))
Evolvability            <- function (beta, cov.matrix) return (t (beta) %*% cov.matrix %*% beta)
ConditionalEvolvability <- function (beta, cov.matrix) return ((t (beta) %*% solve (cov.matrix) %*% beta)^(-1))
Autonomy                <- function (beta, cov.matrix) return (((t (beta) %*% solve (cov.matrix) %*% beta)^(-1)) / (t (beta) %*% cov.matrix %*% beta))
Flexibility             <- function (beta, cov.matrix) return (t (beta) %*% cov.matrix %*% beta / Norm (cov.matrix %*% beta))
Constraints             <- function (beta, cov.matrix) return (abs (t (Normalize (eigen (cov.matrix)$vectors[,1])) %*% Normalize (cov.matrix %*% beta)))
### Marroig, G., Shirai, L. T., Porto, A., de Oliveira, F., & de Conto, V. (2009).
### The evolution of modularity in the mammalian skull II: evolutionary consequences.
### Evolutionary Biology, 36(1), 136–148.
### doi:10.1007/s11692-009-9051-1
MeanSquaredCorrelation <- function (cov.matrix) return (mean (cov2cor (cov.matrix) [lower.tri (diag (nrow (cov.matrix)))]^2))
Pc1Percent             <- function (cov.matrix) return (eigen (cov.matrix)$values [1] / sum (eigen (cov.matrix)$values))

load.of.functions <- list (Norm,
                          Normalize,
                          hansen.houle = list (Respondability,
                                               Evolvability,
                                               ConditionalEvolvability,
                                               Autonomy,
                                               Flexibility,
                                               Constraints),
                          MeanSquaredCorrelation,
                          Pc1Percent)

HansenHouleAverage <- function (mat, nsk = 10000)
{
  with (load.of.functions,
        {
          n.char <- dim (mat) [1]
          beta.mat <- array (rnorm (n.char * nsk), c(n.char, nsk))
          beta.mat <- apply (beta.mat, 2, Normalize)
          iso.vec <- Normalize (rep(1, times = n.char))
          null.dist <- abs (t (iso.vec) %*% beta.mat)
          null.dist <- sort (null.dist)
          crit.value <- null.dist [round (0.95 * nsk)]
          cat ('critical value: ', crit.value, '\n')
          parm.dist <- array (0, c(nsk, 8))
          HansenHouleWrap <- function (hh.func) return (apply (beta.mat, 2, hh.func, cov.matrix = mat))
          parm.dist [,1:6] <- sapply (hansen.houle, HansenHouleWrap)
          parm.dist[,7] <- as.numeric (parm.dist[,5] > crit.value)
          parm.dist[,8] <- as.numeric (parm.dist[,6] > crit.value)
          parm.dist <- cbind (parm.dist, null.dist)
          colnames (parm.dist) <- c('resp','evol','cond.evol', 'auto',
                                    'flex','const','flex.n', 'const.n', 'null.dist')
          parm.av <- colMeans (parm.dist)
          parm.av[7:8] <- parm.av[7:8] * nsk
          pc1 <- eigen (mat)$vectors[,1]
          HansenHouleWrapPc1 <- function (hh.func) return (hh.func (beta = pc1, cov.matrix = mat))
          maximum <- sapply (hansen.houle, HansenHouleWrapPc1)
          integration <- c (MeanSquaredCorrelation (mat), Pc1Percent (mat))
          names (integration) <- c ('MeanSquaredCorrelation', 'pc1%')
          parm.av <- c (integration, parm.av)
          return (list ('dist' = parm.dist, 'mean' = parm.av, 'max.val' = maximum))
        })
}

hh.mod <- function (mat, hip, nsk = 10000)
{
  with (load.of.functions,
        {
          out <- HansenHouleAverage (mat, nsk)
          HansenHouleWrap2 <- function (hh.func) return (apply (hip, 2, hh.func, cov.matrix = mat))
          hip <- apply (hip, 2, Normalize)
          out$mod <- sapply (hansen.houle, HansenHouleWrap2)
          return (out)
        })
}

plot.mod.evol <- function (evo.out, new.dev = TRUE)
{
  require (MASS)
  with (evo.out,
        {
          n.hip <- nrow (mod)
          n.sk <- nrow (dist)
          if (new.dev)
            par (mfrow = c(2,3))
          out <- array (0, dim (mod))
          for (j in 1:6)
          {
            dist[,j] <- dist [,j] / max.val [j]
            mod [,j] <- mod [,j] / max.val [j]
            x.lim <- c (ifelse (min (mod[,j]) < min (dist[,j]), min (mod[,j]), min (dist[,j])),
                       ifelse (max (mod[,j]) > max (dist[,j]), max (mod[,j]), max (dist[,j])))
            truehist (dist[,j], prob = TRUE, border = 'grey', col = 'white', xlim = x.lim,
                      xlab = colnames (mod)[j], main = '')
            lines (density (dist[,j]))
            abline (v = sort (dist[,j]) [round (0.975 * n.sk)], lty = 2, col = 'red')
            abline (v = sort (dist[,j]) [round (0.025 * n.sk)], lty = 2, col = 'red')
            for (i in 1:n.hip)
            {
              abline (v = mod[i,j], lty = 3, col = rgb (.2,.2,.2))
              mtext (side = 3, at = mod[i,j], text = rownames (mod)[i], las = 3, cex = 0.7)
              out [i,j] <- sum (mod [i,j] < dist[,j])/n.sk
            }
          }
          dimnames (out) <- dimnames (mod)
          return (out)
        })
}

plot.group.evol <- function (evol.list)
{
  upper.crit <- function (vec, val = .95)
  {
    n.obs <- length (vec)
    vec <- sort (vec)
    crit.pos <- round (val * n.obs)
    return (vec[crit.pos])
  }
  plot.single <- function (evol.list, which)
  {
    n.obj <- length (evol.list)
    extract.dist <- function (element, which) return (element$dist[,which]/element$max.val[which])
    to.plot <- sapply (evol.list, extract.dist, which)
    colnames (to.plot) <- names (evol.list)
    boxplot (to.plot, cex = 0.2, boxwex = 0.5, border = rgb (.6,.6,.6))
    extract.mod <- function (element, which) return (element$mod[,which]/element$max.val[which])
    extract.mod.names <- function (element) return (rownames (element$mod))
    mod.values <- lapply (evol.list, extract.mod, which)
    mod.names <-  lapply (evol.list, extract.mod.names)
    crits <- apply (to.plot, 2, upper.crit)
    for (i in 1:n.obj)
    {
      segments (x0 = i - 0.25, x1 = i + 0.25, y0 = crits [i],
                col = 'red', lwd = 2)
      segments (x0 = i - 0.25, x1 = i + 0.25, y0 = mod.values [[i]],
                col = 'black', lty = i)
    }
  }
}

ResidualMatrix <- function (model)
  {
    ## Calculates residual matrix from a estimated linear model
    ## Args:
    ##  model: linear model previously estimated
    ## Value:
    ##  cov.mat: residual covariance matrix
    res <- residuals (model)
    res.df <- model $ df.residual
    cov.mat <- t (res) %*% res / res.df
    return (cov.mat)
  }
