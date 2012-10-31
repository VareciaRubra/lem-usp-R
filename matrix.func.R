Norm = function(x){return(sqrt(sum(x*x)))}
Normalize = function(x){return(x/Norm(x))}

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

MultiRsMantel <- function (matrix.list, MatrixCompFunc = RandomSkewers, repeat.vector = NULL, iterations = 1000){
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
  #  if repeat.vector was also passed, values above the diagonal on the correlation matrix
  #  will contain corrected correlation values.
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

RandomSkewers <- function (cov.matrix.1, cov.matrix.2, nsk = 10000){
  # Calculates covariance matrix correlation via random skewers
  # Args:
  #     cov.matrix.(1,2): Two covariance matrices to be compared
  #     nsk: Number of generated random skewers
  # Return:
  #     List with mean value of correlation, p value and standard deviation
  traits <- dim (cov.matrix.1) [1]
  base.vector <- Normalize(rnorm(traits))
  random.vectors <- array (rnorm (nsk * traits, mean = 0, sd = 1), c(traits, nsk))
  random.vectors <- apply (random.vectors, 2, Normalize)
  dist <- abs (base.vector %*% random.vectors)
  dz1 <- apply (cov.matrix.1 %*% random.vectors, 2, Normalize)
  dz2 <- apply (cov.matrix.2 %*% random.vectors, 2, Normalize)
  real <- abs (apply (dz1 * dz2, 2, sum))
  ac <- mean (real)
  stdev <- sd (real)
  prob <- sum (real < dist) / nsk
  output <- c(ac, prob, stdev)
  names(output) <- c("AC","P","SD")
  return(output)
}

MantelCor <- function (cor.matrix.1, cor.matrix.2, nit = 1000, mod = FALSE){
  # Calculates matrix correlation with confidence intervals using mantel permutations
  #
  # Args:
  #     cor.matrix.(1,2): correlation matrices being compared
  #     nit: number of permutations
  #     mod: for when testing binary modularity hipotesis
  # Return:
  #     matrix pearson correelation and significance.
  #     if mod==TRUE also returns average within, between and average ratio correlations
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

SRD <- function (cov.matrix.1, cov.matrix.2, nsk = 1000){
  # Calculates the selection response decomposition comparison between covariance matrices
  #
  # Args:
  #     cov.matrix.(1,2): covariance matrices being compared
  #     nsk: number of RandomSkewers random vectors
  # Return:
  #     SRD scores for each trait and significance using mean and sd of SRD scores
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

PlotSRD <- function (output, matrix.label = ""){
  # Plots the output of the SRD function in standard format
  #
  # Args:
  #     output: the output from the SRD funtion
  #     matrix.label: string with the names of the matrices that were compared in the SRD function
  # Return:
  #     pretty plot
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

KzrCor <- function (cov.matrix.1, cov.matrix.2, ret.dim = 19){
  # Calculates the Kzranowski correlation between matrices
  #
  # Args:
  #     cov.matrix.(1,2): covariance being compared
  #     ret.dim: number of retained dimensions in the comparison
  # Return:
  #     Kzranowski correlation
  func <- function (x) return (eigen(x)$vectors[,1:d])
  A <- func (cov.matrix.1)
  B <- func (cov.matrix.2)
  S <- t(A) %*% B %*% t(B) %*% A
  SL <- sum (eigen(S)$values) / d
  return (SL)
}

CalcRepeatability <- function (ID, ind.data){
  # Calculates Repeatabilities acording to:
  #    Lessels, C. M., & Boag, P. T. (1987).
  #    Unrepeatable repeatabilities: a common mistake.
  #    The Auk, 2(January), 116–121.
  # Args:
  #     ID: indentity of individuals
  #     ind.data: individual measurments
  # Return:
  #     vector of repeatabilities
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

AlphaRep <- function (cor.matrix, tam) {
  # Calculates the matrix repeatability using the equation in Cheverud 1996
  # Quantitative genetic analysis of cranial morphology in the cotton-top
  # (Saguinus oedipus) and saddle-back (S. fuscicollis) tamarins. Journal of Evolutionary Biology 9, 5-42.
  #
  # Args:
  #     cor.matrix: correlation matrix
  #     tam: sample size
  # Return:
  #     matrix repeatability
  vec <- cor.matrix[lower.tri(cor.matrix)]
  var.erro <- (1 - mean(vec)^2)/(tam-2)
  var.vec <- var(vec)
  return((var.vec - var.erro)/var.vec)
}

BootstrapRep <- function (ind.data, nb = 100){
  # Calculates the repeatability of the covariance matrix of the suplied data
  # via bootstrap ressampling
  #
  # Args:
  #     ind.data: original individual data
  #     nb = number of resamples
  # Return:
  #     returns the mean repeatability
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

MonteCarloRep <- function (x.matrix, ind, nit = 100){
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

ExtendMatrix <- function(cov.matrix, cutoff = NULL){
  # Calculates the noise controled covariance matrix using the extension method
  #
  # Args:
  #     cov.matrix: covariance matrixe being extended.
  #                 must be larger then 10x10
  #     cutoff: number of retained eigen values
  #             if is not supplied will be calculated using the gradient variance method
  # Return:
  #     returns the exetended convariance matrix
  if(dim(cov.matrix)[1]<10)
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
