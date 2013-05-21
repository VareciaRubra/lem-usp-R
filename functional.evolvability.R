FuncEvol = function (cov.matrix, func.betas)
  {
    # input:
    #  cov.matrix: a vcv matrix
    #  func.beta: set of modularity hypothesis
    # returns:
    #  evolvabilities of functional betas
    func.betas = apply (func.betas, 2, Normalize)
    evol = apply (func.betas, 2, Evolvability, cov.matrix = cov.matrix)
    return (evol)
  }

RandomCov = function (cov.matrix)
  {
    # input:
    #  cov.matrix: a vcv matrix
    # returns:
    #  a vcv matrix with the same eigenvalue distribution from cov.matrix
    #  but with random eigenvectors
    eig.dec = eigen (cov.matrix)
    random.evec = array (rnorm (prod (dim (cov.matrix))), dim (cov.matrix))
    random.evec = qr.Q (qr (random.evec))
    cov.matrix.revec = random.evec %*% diag (eig.dec$values) %*% t (random.evec)
    return (cov.matrix.revec)
  }

JointModTest = function (cov.matrix, hips, n.randmat = 1000)
  {
    explained.or = sum (FuncEvol (cov.matrix, func.betas = hips))
    total = sum (diag (cov.matrix))
    explained.random = c()
    for (i in 1:n.randmat)
      {
        randmat = RandomCov (cov.matrix)
        explained.random [i] = sum (FuncEvol (randmat, func.betas = hips))
      }
    p = sum (explained.random > explained.or) / n.randmat
    return (c(explained.or, p))
  }

JointModTest2 = function (cov.matrix, hips, n.randhip = 1000)
    {
      explained.or = sum (FuncEvol (cov.matrix, func.betas = hips))
      total = sum (diag (cov.matrix))
      explained.random = c()
      for (i in 1:n.randhip)
        {
          randhip = RandomHip (hips)
          explained.random [i] = sum (FuncEvol (cov.matrix, func.betas = randhip))
        }
      p = sum (explained.random > explained.or) / n.randhip
      return (c(explained.or, p))
    }

RandomHip = function (hip)
  {
    # input:
    #  hip: set of MI hypothesis, coded as functional vectors
    # returns:
    #  random hypothesis with the same number of modules
    out = NULL
    out = t (apply (hip, 1, sample, replace = FALSE))
    dimnames (out) = dimnames (hip)
    return (out)
  }
