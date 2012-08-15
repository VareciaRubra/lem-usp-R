load.of.functions = list (
### Basic Linear Algebra
  norma = function (x) return (sqrt (sum (x * x))),
  normalize = function (x) return (x / norma (x)),
### Hansen, T. F., & Houle, D. (2008).
### Measuring and comparing evolvability and constraint in multivariate characters.
### Journal of Evolutionary Biology, 21(5), 1201–1219.
### doi:10.1111/j.1420-9101.2008.01573.x
  hansen.houle = list (
    Respondability = function (beta, C) return (norma (C %*% beta)),
    Evolvability = function (beta, C) return (t (beta) %*% C %*% beta),
    Conditional.Evolvability = function (beta, C) return ((t (beta) %*% solve (C) %*% beta)^(-1)),
    Autonomy = function (beta, C) return (((t (beta) %*% solve (C) %*% beta)^(-1)) / (t (beta) %*% C %*% beta)),
    Flexibility = function (beta, C) return (t (beta) %*% C %*% beta / norma (C %*% beta)),
    Constraints = function (beta, C) return (abs (t (normalize (eigen (C)$vectors[,1])) %*% normalize (C %*% beta)))),
### Marroig, G., Shirai, L. T., Porto, A., de Oliveira, F., & de Conto, V. (2009).
### The evolution of modularity in the mammalian skull II: evolutionary consequences.
### Evolutionary Biology, 36(1), 136–148.
### doi:10.1007/s11692-009-9051-1
  r2 = function (C) return (mean (cov2cor (C) [lower.tri (diag (nrow (C)))]^2)),
  pc1.percent = function (C) return (eigen (C)$values [1] / sum (eigen (C)$values))
  )

hh.average = function (mat, nsk = 10000)
  {
    with (load.of.functions,
          {
            n.char = dim (mat) [1]
            beta.mat = array (rnorm (n.char * nsk), c(n.char, nsk))
            beta.mat = apply (beta.mat, 2, normalize)
            iso.vec = normalize (rep(1, times = n.char))
            null.dist = abs (t (iso.vec) %*% beta.mat)
            null.dist = sort (null.dist)
            crit.value = null.dist [round (0.95 * nsk)]
            cat ('critical value: ', crit.value, '\n')
            parm.dist = array (0, c(nsk, 8))
            hh.wrap = function (hh.func) return (apply (beta.mat, 2, hh.func, C = mat))
            parm.dist [,1:6] = sapply (hansen.houle, hh.wrap)
            parm.dist[,7] = as.numeric (parm.dist[,5] > crit.value)
            parm.dist[,8] = as.numeric (parm.dist[,6] > crit.value)
            parm.dist = cbind (parm.dist, null.dist)
            colnames (parm.dist) = c('resp','evol','cond.evol', 'auto',
                       'flex','const','flex.n', 'const.n', 'null.dist')
            parm.av = colMeans (parm.dist)
            parm.av[7:8] = parm.av[7:8] * nsk
            pc1 = eigen (mat)$vectors[,1]
            hh.wrap.pc1 = function (hh.func) return (hh.func (beta = pc1, C = mat))
            maximum = sapply (hansen.houle, hh.wrap.pc1)
            integration = c (r2 (mat), pc1.percent (mat))
            names (integration) = c ('r2', 'pc1%')
            parm.av = c (integration, parm.av)
            return (list ('dist' = parm.dist, 'mean' = parm.av, 'max.val' = maximum))
          })
  }

hh.mod = function (mat, hip, nsk = 10000)
  {
    with (load.of.functions,
          {
            out = hh.average (mat, nsk)
            hh.wrap2 = function (hh.func) return (apply (hip, 2, hh.func, C = mat))
            hip = apply (hip, 2, normalize)
            out$mod = sapply (hansen.houle, hh.wrap2)
            return (out)
          })
  }

plot.mod.evol = function (evo.out, new.dev = TRUE)
  {
    require (MASS)
    with (evo.out,
          {
            n.hip = nrow (mod)
            n.sk = nrow (dist)
            if (new.dev)
              par (mfrow = c(2,3))
            out = array (0, dim (mod))
            for (j in 1:6)
              {
                dist[,j] = dist [,j] / max.val [j]
                mod [,j] = mod [,j] / max.val [j]
                x.lim = c (ifelse (min (mod[,j]) < min (dist[,j]), min (mod[,j]), min (dist[,j])),
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
                    out [i,j] = sum (mod [i,j] < dist[,j])/n.sk
                  }
              }
            dimnames (out) = dimnames (mod)
            return (out)
          })
  }
