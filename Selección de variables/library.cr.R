library(survival)
library(compiler)
library(LearnBayes)
library(plyr)
library(numDeriv)
library(doParallel)
library(doSNOW)


# cálculo de m_{\gamma}(y).
lmarg1.laplace.SigmaM.normal.cr = function(y, rel, X1, X2, Weig.matrix1, Weig.matrix2, ne, use.cox = TRUE) {
  # inicializo todo antes de llamar a la función.
  Z <- NULL
  p11 <- NULL
  n <- NULL
  Weig.matrix <- NULL
  SigmaM.hat <- vector("list", 2)
  tt1 <- list()
  Xt.Weig.X <- NULL
  parms.hat1 <- vector("list", 2)
  hessi.prior <- vector("list", 2)
  X_matrix <- NULL
  for (k in 1:2) {
    if (k == 1) {
    X_matrix <- X1
    Weig.matrix <- Weig.matrix1
    } else {
      X_matrix <- X2
      Weig.matrix <- Weig.matrix2
    }
    Z = as.matrix(X_matrix) # matriz
    p11[k] = ncol(Z) # número de columnas.
    n = nrow(Z) # número de filas
    if (p11[k] != 0) { # si hay activas variables para el riesgo k
      Xt.Weig.X = t(Z) %*% Weig.matrix %*% Z
      SigmaM.hat[[k]] = ne[k] * solve(Xt.Weig.X)
      # Hessiano de la log prior, es -Xt.Weig.X/ne, cambiamos de signo porque necesitamos el hessiano de la -log prior:
      hessi.prior[[k]] = Xt.Weig.X / ne[k] # por definicion de Hessiano.
      delta11 <- cbind(rel[, k])
      tt.cox <- coxph(Surv(y, delta11) ~ X_matrix)
      parms.hat1[[k]] = tt.cox$coefficients # coeficiente del modelo
    }
  }  
    parms.hat.12 <- c(parms.hat1[[1]], parms.hat1[[2]]) # la longitud de cada una ya vendrá dado por el vector p11.
    SigmaM.hat.12 <- bloque_matriz(SigmaM.hat) # la matriz diagonal por bloques que corresponde con la matriz de 
    # varianzas-covarianzas resultante para el modelo de riesgos competitivos.
    
    tt1 = optim(par = parms.hat.12, fn = h.parci.SigmaM.normal.cr, 
              gr = grad.anal.cr, Z1 = as.matrix(X1), Z2 = as.matrix(X2), delta = rel, vari = SigmaM.hat.12,
              variinv1 = hessi.prior[[1]], variinv2 = hessi.prior[[2]], p1 = p11[1], p2 = p11[2], n = n, 
              hessian = FALSE, method = "BFGS")
  
    # primer argumento par, que inicializa los valores de los parámetros que
    # optimiza, en este caso, los valores de los coeficientes beta
    
   # Hessiano de la verosimilitud parcial :
    hessi.lk = -hessian(lp.like.u.cr, x = tt1$par, Z1 = as.matrix(X1), Z2 = as.matrix(X2), delta = rel, p1 = p11[1], 
                       p2 = p11[2], n = n)
  
  hessi = hessi.lk + solve(SigmaM.hat.12)
  lm1.lapl = lp.like.u.cr(beta1 = tt1$par, Z1 = as.matrix(X1), Z2 = as.matrix(X2), delta = rel, n = n[1], p1 = p11[1],   p2 = p11[2]) + dmnorm(tt1$par, varcov = SigmaM.hat.12, log=TRUE) + (p11[1] + p11[2]) / 2 * log(2 * pi) - 0.5 * log(det(hessi))
  
  return(list(lm1 = lm1.lapl, betas = tt1$par))
}

h.parci.SigmaM.normal.cr = function(betas, Z1, Z2, delta, vari, variinv1 = NULL, variinv2 = NULL, p1, p2, n) { 
  return(-lp.like.u.cr(beta1 = betas, Z1, Z2, delta, n, p1, p2) - dmnorm(betas, varcov = vari, log = TRUE))
}


# la verosimilitud parcial de un modelo de riesgos competitivos.

lp.like.u.cr = function(beta1, Z1, Z2, delta, n, p1, p2) {
  aux <- NULL # inicializo un vector aux nulo
  termi <- NULL 
  termi_total <- NULL 
  Z.full1 <- list() # inicializo una lista, donde en cada entrada de la lista me guardo
  # los individuos que van a participar porque han tenido un fallo por la causa k.
  for (k in 1:dim(delta)[2]) {
    if (k == 1){
      beta <- beta1[1:p1]
      Z <- Z1
    } else {
      beta <- tail(beta1, p2)
      Z <- Z2
    }
    matriz_evento <- matrix(data = NA, nrow = length(which(delta[, k] == 1)), ncol = dim(Z)[2])
    Z.full1[[k]] <- matriz_evento
    c1 <- 0 # inicializo un número a 0, puedo ir añadiendo más valores a las
    # entradas de la lista.
    for (j in 1:dim(delta)[1]) {
      if(delta[j, k] == 1) {
        c1 <- c1 + 1
        Z.full1[[k]][c1, ] <- Z[j, ] 
      }
    }
    eta_k <- Z.full1[[k]] %*% beta
    e.eta_k <- exp(eta_k)
    eta_total <- Z %*% beta
    e.eta_total = exp(eta_total)
    sumaaa <- NULL
    k1 <- 0
    k2 <- 0
    if (dim(e.eta_k)[2] != 0) {
    for (i in 1:length(e.eta_k)) {
      enter <- TRUE
      for (j1 in (k2 + 1):length(e.eta_total)) { #todo está ordenado gracias a X.full, que viene sacada del
        # dataframe ordenado.
        if (e.eta_k[i] == e.eta_total[j1] & delta[j1, k] == 1 & enter == TRUE) { 
          k1 <- k1 + 1 
          sumaaa[k1] <- sum(e.eta_total[j1:n]) 
          enter <- FALSE
          k2 <- j1
        }
      }
    }
    termi <- eta_k - log(sumaaa) 
    termi_total[k] <- sum(termi)
  } else {
    # bajo el modelo nulo para ese riesgo
    n = length(delta_dummy[,1])
    n.i = n:1
    Bs = 1/n.i
    termi_total[k] = sum(delta_dummy[,k] * log(Bs))
  }
}
  aux <- sum(termi_total) # lo sumo para todos los k, todos los riesgos
  return(aux)
}


# marginal para \beta=0
l.margi0.parcial.cr = function(delta_dummy) {
  term <- NULL
  for (k in 1:dim(delta_dummy)[2]) {
    n = length(delta_dummy[, k])
    n.i = n:1
    Bs = 1 / n.i
    term[k] = sum(delta_dummy[, k] * log(Bs))
  }
  aux <- sum(term)
  return(aux)
}

# crea una matriz diagonal por bloques.
bloque_matriz = function(mlist) {
  a1 <- NULL
  a2 <- NULL
  b1 <- NULL
  b2 <- NULL
  dim_matrix <- lapply(mlist, dim)   # dimensiones de la matriz
  for (k in 1:2) {
  if(is.null(dim_matrix[[k]])) {
    dim_matrix[[k]] <- c(0, 0)
  }
  }
  mb <- matrix(0, nrow = sum(as.numeric(lapply(dim_matrix, function(x) x[[1]]))), 
               ncol=sum(as.numeric(lapply(dim_matrix, function(x) x[[2]]))))
  for (k in 1:2) {
  if (dim_matrix[[k]][1] != 0) {
  for (i in k) {
    a1[i] <- ifelse(i==1, 1, (cumsum(sapply(dim_matrix, function(x) sum(x[[1]]))) + 1)[i - 1])
    a2[i] <- cumsum(sapply(dim_matrix, function(x) sum(x[[1]])))[i]
    b1[i] <- ifelse(i==1, 1, (cumsum(sapply(dim_matrix, function(x) sum(x[[2]]))) + 1)[i - 1])
    b2[i] <- cumsum(sapply(dim_matrix, function(x) sum(x[[2]])))[i]
    mb[a1[i]:a2[i], b1[i]:b2[i]] <- mlist[[i]]
   } 
  }
 }
  return(mb)  # resultado
}

grad.anal.c.cr = function(beta1, Z1, Z2, delta, variinv1, variinv2, p1, p2, n, vari = NULL) {
  grad2 <- NULL
  grad1 <- matrix(NA, ncol = 2, nrow = p)
  for (k in 1:2) {
    if (k == 1) {
      beta <- beta1[1:p1]
      Z <- Z1
      p <- p1
      variinv <- variinv1
    } else {
      beta <- tail(beta1, p2)
      Z <- Z2
      p <- p2
      variinv <- variinv2
    }
    if (dim(Z)[2] != 0) {
    tX = t(Z)
    eta = exp(Z %*% beta)
    Xk = array(NA, dim = c(p, n))
    for(i in 1:(n - 1)) { 
      psi.ki = sum(eta[i:n])
      Xk[, i] = (tX[, i:n] %*% eta[i:n] / psi.ki)
    }
    Xk[, n]=tX[, n]
    
    grad = (Xk - t(Z)) %*% as.numeric(delta_dummy[, k] == 1)
    grad1 = grad + variinv %*% beta # El último término es el gradiente de la prior.
    grad2 <- c(grad2, grad1)
    } else {
    grad1 <- NULL
    grad2 <- c(grad2, grad1)
    }
  }
  return(c(grad2))
}

grad.anal.cr = cmpfun(grad.anal.c.cr)
