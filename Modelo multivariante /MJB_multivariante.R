load("Resultados/datosfinalesN500-burnin50.RData")
load("Resultados/matrizdelmodelo.RData")
# la matriz de diseño para cada una de los riesgos será.

X.full = dat[, !colnames(dat) %in% c("delta_1", "delta_2", "time", "order_time")]
X <- list()
X1 <- X.full[, as.vector(which(res.bvs.gibbs$hpm[1:13] == 1))]
X2 <- X.full[, as.vector(which(res.bvs.gibbs$hpm[14:26] == 1))]
X2 <- X2[, -c(1, 2)]

set.seed(123) # semilla para que los datos sean reproducibles.
delta_dummy <- matrix(c(dat$delta_1, dat$delta_2), ncol = 2)
cat(file = "modelo.multivariable", "
model {
  for(i in 1:n) {
    for(k in 1:Nrisks) {
      # función de riesgo basal Weibull
      base[i, k] <- lambda[k] * alpha[k] * pow(t[i], alpha[k] - 1)
      
      # Covariables y coeficientes diferentes para cada riesgo k
       elinpred[i, k] <- ifelse (k > 1, exp(X2[i, ] %*% beta2), exp(X1[i, ] %*% beta1))
      
      # Funcion logHaz
      logHaz[i, k] <- log(base[i, k] * elinpred[i, k])
      
      # Funcion logarítmica de supervivencia
      logSurv[i, k] <- -lambda[k] * pow(t[i], alpha[k]) * elinpred[i, k]
    }
    
    # Definicion de la verosimilitud utilizando el truco de ceros
    phi[i] <- 100000 - inprod(delta[i, ], logHaz[i, ]) - sum(logSurv[i, ])
    zeros[i] ~ dpois(phi[i])
    
  }
  
  # Distribuciones a priori
  for(l in 1:Nbetas1) {
    beta1[l] ~ dnorm(0, 0.001)
  }
  for(l in 1:Nbetas2) {
    beta2[l] ~ dnorm(0, 0.001)
  }
  for(k in 1:Nrisks) {
    alpha[k] ~ dunif(0, 10)
    lambda[k] ~ dgamma(0.01, 0.01)
  }
  
  # Medidas de asociación
  for (i in 1:Nbetas1){
    Pbeta1[i] <- step(beta1[i])
  }
  for (i in 1:Nbetas2){
    Pbeta2[i] <- step(beta2[i])
  }
for (k in 1:Nrisks){
  Palpha[k] <- step(alpha[k])
  Plambda[k] <- step(lambda[k])
  }
}
")

# datos
d.jags <- list(n = dim(dat)[1], t = as.vector(dat$time), X1 = X1, X2 = X2,
              delta = delta_dummy, zeros = rep(0, dim(dat)[1]), Nbetas1 = ncol(X1), Nbetas2 = ncol(X2),
              Nrisks = ncol(delta_dummy))

# valores iniciales
i.jags <- function() {
  list(beta1 = rnorm(ncol(X1)), 
       beta2 = rnorm(ncol(X2)),
       lambda = runif(ncol(delta_dummy)),
       alpha = runif(ncol(delta_dummy))
      )
}
p.jags <- c("beta1", "beta2", "Pbeta1", "Pbeta2", "alpha", "lambda", "Palpha", "Plambda")
library(jagsUI)
n.thin <- 10 #tasade'thining'
Resul.j1 <- jags(data = d.jags, inits = i.jags, parameters.to.save = p.jags,
                 model.file = "modelo.multivariable", n.chains = 3, n.thin = n.thin,
                 n.iter = 102000, n.burnin = 2000)
save(Resul.j1, file = "Resultados/Modelomultivariante.RData")

