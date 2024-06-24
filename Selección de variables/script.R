rm(list = ls()) # borramos todo lo que haya en el environment

# se inicializan algunos paquetes necesarios para su funcionamiento.
library(compiler)
library(mvtnorm)
library(plyr)
library(BVSNLP)
library(MASS)
source('Selección de variables/library.cr.R') # cargamos la librería donde tenemos todas las funciones necesarias
# para realizar la selección de variables en el contexto bayesiano con riesgos competitivos. 

load(file = "Resultados/datos_finales.RData") # los datos, ya listos para su uso, se descargan de un RData.
load(file = "Resultados/matriz_diseño.RData") # matriz de diselo que al final volvemos a realizarla, y esta nos sirve 
# como guía para ver que todo está bien.
 
Xin <- datos_tfm11_def
# Cálculo de la matriz de diseño Xn, necesaria para poner como variables dummy aquellas 
# variables que son  categóricas.
Xn = matrix(c(as.integer(Xin$EtiologiaCat == 1), 
            as.integer(Xin$EtiologiaCat == 2), 
            as.integer(Xin$EtiologiaCat == 3), 
            as.integer(Xin$EtiologiaCat == 4), 
            Xin$Edad, Xin$AlbgDL, 
            as.integer(Xin$Bblockers == 1), 
            Xin$Br,Xin$Cr, 
            as.integer(Xin$Diabetes == 1), Xin$INN,  
            as.integer(Xin$Sexo == 1), Xin$ALS), ncol = 13)

vec_names <- c("Etiologia_Viral", "EtiologiaCat_Autoinmune", "EtiologiaCat_NAFLD", "EtiologiaCat_Otras",
             "Edad", "AlbgDL", "Bblocker_Si", "Br", "Cr", "Diabetes_Si", "INN", "Sexo_Mujer", "ALS") 
  
colnames(Xn) <- vec_names 
sum(Xn1 == Xn) == dim(datos_tfm11_def)[1] * (dim(datos_tfm11_def)[2] - 2) # la matriz de diseño está bien creada.

# variable dummy para la variable COT.
delta_dummy <- matrix(NA, ncol = 2, nrow = dim(Xn)[1]) 

for (i in 1:length(datos_tfm11_def$COT)){
  delta_dummy[i, 1] <- ifelse((datos_tfm11_def$COT[i] == 2 || datos_tfm11_def$COT[i] == 0), 0, 1) # para la primera 
  # columna de delta_dummy en el que será 0 siempre que sea 0 o 2 y 1 si es 1.
  delta_dummy[i, 2] <- ifelse((datos_tfm11_def$COT[i] == 1 || datos_tfm11_def$COT[i] == 0), 0, 1) # para la segunda
  # columna de delta_dummy, que será 0 siempre que sea 0 o 1 y 1 siempre que sea 2.
}
# Se crea un dataframe ordenado por tiempos, ya que es interesenta para el trabajo.
timee <- Xin$TECOT

# Se crea una variable nueva para poder ordenar el dataframe según tiempos. Esta variable determinará los tiempos 
# diferentes para cada indiviudo, sabiendo que si hay coincidencias de tiempos en dos individuos, no tienen porque 
# ocurrir en el mismo momento específico. 

table(timee)
sum(duplicated(timee))

datos_tfm11_def$timee_new <- matrix(data = NA, nrow = 734, ncol = 1)
for (i in 1:length(table(datos_tfm11_def$TECOT))) {#recorro todos los elementos de la tabla
  if (table(datos_tfm11_def$TECOT)[[i]] > 1) { #si el elemento se repetia más de una vez (>1)
    if(table(datos_tfm11_def$TECOT)[[i]] < 10) { # y si no es <10 
      num_rep <- as.numeric(names(table(datos_tfm11_def$TECOT))[i]) # búsqueda del número que se repite
      position <- which(datos_tfm11_def$TECOT == num_rep) # búsqueda de la posición en mi dataframe del que se repite.
      for (k in 1:length(position)) {
        datos_tfm11_def$timee_new[position[k]] <- num_rep + 0.1 * (k-1) # cálculo de tiempos diferentes
      } 
    } else {
      if (table(datos_tfm11_def$TECOT)[[i]] > 10) { # si el elemento repetido es (>10)
        num_rep<-as.numeric(names(table(datos_tfm11_def$TECOT))[i]) # búsqueda del número repetido
        position<-which(datos_tfm11_def$TECOT == num_rep) # búsqueda de las posiciones en el dataframe inicial
        for (k in 1:length(position)) {
          datos_tfm11_def$timee_new[position[k]] <- num_rep + 0.05 * (k-1) # cálculo de tiempos diferentes
        }
      }
    }
    
  }
  if (table(datos_tfm11_def$TECOT)[[i]] == 1) { # el número no tiene empates
    num_rep_unique <- as.numeric(names(table(datos_tfm11_def$TECOT))[i]) #el numero que se repite (que es único).
    position <- which(datos_tfm11_def$TECOT == num_rep_unique)
    datos_tfm11_def$timee_new[position] <- num_rep_unique # cálculo de tiempos diferentes (en este caso, el mismo número)
  } 
}

sum(duplicated(datos_tfm11_def$timee_new)) # comprobación de tiempos duplicados
p = dim(Xn)[2]
n = dim(Xn)[1]
Xn <- as.data.frame(Xn)
dat = data.frame(Xn, timee, delta_dummy, datos_tfm11_def$timee_new) # creación de un dataframe en el que se 
# incorpora toda la información hasta el momento. 
dat = dat[order(datos_tfm11_def$timee_new), ] # se ordena la base de datos respecto a los tiempos nuevos determinados
names(dat)[p + 1] = "time" 
names(dat)[p + 2] = "delta_1" 
names(dat)[p + 3] = "delta_2"
names(dat)[p + 4] = "order_time" 
save(dat, file = "Resultados/matrizdelmodelo.RData")
cGibbsBvs = function(dat = dat, mat.full.model = NULL, 
                     model.ini = NULL, d.ini = 4, N = 500,burn.in = 50,
                     use.posterior, verbose = FALSE) {
  W.cox <- list()
  neff.cox <- NULL
  X.full <- NULL
  X <- list()
  p_total <- NULL
  
  use.posterior = get(use.posterior)

  time = as.numeric(dat$time)

  if(is.null(mat.full.model)) {
    X.full1 = dat[, !colnames(dat) %in% c("delta_1", "delta_2", "time", "order_time")]
    } else {
    X.full1 = dat[mat.full.model]
    }
    n = nrow(X.full1)
    p = ncol(X.full1)
    Uno = matrix(1, nrow = n, ncol = 1)
    Id = diag(1, n)
    X.full1 = as.matrix(X.full1)
    
    # Centramos la matriz X.full:
    X.full = (Id - Uno %*% t(Uno) / n) %*% X.full1
    

  for (k in 1:2) {
    if (k == 1) {
      delta = cbind(dat$delta_1)
    } else {
      delta = cbind(dat$delta_2)
    }
  # Estimamos el hazar acumulado, bajo H0
  res.cox0 = coxph(Surv(time, delta) ~ 1)
  hazard.cum <- basehaz(res.cox0)
  hazard.cum = hazard.cum[match(time, hazard.cum[, "time"]), "hazard"]
    
  # Definicion matriz var-covar de la prior:
  # en este caso, los tiempos de censura y los tiempos de supervivencia coinciden por lo que  obtenemos: 
  
  H.0.i <- hazard.cum
  S.u = diag(exp(-H.0.i), nrow = n)
  D.u = diag(exp(-H.0.i) * H.0.i, nrow = n)
  W = (Id - S.u + D.u)
  
  # Factor por el que multiplicamos para obtener una matriz unitaria:
  neff.cox[k] = sum(diag(W)) 
  # Usamos la otra definicion de Weig (W_h_0):
  W.cox[[k]] = W - W %*% (Uno %*% t(Uno) / neff.cox[k]) %*% W
  }
    
  p_total = 2 * p # el número total de covariables
  modelsB.PM <- array(rep(0, p_total + 1), dim = c(1, p_total + 1)) # la última columna contiene el producto de BF_{\gamma0} * P(M_\gamma)
  modelsB.PM[1, p_total + 1] <- 1 # B.PM(rep(0,p)) Null vs Null
  set.seed(123) # semilla para que el proyecto sea reproducible.
  model.ini <- NULL
    if(is.null(model.ini)) {
      if(is.null(d.ini)) d.ini = round(p_total / 2)
         model.ini <- sample(c(rep(1, d.ini), rep(0, p_total - d.ini)))
    }
    current.model = model.ini 

   delta_dummy <- matrix(c(dat$delta_1, dat$delta_2), ncol = 2) # variable dummy para la variable COT.
  
  # Marginal para el nulo
  lm0.parci = l.margi0.parcial.cr(delta_dummy)
  
  # Cálculo de m_{\gamma}(y) :
  aux<-use.posterior(y = time,rel = delta_dummy,
                     X1 = X.full[, current.model[1:p] == 1],
                     X2 = X.full[, tail(current.model, p) == 1],
                     Weig.matrix1 = W.cox[[1]], Weig.matrix2 = W.cox[[2]], ne = neff.cox)
  
  B.PMcurrent <- exp(aux$lm1 - lm0.parci - lchoose(p_total, sum(current.model)))
  proposal.model <- current.model

  modelsB.PM <- rbind(modelsB.PM, c(current.model, B.PMcurrent)) #se añade una de las filas
  # al array
  
  visitedmodels.PM <- modelsB.PM # todos los modelos visitados hasta el momento
  
  # Rao-Blackwellized inclusion probabilities :
  incl.probRB <- array(0, dim = c(N + burn.in, p_total))
  for (i in 1:(N + burn.in)) {
    if(verbose) cat("It:", i, "\n")
    for (j in 1:p_total){ #recorre todas las variables y escoge cuales cambiar.
      proposal.model <- current.model
      proposal.model[j] <- 1 - current.model[j]
      
    if (sum(proposal.model) > 0)  { # si al menos hay una covariable activa.
        aux<-use.posterior(y = time, rel = delta_dummy,
                           X1 = X.full[, proposal.model[1:p] == 1], 
                           X2 = X.full[, tail(proposal.model,p) == 1],
                           Weig.matrix1 = W.cox[[1]], Weig.matrix2 = W.cox[[2]], ne = neff.cox)
          } else {
          aux$lm1 = lm0.parci
      }
      #si todos son 0 el sum saldrá 0, y entonces estaremos en el caso de la parcial de lm0.
      B.PMproposal <- exp(aux$lm1 - lm0.parci - lchoose(p_total, sum(proposal.model)))
    
      ratio <- B.PMproposal / (B.PMproposal + B.PMcurrent)
      
      if (runif(1) < ratio) {
        current.model[j] <- proposal.model[j]
        B.PMcurrent <- B.PMproposal
      }
      if(i>1) {
        incl.probRB[i, j] <- incl.probRB[i-1, j] + proposal.model[j] * ratio + (1 - proposal.model[j]) * (1 - ratio)
      }	
    }
    visitedmodels.PM <- rbind(visitedmodels.PM, c(current.model, B.PMcurrent))
  }
  
  for(j in 1:p_total) incl.probRB[, j] <- incl.probRB[, j] / seq(1, (N + burn.in))
  
  visitedmodels.PM = visitedmodels.PM[- (1:(burn.in+2)), ]
  inc.prob = colMeans(visitedmodels.PM[ ,- (p_total + 1)]) #  media de las columnas
  # menos la ultima (cada columna tiene un 0 o un 1), cuantos mas 1, mas mayor la prob de
  # inclusion.
  
  xx = unique(visitedmodels.PM) #se van a repetir y por eso se guardan solo los modelos únicos.
  xx[, p_total + 1] = xx[, p_total + 1] / sum(xx[, p_total + 1]) #normalizamos los valores para que esten entre (0,1)
  xx = xx[order(xx[, p_total + 1],decreasing = TRUE), ] 
  colnames(xx) = c(rep(colnames(X.full), 2), "pp")
  hpm = xx[1, ]
  
  return(list(complete = visitedmodels.PM, hpm = hpm, 
              inc.prob = inc.prob, prob.post = xx, incl.probRB = incl.probRB))
}

# con el hpm tendré el modelo con la probabilidad a posteriori mayor y las probabilidades de inclusión de las covariables indican qué modelo es mpm, aquellas que tengan probabilidad de inclusión > 0,5.

res.bvs.gibbs<-cGibbsBvs(dat = dat, N = 500, burn.in = 50, 
                         use.posterior = "lmarg1.laplace.SigmaM.normal.cr", 
                         mat.full.model = NULL, d.ini = 4, verbose = TRUE)
save(res.bvs.gibbs, file = "Resultados/datosfinalesN500-burnin50.RData")

