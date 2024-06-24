# Selección de variables en modelo de riesgos competitivos
INTRODUCCIÓN

Bienvenido a este repositorio de GitHub, que contiene los ficheros R que corresponden con el trabajo final de máster titualdo "Selección de variables en modelo de riesgos competitivos". Este trabajo fue desarrolado junto a Anabel Forte, profesora de la Universidad de Valencia y en colaboración con INCLIVA, para el Máster de Bioestadística de la Universidad de Valencia, España. 

Si tiene alguna duda, no dude en contactar conmigo mediante el correo ancruzhe@alumni.uv.es. Para obtener el banco de datos original, contactad con el autor Juan Antonio Carbonell Asins mediante el correo electrónico jacarbonell@incliva.es.

Este GitHub contiene diferentes ficheros:

1.- La carpeta Selección de variables contiene dos archivos: script.cr.R que contiene la implementación del algoritmo de Gibbs Sampling para el muestreo de diferentes modelos y otro archivo llamado library.cr.R, que contiene 7 funciones para poder realizar los cálculos necesarios para que el algoritmo Gibbs Sampling desarrollado en el script.cr.R funcione correctamente. Bajo las condiciones establecidas en el documento pdf y en el script.cr.R, el código tardó 5 días en mostrar resultados. Si prefiere no ejecutar el script.cr.R, la carpeta Resultados contiene el archivo RData con el resultados llamado: "datosfinalesN500-burnin50.RData" correspondiente al resultado del script.cr.R.

2.- La carpeta Modelo Multivariante contiene el archivo de R llamado MJB_multivariante.R que contiene el modelo multivariante. Este modelo no tardó mucho, alrededor de 1-2 horas. Aun así, la carpeta Resultados contiene el archivo RData correspondiente al resultado del modelo multivariante llamado: "Modelomultivariante.RData".

3.- La carpeta Resultados contiene varios ficheros RData, algunos ya nombrados anteriormente. Los archivos "datosfinalesN500-burnin50.RData" y "Modelomultivariante.RData" contienen los resultados de la selección de variables y del modelo multivariante, respectivamente. No obstante, también contiene otros elementos claves para el desarrollo de los archivos .R nombrados en los puntos 1 y 2. Concretamente, tenemos la matriz de diseño guardada en "matriz_diseño.RData", el banco de datos final depurado "datos_finales.RData" y el archivo "matrizdelmodelo.RData" que corresponde con el banco de datos generado que es una de las entradas para el modelo generar el algoritmo Gibbs Sampling.

4.-
