# Selección de variables en modelo de riesgos competitivos

Bienvenido a este repositorio de GitHub, que contiene los ficheros R que corresponden con el trabajo final de máster titualdo "Selección de variables en modelo de riesgos competitivos". Este trabajo fue desarrolado junto a Anabel Forte, profesora de la Universidad de Valencia y en colaboración con INCLIVA, con Juan Antonio Carbonell Asins, para el Máster de Bioestadística de la Universidad de Valencia, España. 

Si tiene alguna duda, no dude en contactar conmigo mediante el correo ancruzhe@alumni.uv.es. Para obtener el banco de datos original, contactad con el autor Juan Antonio Carbonell Asins mediante el correo electrónico jacarbonell@incliva.es.

El contenido de este GitHub es el siguiente:

1.- La carpeta Selección de variables contiene dos archivos: script.cr.R que contiene la implementación del algoritmo de Gibbs Sampling para el muestreo de diferentes modelos y otro archivo llamado library.cr.R, que contiene 7 funciones necesarias que realiza los cálculos pertinentes para que el algoritmo Gibbs Sampling desarrollado en el script.cr.R funcione correctamente. Bajo las condiciones establecidas en el documento pdf y en el script.cr.R, el código tardó 5 días en mostrar resultados. Si prefiere no ejecutar el script.cr.R, la carpeta Resultados contiene el archivo RData con el resultados llamado: "datosfinalesN500-burnin50.RData" correspondiente al resultado del script.cr.R.

2.- La carpeta Modelo Multivariante contiene el archivo de R llamado MJB_multivariante.R que contiene el código desarrollado para el modelo multivariante. Este código no tardó mucho en ejecutarse, alrededor de 1-2 horas. Aun así, la carpeta Resultados contiene el archivo RData correspondiente al resultado del modelo multivariante llamado: "Modelomultivariante.RData".

3.- La carpeta Resultados contiene varios ficheros RData, ya nombrados anteriormente, que son los archivos "datosfinalesN500-burnin50.RData" y "Modelomultivariante.RData", que contienen los resultados de la selección de variables y del modelo multivariante, respectivamente. 

4.- La carpeta Archivos_RData contiene los banco de datos que se necesitan para poder llevar a cabo el trabajo. Concretamente, contiene la matriz de diseño, que está guardada en "matriz_diseño.RData", el banco de datos final depurado denotado por "datos_finales.RData", el archivo "matrizdelmodelo.RData" correspondiente con el banco de datos generado, siendo este una de las entradas para el script donde se desarrolla el algoritmo Gibbs Sampling antes nombrado "script.cr.R" y "datos_finales_NA.RData", que contiene la base de datos con los datos faltantes, que nos permite realizar una de las figuras.

5.- La carpeta Figuras contiene todas las figuras incluidas en el trabajo.

6.- El archivo Descriptivo.Rmd contiene el código que genera las figuras que se muestran a lo largo del trabajo, y quedan guardadas en la carpeta Figuras, antes nombrada.
