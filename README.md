# Script de R para la Prueba Q de Normalidad Multivariante
Se presenta un script de R en español para Cuatro Variantes de la Prueba Q de Normalidad Multivariante (Moral, 2023, 2025)

En la primera parte (Apéndice 1) del script, se encuentran los códigos para implementar las dos variantes (aproximación chi-cuadrado y Bootstrap) de la prueba Q de normalidad multivariante de Moral (2023) calculadas con las estadísticas de normalidad univariante de Shapiro y Wilk (1965) estandarizadas por el método de Royston (1992) para muestras de tamaño de 12 a 2000 con k mediciones (QSW).
En la segunda parte (Apéndice 2) del script, se hallan los códigos para las dos variantes (aproximación chi-cuadrado y Bootstrap) de la prueba Q de Moral (2023) calculadas con las estadísticas de normalidad univariante de Shapiro y Francia (1972) estandarizadas por el método de Royston (1993) para muestras de tamaño de 5 a 5000 con k mediciones (QSF).
Las dos variantes de la prueba QSW o QSF de Moral (2023) se complementan con la prueba de normalidad multivariante H de Royston (1983), que se encuentra en el Apéndice 3 del script, para muestras de tamaño de 10 a 2000 con k mediciones.
Se usa un nivel de significancia (α) de 0.05 para la prueba H y la aproximación chi-cuadrado de la prueba Q y 0.1 para el enfoque Bootstrap de la prueba Q. El remuestreo se implementa con 1000 muestras y una semilla (123) para que sea replicable. Estos valores destacados en verde pueden cambiarse.
Finalmente, se incluye el cálculo de la potencia relativa entre la prueba Q (evaluada mediante la aproximación chi-cuadrado o el enfoque Bootstrap) y la prueba H (Royston, 1983). Para implementar correctamente este último análisis del Apéndice 3, es necesario ejecutar previamente la primera parte del script (Apéndice 1: prueba QSW con las estadísticas W de Shapiro y Wilk, 1965, estandarizadas por Royston, 1992) o la segunda parte del script (Apéndice 2: prueba QSF con las estadísticas W’ de Shapiro y Francia, 1972, estandarizadas por Royston, 1993). Dichos scripts generan los valores de: power (potencia bajo la aproximación chi-cuadrado) y boot_power (potencia obtenida mediante Bootstrap). Después de ello, debe correrse la tercera parte del script (Apéndice 3), que calcula power_H para la prueba H. Alternativamente, es posible introducir manualmente los valores de potencia estadística obtenidos en una ejecución previa e independiente de los Apéndices 1 o 2.
Las partes 1 o 2 y 3 del script pueden ejecutarse localmente en una computadora personal o portátil, previa instalación del programa R (disponible en https://cran.r-project.org/bin/windows/base/), preferentemente en combinación con el entorno de desarrollo integrado RStudio (https://posit.co/download/rstudio-desktop/), que facilita la modificación (valores en azul o verde) y ejecución de los códigos. Para su funcionamiento, es necesario tener instalados y activados los paquetes nortest (Gross &, Ligges, 2022), MASS (Ripley, et al., 2025), DescTools (Signorell, 2025), stats (R Core Team, 2025a) y MVN (Korkmaz, 2025) del programa R (R Core Team, 2025b).
Como alternativa, las partes 1 o 2 y 3 del script pueden ejecutarse en línea a través de la plataforma rdrr.io (https://rdrr.io/snippets/), la cual cuenta con los paquetes requeridos preinstalados. Si se quiere guardar la gráfica como un archivo JPG, se debe quitar los símbolos del numeral (hashtag) de las instrucciones: # jpg(), # par() y # dev.off(). Véase Apéndices 1‑3.
Se recomienda utilizar la última versión de R. El paquete 'nortest' no especifica una versión mínima clara; la versión actual 1.0-4 funciona en R 4.4 y versiones anteriores. La versión 7.3-65 del paquete 'MASS' requiere R (≥ 4.4.0). La versión 0.99.60 del paquete 'DescTools' requiere R (≥ 4.2.0). El paquete 'stats' es un paquete base de R y se actualiza junto con R. El paquete 'MVN' necesita versiones de R de 3.5 en adelante.
El script se encuentra aplicado a una muestra simulada con 36 observaciones en seis variables correlacionadas (Cuadro 1). No obstante, puede adaptarse a otros conjuntos de datos modificando la lista de variables, que se destaca en color azul.

Cuadro 1
Condiciones de simulación
Características	Especificación
Tamaño de la muestra original (n)	36
Número de variables (p)	6
Matrices de correlaciones	Correlaciones producto-momento de Pearson empíricas. Una matriz generada de tamaño 36×6 compuesta por vectores gaussianos independientes multiplicada por la matriz diagonal inferior de la factorización de Cholesky de la matriz de correlaciones de las variables originales proporciona la muestra normativa de 6 vectores gaussianos (con 36 dimensiones) correlacionados
Número de réplicas (B)	1000
Tipo de distribución generadora	Muestreo con reposición por fila para preservar las correlaciones de la muestra original. Por una parte, se realiza desde la muestra original para el cálculo de la potencia estadística, generando la distribución Bootstrap empírica. Por otra parte, se ejecuta desde la muestra normativa para calcular la probabilidad Bootstrap y obtener el valor crítico Bootstrap (cuantil de orden .95) para el cálculo de la potencia estadística, generando la distribución Bootstrap normativa.

Para la aplicación de las pruebas Q y H, se recomienda:
1.	Ejecutar el script con las cuatro variantes de la prueba Q y con la prueba H, para obtener una evaluación más completa y robusta del ajuste a la normalidad multivariante.
2.	Elegir el enfoque Bootstrap de la prueba Q cuando no se cumple el supuesto de independencia serial, ya que ofrece una distribución muestral más adecuada para este contexto.
3.	Utilizar las estadísticas W' de Shapiro y Francia al implementar la variante Bootstrap, dado su mejor rendimiento en muestras pequeñas.
4.	Ajustar el nivel de significancia al 10% en el enfoque Bootstrap, considerando su naturaleza conservadora frente a la hipótesis nula de normalidad multivariante.
5.	Evitar restar el número de estadísticas truncadas a cero al calcular los grados de libertad en la variante basada en la aproximación chi-cuadrado de la prueba Q.
6.	Emplear la aproximación chi-cuadrado de la prueba Q únicamente cuando se confirme la independencia serial, mediante pruebas de contraste y análisis gráfico.

Referencias
Gross, J., & Ligges, U. (2022). Package ‘nortest’, Tests for Normality. Version 1.0-4. https://cran.r-project.org/web/packages/nortest/nortest.pdf
Korkmaz, S. (2025). Package ‘MVN’. Multivariate normality tests. https://cran.r-project.org/web/packages/MVN//MVN.pdf
Moral, J. (2023). Proposal and pilot study: a generalization of the W or W' statistic for multivariate normality. Open Journal of Statistics, 13(1), 119‑169. https://doi.org/10.4236/ojs.2023.131008
R Core Team. (2025a). The R Stats package (version 4.6.0). https://stat.ethz.ch/R-manual/R-devel/library/stats/html/00Index.html
R Core Team. (2025b). R: A language and environment for statistical computing. Vienna, Austria: R Foundation for Statistical Computing. https://www.R-project.org/
Ripley, B., Venables, B. Bates, D. M., Hornik, K., Gebhardt, A., & Firth, D. (2025). Support functions and datasets for Venables and Ripley's MASS. Versión 7.3-65. https://doi.org/10.32614/CRAN.paquete.MASS
Royston, J. P. (1983). Some techniques for assessing multivariate normality based on the Shapiro-Wilk W. Journal of the Royal Statistical Society. Series C (Applied Statistics), 32(2), 121‑133. https://doi.org/10.2307/2347291
Royston, J. P. (1992). Approximating the Shapiro-Wilk W-test for non-normality. Statistics and Computing, 2(3), 117‑119. https://doi.org/10.1007/BF01891203
Royston, J. P. (1993). A toolkit of testing for non-normality in complete and censored samples. Journal of the Royal Statistical Society, Series D (The Statistician), 42(1), 37‑43. https://doi.org/10.2307/2348109
Shapiro, S. S., & Francia, R. S. (1972). An approximate analysis of variance test for normality. Journal of the American Statistical Association, 67(337), 215‑216. https://doi.org/10.1080/01621459.1972.10481232
Shapiro, S. S., & Wilk, M. B. (1965) An analysis of variance test for normality (complete samples). Biometrika, 52(3/4), 591‑611. https://doi.org/10.2307/2333709
Signorell, A. (2025). DescTools: Tools for descriptive statistics. Version 0.99.6. https://doi.org/10.32614/CRAN.package.DescTools

