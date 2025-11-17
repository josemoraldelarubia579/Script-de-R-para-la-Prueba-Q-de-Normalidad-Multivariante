# Script de R para la Prueba Q de Normalidad Multivariante
Se presenta un script de R en español para Cuatro Variantes de la Prueba Q de Normalidad Multivariante (Moral, 2023, 2025)
#------------------------------------------------------------------------------------------------------ 
# Apéndice 1
# Prueba Q de normalidad multivariante desde las estadísticas W de normalidad univariante de Shapiro y Wilk (1965) estandarizadas por el método de Royston (1992) para muestras de tamaño de 12 a 2000 con k mediciones.
#-------------------------------------------------------------------------------------------------------

# Cargar cuatro paquetes requeridos.
library(nortest)
library(MASS)
library(DescTools)
library(stats)

# Se recomienda la última versión de R.
# El paquete 'nortest' no especifica una versión mínima clara; la versión actual 1.0-4 funciona en R 4.4 y versiones anteriores.
# La versión 7.3-65 del paquete 'MASS' requiere R (≥ 4.4.0).
# La versión 0.99.60 del paquete 'DescTools' requiere R (≥ 4.2.0).
# El paquete 'stats' es un paquete base de R y se actualiza junto con R.
# Si no se tienen los paquetes, seleccionar espejo más próximo: chooseCRANmirror()
# install.packages("nortest")
# install.packages("MASS")
# install.packages("DescTools")

# Lista de variables originales.
x1 <- c(7, 6, 4, 5, 3, 9, 5, 3, 3, 5, 5, 6, 5, 9, 6, 10, 8, 7, 7, 3, 6, 5, 8, 7, 1, 4, 7, 2, 6, 4, 4, 8, 6, 2, 1, 4)
x2 <- c(6, 7, 4, 4, 5, 7, 6, 5, 3, 4, 7, 5, 1, 8, 5, 10, 8, 6, 8, 2, 4, 2, 7, 9, 3, 3, 5, 3, 5, 4, 6, 7, 6, 6, 1, 4)
x3 <- c(7, 9, 5, 5, 6, 9, 5, 5, 6, 2, 6, 7, 2, 10, 6, 9, 9, 3, 8, 2, 4, 4, 9, 6, 2, 3, 7, 4, 3, 6, 4, 7, 6, 5, 1, 4)
x4 <- c(8, 7, 6, 7, 3, 9, 4, 2, 3, 4, 5, 5, 1, 10, 6, 8, 9, 2, 6, 2, 6, 4, 7, 7, 5, 5, 7, 5, 5, 4, 6, 6, 8, 3, 1, 4)
x5 <- c(6, 7, 4, 6, 3, 6, 8, 2, 4, 4, 8, 8, 4, 10, 6, 8, 7, 4, 6, 2, 2, 4, 9, 9, 3, 4, 4, 6, 4, 6, 6, 7, 7, 3, 1, 3)
x6 <- c(8, 6, 4, 7, 5, 9, 8, 3, 5, 7, 6, 6, 5, 10, 4, 10, 9, 4, 6, 3, 6, 2, 8, 8, 5, 3, 7, 7, 4, 5, 6, 7, 5, 3, 2, 4)
datos_orig <- data.frame(x1, x2, x3, x4, x5, x6)

# Tamaño muestral, número de variables y nivel de significancia.
n <- length(x1)
cat("\nTamaño muestral: n =", n, ".\n")
k <- length(datos_orig)
cat("Número de variables en la muestra original: k =", k, ".\n")
R <- cor(datos_orig)
m_R <- mean(R[lower.tri(R)])
cat("Media de las correlaciones entre las variables: m(R) =", round(m_R, 4), ".\n")
alpha <- 0.05

# Generar las sumas de las combinaciones sin repetición entre las variables.
generar_combinaciones_lineales <- function(df) {
lista_combinaciones <- list()
k <- ncol(df)
# Agregar las variables originales (combinaciones de tamaño 1).
for (i in 1:k) {
lista_combinaciones[[paste0("c", i)]] <- df[, i]}
# Agregar combinaciones de tamaño 2 hasta k.
contador <- k
for (i in 2:k) {
combinaciones <- combn(k, i)
for (j in 1:ncol(combinaciones)) {
idx <- combinaciones[, j]
contador <- contador + 1
lista_combinaciones[[paste0("c", contador)]] <- rowSums(df[, idx])}}
return(lista_combinaciones)}
lista_x <- generar_combinaciones_lineales(datos_orig)

# Cálculo de los parámetros comunes para la estandarización para tamaños muestrales de 12 a 2000.
mu <- 0.0038915 * log(n)^3 - 0.083751 * log(n)^2 - 0.31082 * log(n) - 1.5861
sigma <- exp(0.0030302 * log(n)^2 - 0.082676 * log(n) - 0.4803)
cat("\nParámetros comunes para la estandarización de las", 2^k - 1, "variables creadas por combinación lineal de las", k, "variables originales:\n")
cat("\nValor esperado de ln(1-W): μ =", mu, ".\n")
cat("Desviación estándar de ln(1-W): σ =", sigma, ".\n")

# Aplicar prueba y guardar en tabla.
resultados <- lapply(names(lista_x), function(nombre) {
datos <- lista_x[[nombre]]
resultado <- shapiro.test(datos)
z <- (log(1 - resultado$statistic) - mu) / sigma
p <- 1 - pnorm(z)
Normalidad <- if (resultado$p.value < alpha) "No" else "Si"
return(c(nombre, round(resultado$statistic, 4), round(z, 4), round(p, 4), Normalidad))
})

# Convertir a marco de datos (data frame).
tabla_resultados <- as.data.frame(do.call(rbind, resultados), stringsAsFactors = FALSE)
colnames(tabla_resultados) <- c("Variable", "W", "z", "p", "Normalidad")

# Mostrar tabla.
cat("\nTabla. Prueba de normalidad univariante de Shapiro y Wilk con la estandarización de Royston para las", 2^k-1, "combinaciones lineales de variables\n")
print(tabla_resultados, row.names = FALSE)
cat("\nNota. Variable = suma de la correspondiente combinación sin repetición de las variables; W = estadística W de Shapiro y Wilk (1965); z = valor estandarizado de W mediante la estandarización de Royston (1992); p = valor de probabilidad a la cola derecha de una distribución normal estándar; y Normalidad univariante: Sí = p ≥ α =", alpha, ", se mantiene la hipótesis nula de normalidad; y No = p <", alpha, ", se rechaza.\n")

# Función para calcular q.
calcular_q <- function(lista_datos) {
z_vals <- sapply(lista_datos, function(x) {
resultado_sw <- shapiro.test(x)
z <- (log(1 - resultado_sw$statistic) - mu) / sigma
return(as.numeric(z))
})
z_values <- ifelse(z_vals < 0, 0, z_vals)
a <- sum(z_vals < 0)
q <- sum(z_values^2)
return(list(q = q, a = a, z_vals = z_vals, z_values = z_values))
}

# Valor de probabilidad con la distribución chi-cuadrado estándar.
resultado_q <- calcular_q(lista_x)
q_stat <- resultado_q$q
a <- resultado_q$a
nc <- length(lista_x)
# gl_chisq <- nc - a # con corrección por valores truncados.
gl_chisq <- nc # sin corrección por valores truncados.
q_crit <- qchisq(1 - alpha, df = gl_chisq)
p_value <- pchisq(q_stat, df = gl_chisq, lower.tail = FALSE)
power <- 1 - pchisq(q_crit, df = gl_chisq, ncp = q_stat, lower.tail = TRUE)

# Resultados para QSWa.
cat("\nNúmero de combinaciones entre las", k, "variables originales: nc =", nc, ".\n")
cat("Número de valores z anulados por ser negativos: a =", a, ".\n")
cat("Grados de libertad: gl =", gl_chisq, ".\n")
cat("\nEstadística de contraste: q =", round(q_stat, 4), ".\n")
cat("\nValor de probabilidad y valor crítico desde una distribución chi-cuadrado estándar:\n")
cat("\nValor crítico o cuantil de orden", 1 - alpha," de una distribución chi-cuadrado estándar con", gl_chisq, "grados de libertad: 1-αχ²[gl] =", round(q_crit, 4), ".\n")
if (q_stat <= q_crit) {cat("La hipótesis nula de normalidad multivariante se mantiene a un nivel de significancia de", alpha, "con base en el valor crítico desde la aproximación chi-cuadrado.\n")
} else {cat("La hipótesis nula de normalidad multivariante se rechaza a un nivel de significancia de", alpha, "con base en el valor crítico desde la aproximación chi-cuadrado.\n")}
cat("Valor de probabilidad a la cola derecha en una distribución chi-cuadrado estándar con", gl_chisq, "grados de libertad: p =", round(p_value, 4), ".\n")
if (p_value < alpha) {cat("La hipótesis nula de normalidad multivariante se rechaza a un nivel de significancia de", alpha, "con base en el valor de probabilidad desde la aproximación chi-cuadrado.\n")
} else {cat("La hipótesis nula de normalidad multivariante se mantiene a un nivel de significancia de", alpha, "con base en el valor de probabilidad desde la aproximación chi-cuadrado.\n")}
cat("Potencia estadística desde la aproximación chi-cuadrado para un nivel de significancia de", alpha, ": phi =", power, ".\n")
cat("Si se rechaza la hipótesis nula, la potencia estadística debería ser superior a 0.5, preferiblemente mayor que 0.8; mientras que, si no se rechaza la hipótesis nula, debería ser inferior a 0.5, preferiblemente menor que 0.2. En otro caso, el resultado es contradictorio o cuestionable.\n")

# Contraste de la independencia serial en los valores z.
z_w <- resultado_q$z_values
cat("\nContraste de la independencia entre los", nc, "valores de W logarítmicamente transformados, estandarizados y truncados (a 0 cuando son negativos), usados para calcular q(z_w).\n")
RunsTest(z_w, y = NULL, alternative = "two.sided", exact = TRUE, correct = FALSE)
lag_max <- min(10, round(nc/5, 0)) # Regla de Hyndman y Athanasopoulos (2021) para determinar el rezago máximo en series temporales no estacionarias.
# lag_max <- 12 * (n/100)^0.25 # Regla de Schwert para determinar el rezago máximo.
ljung_box_test <- lapply(1: lag_max, function(k) Box.test(z_w, lag = k, type = "Ljung-Box", fitdf = 0))
ljung_box_tabla <- data.frame(
Rezago = 1:lag_max,
Chi_cuadrado = sapply(ljung_box_test, function(x) round(x$statistic, 4)),
gl = sapply(ljung_box_test, function(x) round(x$parameter, 4)),
p = sapply(ljung_box_test, function(x) round(x$p.value, 4))
)
cat("\nTabla. Prueba de Ljung y Box de independencia serial en la secuencia de valores z_w\n")
print(ljung_box_tabla)
cat("\nNota. Rezago = intervalo entre los valores que conforman los pares de datos en la autocorrelación, X² = estadística de contraste, gl = grados de libertad y p = probabilidad a la cola derecha en una distribución chi-cuadrado estándar.\n")

# Función Bootstrap preservando las correlaciones de la muestra original.
resample_preserving_corr <- function(df, B) {
replicate(B, {
df[sample(1:nrow(df), size = nrow(df), replace = TRUE), ]
}, simplify = FALSE)
}

# Función para generar lista de combinaciones.
generar_lista_boot <- function(df) {
lista_boot <- setNames(as.list(df), paste0("c", 1:ncol(df)))
for (i in 2:ncol(df)) {
combinaciones <- combn(ncol(df), i)
for (j in 1:ncol(combinaciones)) {
idx <- combinaciones[, j]
nombre <- paste0("c", length(lista_boot) + 1)
lista_boot[[nombre]] <- rowSums(df[, idx])
}
}
return(lista_boot)}

# Simulaciones Bootstrap preservando correlaciones.
set.seed(123) # Se fija una semilla para reproducibilidad
B <- 1000 # Número de remuestreos con reposición desde la muestra original
boot_datasets <- resample_preserving_corr(datos_orig, B)
q_boot_vals <- sapply(boot_datasets, function(df_boot) {
lista_boot <- generar_lista_boot(df_boot)
calcular_q(lista_boot)$q
})

# Valor de probabilidad en la distribución Bootstrap empírica.
boot_p_value <- mean(q_boot_vals >= q_stat)
cat("\nDistribución Bootstrap (con correlación conservada) empírica o desde la muestral original:\n")
cat("Número de muestras Bootstrap: B =", B, ".\n")
cat("Valor de probabilidad Bootstrap a la cola derecha en la distribución Bootstrap empírica: p_boot =", round(boot_p_value, 4), ".\n")
if (boot_p_value < alpha) {cat("La hipótesis nula de normalidad multivariante se rechaza a un nivel de significancia de", alpha, "con base en el valor de probabilidad Bootstrap calculado mediante la distribución Bootstrap empírica.\n")
} else {cat("La hipótesis nula de normalidad multivariante se mantiene a un nivel de significancia de", alpha, "con base en el valor de probabilidad Bootstrap calculado mediante la distribución Bootstrap empírica.\n")}

# Generar una muestra normativa con probabilidades equiespaciadas y la transformación Cholesky de la matriz de correlaciones de la muestra original.
p <- seq(0.5 / n, 1 - 0.5 / n, length.out = n)
datos_indep <- matrix(NA, nrow = n, ncol = k)
set.seed(123)
for (j in 1:k) {datos_indep[, j] <- qnorm(sample(p))}
mat_cor <- cor(datos_orig)
mat_chol <- chol(mat_cor)
datos_normativos <- datos_indep %*% mat_chol

# Definir las variables-suma de todas las posibles combinaciones entre las variables originales.
lista_x_norm <- generar_combinaciones_lineales(datos_normativos)
result_q_norm <- calcular_q(lista_x_norm)

# Generar distribución Bootstrap normativa con remuestreo de filas completas.
B <- 1000 # Número de remuestreos con reposición desde la muestra normativa
q_boot_norm_vals <- replicate(B, {
indices <- sample(1:n, size = n, replace = TRUE)
muestra_boot <- datos_normativos[indices, ]
lista_boot_n <- generar_combinaciones_lineales(muestra_boot)
calcular_q(lista_boot_n)$q
})

# Valor de probabilidad y valor crítico en la distribución Bootstrap normativa.
boot_norm_p_value <- mean(q_boot_norm_vals >= q_stat)
# orden_cuantil <- 1-alpha # Sin compensación por cola derecha muy pesada.
orden_cuantil <- 1-2*alpha # Con compensación por cola derecha muy pesada.
valor_critico_boot_norm <- quantile(q_boot_norm_vals, probs = orden_cuantil)
cat("\nDistribución Bootstrap normativa (con correlación conservada):\n")
cat("Número de muestras Bootstrap: B =", B, ".\n")
cat("Valor de la estadística q en la muestra normativa: q_mvn_boot =", result_q_norm$q, ".\n")
cat("Media de la distribución Bootstrap normativa: m(mvn_boot_dist) =", mean(q_boot_norm_vals), ".\n")
cat("Mediana de la distribución Bootstrap normativa: mdn(mvn_boot_dist) =", median(q_boot_norm_vals), ".\n")
p_value_mdn_norm <- 2 * min(mean(q_boot_vals >= median(q_boot_norm_vals)), mean(q_boot_vals <= median(q_boot_norm_vals)))
cat("Probabilidad de que la distribución Bootstrap empírica esté centrada en la mediana de la distribución Bootstrap normativa: p_2 colas =", round(p_value_mdn_norm, 4), ".\n")
if (p_value_mdn_norm < alpha) {cat("Se rechaza hipótesis nula de medianas equivalentes entre las distribuciones Bootstrap normativa y empírica en un contraste a dos colas con un nivel de significancia de", 1-orden_cuantil, ".\n")
} else {cat("Se mantiene hipótesis nula de medianas equivalentes entre las distribuciones Bootstrap normativa y empírica en un contraste a dos colas con un nivel de significancia de ", 1-orden_cuantil, ".\n")}
cat("Valor crítico Bootstrap o cuantil de orden", orden_cuantil, "de la distribución Bootstrap normativa: q_crit =", round(valor_critico_boot_norm, 4), ".\n")
if (q_stat <= valor_critico_boot_norm) {cat("La hipótesis nula de normalidad multivariante se mantiene a un nivel de significancia de", 1-orden_cuantil, "con base en el valor crítico Bootstrap (normativo).\n")
} else {cat("La hipótesis nula de normalidad multivariante se rechaza a un nivel de significancia de", 1- orden_cuantil, "con base en el valor crítico Bootstrap (normativo).\n")}
cat("Valor de probabilidad a la cola derecha en la distribución Bootstrap normativa: p_mvn_boot =", round(boot_norm_p_value, 4), ".\n")
if (boot_norm_p_value < 1-orden_cuantil) {cat("La hipótesis nula de normalidad multivariante se rechaza a un nivel de significancia de", 1-orden_cuantil, "con base en el valor de probabilidad Bootstrap calculado mediante la distribución Bootstrap normativa.\n")
} else {cat("La hipótesis nula de normalidad multivariante se mantiene a un nivel de significancia de", 1-orden_cuantil, "con base en el valor de probabilidad Bootstrap calculado mediante la distribución Bootstrap normativa.\n")}
# Potencia Bootstrap.
boot_power <- mean(q_boot_vals > valor_critico_boot_norm)
cat("Potencia estadística Bootstrap para un nivel de significancia de", 1-orden_cuantil, ": phi_boot =", boot_power, ".\n")

# Representación de las distribuciones muestral Bootstrap (empírica y normativa).
# Quitar los símbolos hashtag precedentes para grabar como un archivo JPG o TIFF
# jpeg("Hist_curva_dens1.jpg", width = 1200, height = 900, units = "px", res = 300)
# tiff("Hist_curva_dens1.tiff", width = 1200, height = 900, units = "px", res = 300)
# par(mar = c(4.5, 4.5, 0.5, 0.5), cex.axis = 0.8)

# Curva de densidad de la distribución chi-cuadrado estándar como gráfico base.
x_seq <- seq(0, max(qchisq(0.999, df = gl_chisq), max(q_boot_vals) +1, q_stat +1, max(q_boot_norm_vals) +1), length.out = 1000)
y_chi_sq <- dchisq(x_seq, df = gl_chisq)
plot(x_seq, y_chi_sq, type = "l", lwd = 3, col = "darkgreen",
main = "", xlab = "Valores q", ylab = "Densidad",
xlim = c(0, max(qchisq(0.999, df = gl_chisq), max(q_boot_vals) +1, q_stat +1, max(q_boot_norm_vals) +1)),
ylim = c(0, max(density(q_boot_vals)$y, y_chi_sq)))
# Histograma amarillo de la distribución Bootstrap empírica de la estadística de contraste q.
hist(q_boot_vals, breaks = 20, freq = FALSE, col = rgb(1, 1, 0, 0.5), border = "yellow2", add = TRUE)
# Curva amarilla de densidad de la distribución Bootstrap empírica de la estadística de contraste q.
lines(density(q_boot_vals), col = "yellow3", lwd = 2)
# Curva roja de densidad bootstrap normativa.
lines(density(q_boot_norm_vals), col = "red", lwd = 2)
# Línea vertical violeta para el valor de la estadística de contraste q observado.
abline(v = q_stat, col = "purple", lwd = 2, lty = 2)
# Línea vertical roja para el valor crítico normativo bootstrap.
abline(v = valor_critico_boot_norm, col = "red", lwd = 2, lty = 2)
# Línea vertical verde para el valor crítico de la distribución chi-cuadrado estándar.
abline(v = qchisq(1-alpha, df = gl_chisq), col = "darkgreen", lwd = 2, lty = 2)
# dev.off() # Quitar el símbolo hashtag precedente cuando se graba la figura.

cat("\nFigura 1. Histograma (amarillo) con las curvas de densidad de las distribuciones Bootstrap empírica (amarilla), Bootstrap normativa (roja) y chi-cuadrado estándar (verde). El valor de la estadística de contraste q se representa mediante una línea vertical violeta (q =", q_stat, "), el valor crítico de la distribución Bootstrap normativa mediante una línea vertical roja (", orden_cuantil, "_q_mvn_boot =", valor_critico_boot_norm, ") y el valor crítico de la distribución chi-cuadrado estándar con", gl_chisq, "grados de libertad mediante una línea vertical verde (", 1-alpha, "χ²[", gl_chisq, "] =", q_crit, ").\n")
cat("\nNota. Histograma amarillo (regla de la Universidad de Rice) = densidades de la distribución Bootstrap empírica de la estadística de contraste q,
curva amarilla = curva de densidad de la distribución Bootstrap empírica de la estadística de contraste q,
curva roja = curva de densidad de la distribución Bootstrap normativa (generada desde variables normalmente distribuidas con la misma estructura correlacional que las originales) de la estadística de contraste q,
curva verde = distribución chi-cuadrado estándar con", gl_chisq, "grados de libertad,
línea vertical violeta = valor de la estadística de contraste q en la muestra original,
línea vertical roja = valor crítico o cuantil de orden", orden_cuantil, "de la distribución Bootstrap normativa,
línea vertical verde = valor crítico o cuantil de orden", 1-alpha, "de una distribución chi-cuadrado estándar con", gl_chisq, "grados de libertad.\n")

#------------------------------------------------------------------------------------------------------ 
# Apéndice 2
# Prueba Q' de normalidad multivariante desde las estadísticas W' de normalidad univariante de Shapiro y Francia (1972) estandarizadas por el método de Royston (1993) para muestras de tamaño de 5 a 5000 con k mediciones.
#-------------------------------------------------------------------------------------------------------


# Cargar cuatro paquetes requeridos.
library(nortest)
library(MASS)
library(DescTools)
library(stats)

# Lista de variables originales.
x1 <- c(7, 6, 4, 5, 3, 9, 5, 3, 3, 5, 5, 6, 5, 9, 6, 10, 8, 7, 7, 3, 6, 5, 8, 7, 1, 4, 7, 2, 6, 4, 4, 8, 6, 2, 1, 4)
x2 <- c(6, 7, 4, 4, 5, 7, 6, 5, 3, 4, 7, 5, 1, 8, 5, 10, 8, 6, 8, 2, 4, 2, 7, 9, 3, 3, 5, 3, 5, 4, 6, 7, 6, 6, 1, 4)
x3 <- c(7, 9, 5, 5, 6, 9, 5, 5, 6, 2, 6, 7, 2, 10, 6, 9, 9, 3, 8, 2, 4, 4, 9, 6, 2, 3, 7, 4, 3, 6, 4, 7, 6, 5, 1, 4)
x4 <- c(8, 7, 6, 7, 3, 9, 4, 2, 3, 4, 5, 5, 1, 10, 6, 8, 9, 2, 6, 2, 6, 4, 7, 7, 5, 5, 7, 5, 5, 4, 6, 6, 8, 3, 1, 4)
x5 <- c(6, 7, 4, 6, 3, 6, 8, 2, 4, 4, 8, 8, 4, 10, 6, 8, 7, 4, 6, 2, 2, 4, 9, 9, 3, 4, 4, 6, 4, 6, 6, 7, 7, 3, 1, 3)
x6 <- c(8, 6, 4, 7, 5, 9, 8, 3, 5, 7, 6, 6, 5, 10, 4, 10, 9, 4, 6, 3, 6, 2, 8, 8, 5, 3, 7, 7, 4, 5, 6, 7, 5, 3, 2, 4)
datos_orig <- data.frame(x1, x2, x3, x4, x5, x6)

# Tamaño muestral, número de variables y nivel de significancia.
n <- length(x1)
k <- length(datos_orig)
cat("\nTamaño muestral: n =", n, ".\n")
cat("Número de variables en la muestra original: k =", k, ".\n")
R <- cor(datos_orig)
m_R <- mean(R[lower.tri(R)])
cat("Media de las correlaciones entre las variables: m(R) =", round(m_R, 4), ".\n")
alpha <- 0.05

# Generar las sumas de las combinaciones sin repetición entre las variables.
generar_combinaciones_lineales <- function(df) {
lista_combinaciones <- list()
k <- ncol(df)
# Agregar las variables originales (combinaciones de tamaño 1).
for (i in 1:k) {
lista_combinaciones[[paste0("c", i)]] <- df[, i]}
# Agregar combinaciones de tamaño 2 hasta k.
contador <- k
for (i in 2:k) {
combinaciones <- combn(k, i)
for (j in 1:ncol(combinaciones)) {
idx <- combinaciones[, j]
contador <- contador + 1
lista_combinaciones[[paste0("c", contador)]] <- rowSums(df[, idx])}}
return(lista_combinaciones)}
lista_x <- generar_combinaciones_lineales(datos_orig)

# Cálculo de los parámetros comunes para la estandarización.
u <- log(log(n)) - log(n)
mu <- 1.0521 * u - 1.2725
v <- log(log(n)) + 2 / log(n)
sigma <- -0.26758 * v + 1.0308
cat("\nParámetros comunes para la estandarización de las", 2^k - 1, "variables creadas por combinación lineal de las", k, "variables originales:\n")
cat("\nValor esperado de ln(1-W'): μ =", mu, ".\n")
cat("Desviación estándar de ln(1-W'): σ =", sigma, ".\n")

# Aplicar prueba y guardar en tabla.
resultados <- lapply(names(lista_x), function(nombre) {
datos <- lista_x[[nombre]]
resultado <- shapiro.test(datos)
z <- (log(1 - resultado$statistic) - mu) / sigma
p <- 1 - pnorm(z)
Normalidad <- if (resultado$p.value < alpha) "No" else "Si"
return(c(nombre, round(resultado$statistic, 4), round(z, 4), round(p, 4), Normalidad))
})

# Convertir a marco de datos (data frame).
tabla_resultados <- as.data.frame(do.call(rbind, resultados), stringsAsFactors = FALSE)
colnames(tabla_resultados) <- c("Variable", "W'", "z'", "p", "Normalidad")

# Mostrar tabla.
cat("\nTabla. Prueba de normalidad univariante de Shapiro y Wilk con la estandarización de Royston para las", 2^k-1, "combinaciones lineales de variables\n")
print(tabla_resultados, row.names = FALSE)
cat("\nNota. Variable = suma de la correspondiente combinación sin repetición de las variables; W’ = estadística W' de Shapiro y Francia, z' = valor estandarizado de W’ mediante la estandarización de Royston; p = valor de probabilidad a la cola derecha de una distribución normal estándar; y Normalidad univariante: Sí = p ≥ α =", alpha, ", se mantiene la hipótesis nula de normalidad, y No = p <", alpha, ", se rechaza.\n")

# Función para calcular q'.
calcular_q <- function(lista_datos) {
z_vals <- sapply(lista_datos, function(x) {
resultado_sf <- sf.test(x)
z <- (log(1 - resultado_sf$statistic) - mu) / sigma
return(as.numeric(z))
})
z_values <- ifelse(z_vals < 0, 0, z_vals)
a <- sum(z_vals < 0)
q <- sum(z_values^2)
return(list(q = q, a = a, z_vals = z_vals, z_values = z_values))
}

# Valor de probabilidad con la distribución chi-cuadrado estándar.
resultado_q <- calcular_q(lista_x)
q_stat <- resultado_q$q
a <- resultado_q$a
nc <- length(lista_x)
# gl_chisq <- nc - a # con corrección por valores truncados.
gl_chisq <- nc # sin corrección por valores truncados.
q_crit <- qchisq(1 - alpha, df = gl_chisq)
p_value <- pchisq(q_stat, df = gl_chisq, lower.tail = FALSE)
power <- 1 - pchisq(q_crit, df = gl_chisq, ncp = q_stat, lower.tail = TRUE)

# Resultados QSFa.
cat("\nNúmero de combinaciones entre las", k, "variables originales: nc =", nc, ".\n")
cat("Número de valores z' anulados por ser negativos: a =", a, ".\n")
cat("Grados de libertad: gl =", gl_chisq, ".\n")
cat("\nEstadística de contraste: q' =", round(q_stat, 4), ".\n")
cat("\nValor de probabilidad y valor crítico desde una distribución chi-cuadrado estándar:\n")
cat("\nValor crítico o cuantil de orden", 1 - alpha," de una distribución chi-cuadrado estándar con", gl_chisq, "grados de libertad: 1-αχ²[gl] =", round(q_crit, 4), ".\n")
if (q_stat <= q_crit) {cat("La hipótesis nula de normalidad multivariante se mantiene a un nivel de significancia de", alpha, "con base en el valor crítico desde la aproximación.\n")
} else {cat("La hipótesis nula de normalidad multivariante se rechaza a un nivel de significancia de", alpha, "con base en el valor crítico desde la aproximación chi-cuadrado.\n")}
cat("Valor de probabilidad a la cola derecha en una distribución chi-cuadrado estándar con", gl_chisq, "grados de libertad: p =", round(p_value, 4), ".\n")
if (p_value < alpha) {cat("La hipótesis nula de normalidad multivariante se rechaza a un nivel de significancia de", alpha, "con base en el valor de probabilidad desde la aproximación.\n")
} else {cat("La hipótesis nula de normalidad multivariante se mantiene a un nivel de significancia de", alpha, "con base en el valor de probabilidad desde la aproximación chi-cuadrado.\n")}
cat("Potencia estadística desde la aproximación chi-cuadrado para un nivel de significancia de", alpha, ": phi =", power, ".\n")
cat("Si se rechaza la hipótesis nula, la potencia estadística debería ser superior a 0.5, preferiblemente mayor que 0.8; mientras que, si no se rechaza la hipótesis nula, debería ser inferior a 0.5, preferiblemente menor que 0.2. En otro caso, el resultado es contradictorio o cuestionable.\n")

# Contraste de la independencia serial en los valores z'.
z_w <- resultado_q$z_values
cat("\nContraste de la independencia entre los", nc, "valores de W’ logarítmicamente transformados, estandarizados y truncados (a 0 cuando son negativos), usados para calcular q'(z_w).\n")
RunsTest(z_w, y = NULL, alternative = "two.sided", exact = TRUE, correct = FALSE)
lag_max <- min(10, round(nc/5, 0)) # Regla de Hyndman y Athanasopoulos (2021) para determinar el rezago máximo en series temporales no estacionarias.
# lag_max <- 12 * (n/100)^0.25 # Regla de Schwert para determinar el rezago máximo.
ljung_box_test <- lapply(1: lag_max, function(k) Box.test(z_w, lag = k, type = "Ljung-Box", fitdf = 0))
ljung_box_tabla <- data.frame(
Rezago = 1:lag_max,
Chi_cuadrado = sapply(ljung_box_test, function(x) round(x$statistic, 4)),
gl = sapply(ljung_box_test, function(x) round(x$parameter, 4)),
p = sapply(ljung_box_test, function(x) round(x$p.value, 4))
)
cat("\nTabla. Prueba de Ljung y Box de independencia serial en la secuencia de valores z_w_prima\n")
print(ljung_box_tabla)
cat("\nNota. Rezago = intervalo entre los valores que conforman los pares de datos en la autocorrelación, X² = estadística de contraste, gl = grados de libertad y p = probabilidad a la cola derecha en la distribución chi-cuadrado estándar.\n")

# Función Bootstrap preservando las correlaciones de la muestra original.
resample_preserving_corr <- function(df, B) {
replicate(B, {
df[sample(1:nrow(df), size = nrow(df), replace = TRUE), ]
}, simplify = FALSE)
}

# Función para generar lista de combinaciones.
generar_lista_boot <- function(df) {
lista_boot <- setNames(as.list(df), paste0("c", 1:ncol(df)))
for (i in 2:ncol(df)) {
combinaciones <- combn(ncol(df), i)
for (j in 1:ncol(combinaciones)) {
idx <- combinaciones[, j]
nombre <- paste0("c", length(lista_boot) + 1)
lista_boot[[nombre]] <- rowSums(df[, idx])
}
}
return(lista_boot)}

# Simulaciones Bootstrap preservando correlaciones.
set.seed(123) # Se fija una semilla para reproducibilidad
B <- 1000 # Número de remuestreos con reposición desde la muestra normativa
boot_datasets <- resample_preserving_corr(datos_orig, B)

q_boot_vals <- sapply(boot_datasets, function(df_boot) {
lista_boot <- generar_lista_boot(df_boot)
calcular_q(lista_boot)$q
})

# Probabilidad desde la distribución Bootstrap empírica.
boot_p_value <- mean(q_boot_vals >= q_stat)
cat("\nDistribución Bootstrap (con correlación conservada) empírica o desde la muestral original:\n")
cat("Número de muestras Bootstrap: B =", B, ".\n")
cat("Valor de probabilidad Bootstrap a la cola derecha en la distribución Bootstrap empírica: p_boot =", round(boot_p_value, 4), ".\n")
if (boot_p_value < alpha) {cat("La hipótesis nula de normalidad multivariante se rechaza a un nivel de significancia de", alpha, "con base en el valor de probabilidad Bootstrap calculado mediante la distribución Bootstrap empírica.\n")
} else {cat("La hipótesis nula de normalidad multivariante se mantiene a un nivel de significancia de", alpha, "con base en el valor de probabilidad Bootstrap calculado mediante la distribución Bootstrap empírica.\n")}

# Generar una muestra normativa con probabilidades equiespaciadas y la transformación Cholesky de la matriz de correlaciones de la muestra original.
p <- seq(0.5 / n, 1 - 0.5 / n, length.out = n)
datos_indep <- matrix(NA, nrow = n, ncol = k)
set.seed(123)
for (j in 1:k) {datos_indep[, j] <- qnorm(sample(p))}
mat_cor <- cor(datos_orig)
mat_chol <- chol(mat_cor)
datos_normativos <- datos_indep %*% mat_chol

# Definir las variables-suma de todas las posibles combinaciones entre las variables originales.
lista_x_norm <- generar_combinaciones_lineales(datos_normativos)
result_q_norm <- calcular_q(lista_x_norm)

# Generar distribución Bootstrap normativa con remuestreo de filas completas.
B <- 1000 # Número de remuestreos con reposición desde la muestra normativa
q_boot_norm_vals <- replicate(B, {
indices <- sample(1:n, size = n, replace = TRUE)
muestra_boot <- datos_normativos[indices, ]
lista_boot_n <- generar_combinaciones_lineales(muestra_boot)
calcular_q(lista_boot_n)$q
})

# Valor de probabilidad y valor crítico en la distribución Bootstrap normativa.
boot_norm_p_value <- mean(q_boot_norm_vals >= q_stat)
# orden_cuantil <- 1-alpha # Sin compensación por cola derecha muy pesada.
orden_cuantil <- 1-2*alpha # Con compensación por cola derecha muy pesada.
valor_critico_boot_norm <- quantile(q_boot_norm_vals, probs = orden_cuantil)
cat("\nDistribución Bootstrap normativa (con correlación conservada):\n")
cat("Número de muestras Bootstrap: B =", B, ".\n")
cat("Valor de la estadística q' en la muestra normativa: q_mvn_boot =", result_q_norm$q, ".\n")
cat("Media de la distribución Bootstrap normativa: m(mvn_boot_dist) =", mean(q_boot_norm_vals), ".\n")
cat("Mediana de la distribución Bootstrap normativa: mdn(mvn_boot_dist) =", median(q_boot_norm_vals), ".\n")
p_value_mdn_norm <- 2 * min(mean(q_boot_vals >= median(q_boot_norm_vals)), mean(q_boot_vals <= median(q_boot_norm_vals)))
cat("Probabilidad de que la distribución Bootstrap empírica esté centrada en la mediana de la distribución Bootstrap normativa: p_2 colas =", round(p_value_mdn_norm, 4), ".\n")
if (p_value_mdn_norm < alpha) {cat("Se rechaza hipótesis nula de medianas equivalentes entre las distribuciones Bootstrap normativa y empírica en un contraste a dos colas con un nivel de significancia de", 1-orden_cuantil, ".\n")
} else {cat("Se mantiene hipótesis nula de medianas equivalentes entre las distribuciones Bootstrap normativa y empírica en un contraste a dos colas con un nivel de significancia de ", 1-orden_cuantil, ".\n")}
cat("Valor crítico Bootstrap o cuantil de orden", orden_cuantil, "de la distribución Bootstrap normativa: q'_crit =", round(valor_critico_boot_norm, 4), ".\n")
if (q_stat <= valor_critico_boot_norm) {cat("La hipótesis nula de normalidad multivariante se mantiene a un nivel de significancia de", 1-orden_cuantil, "con base en el valor crítico Bootstrap (normativo).\n")
} else {cat("La hipótesis nula de normalidad multivariante se rechaza a un nivel de significancia de", 1- orden_cuantil, "con base en el valor crítico Bootstrap (normativo).\n")}
cat("Valor de probabilidad a la cola derecha en la distribución Bootstrap normativa: p_mvn_boot =", round(boot_norm_p_value, 4), ".\n")
if (boot_norm_p_value < 1-orden_cuantil) {cat("La hipótesis nula de normalidad multivariante se rechaza a un nivel de significancia de", 1-orden_cuantil, "con base en el valor de probabilidad Bootstrap calculado mediante la distribución Bootstrap normativa.\n")
} else {cat("La hipótesis nula de normalidad multivariante se mantiene a un nivel de significancia de", 1-orden_cuantil, "con base en el valor de probabilidad Bootstrap calculado mediante la distribución Bootstrap normativa.\n")}
boot_power <- mean(q_boot_vals > valor_critico_boot_norm)
cat("Potencia estadística Bootstrap para un nivel de significancia de", 1-orden_cuantil, ": phi_boot =", boot_power, ".\n")

# Representación de las distribuciones muestral Bootstrap (empírica y normativa).
# Quitar el símbolo hashtag precedente para grabar como un archivo JPG o TIFF.
# jpeg("Hist_curva_dens2.jpg", width = 1200, height = 900, units = "px", res = 300)
# tiff("Hist_curva_dens2.tiff", width = 1200, height = 900, units = "px", res = 300)
# par(mar = c(4.5, 4.5, 0.5, 0.5), cex.axis = 0.8)

# Curva de densidad de la distribución chi-cuadrado estándar como gráfico base.
x_seq <- seq(0, max(qchisq(0.999, df = gl_chisq), max(q_boot_vals) + 1, q_stat +1, max(q_boot_norm_vals) + 1), length.out = 1000)
y_chi_sq <- dchisq(x_seq, df = gl_chisq)
plot(x_seq, y_chi_sq, type = "l", lwd = 3, col = "darkgreen",
main = "", xlab = "Valores q_prima", ylab = "Densidad",
xlim = c(0, max(qchisq(0.999, df = gl_chisq), max(q_boot_vals), q_stat +1, max(q_boot_norm_vals))),
ylim = c(0, max(density(q_boot_vals)$y, y_chi_sq)))
# Histograma amarillo de la distribución Bootstrap empírica de la estadística de contraste q'.
hist(q_boot_vals, breaks = 20, freq = FALSE, col = rgb(1, 1, 0, 0.5), border = "yellow2", add = TRUE)
# Curva amarilla de densidad de la distribución Bootstrap empírica de la estadística de contraste q'.
lines(density(q_boot_vals), col = "yellow3", lwd = 2)
# Curva roja de densidad Bootstrap normativa.
lines(density(q_boot_norm_vals), col = "red", lwd = 2)
# Línea vertical violeta para el valor de la estadística de contraste q' observado.
abline(v = q_stat, col = "purple", lwd = 2, lty = 2)
# Línea vertical roja para el valor crítico normativo Bootstrap.
abline(v = valor_critico_boot_norm, col = "red", lwd = 2, lty = 2)
# Línea vertical verde para el valor crítico de la distribución chi-cuadrado estándar.
abline(v = qchisq(1-alpha, df = gl_chisq), col = "darkgreen", lwd = 2, lty = 2)
# dev.off() # Quitar el símbolo hashtag cuando se graba la figura.

cat("\nFigura 2. Histograma (amarillo) con las curvas de densidad de las distribuciones Bootstrap empírica (amarilla), Bootstrap normativa (roja) y chi-cuadrado estándar (verde). El valor de la estadística de contraste q se representa mediante una línea vertical violeta (q' =", q_stat, "), el valor crítico de la distribución Bootstrap normativa mediante una línea vertical roja (", orden_cuantil, "_q_mvn_boot =", valor_critico_boot_norm, ") y el valor crítico de la distribución chi-cuadrado estándar con", gl_chisq, "grados de libertad mediante una línea vertical verde (", 1-alpha, "χ²[", gl_chisq, "] =", q_crit, ").\n")
cat("\nNota. Histograma amarillo (regla de la Universidad de Rice) = densidades de la distribución Bootstrap empírica de la estadística de contraste q',
curva amarilla = curva de densidad de la distribución Bootstrap empírica de la estadística de contraste q',
curva roja = curva de densidad de la distribución Bootstrap normativa (generada desde variables normalmente distribuidas con la misma estructura correlacional que las originales) de la estadística de contraste q',
curva verde = distribución chi-cuadrado estándar con", gl_chisq, "grados de libertad,
línea vertical amarilla = valor de la estadística de contraste q' en la muestra original,
línea vertical roja = valor crítico o cuantil de orden", orden_cuantil, "de la distribución Bootstrap normativa,
línea vertical verde = valor crítico o cuantil de orden", 1-alpha, "de una distribución chi-2 estándar con", gl_chisq, "grados de libertad.\n")

#------------------------------------------------------------------------------------------------------ 
# Apéndice 3
# Prueba H de normalidad multivariante de Royston (1983) para muestras de tamaño de 10 a 2000 con k mediciones.
#-------------------------------------------------------------------------------------------------------

# Paquete requerido.
library(MVN)
# El paquete MVN requiere R ≥ 3.5.0
# Si no se tiene el paquete, seleccionar espejo más próximo: chooseCRANmirror()
# install.packages("MVN")

# Lista de variables originales.
x1 <- c(7, 6, 4, 5, 3, 9, 5, 3, 3, 5, 5, 6, 5, 9, 6, 10, 8, 7, 7, 3, 6, 5, 8, 7, 1, 4, 7, 2, 6, 4, 4, 8, 6, 2, 1, 4)
x2 <- c(6, 7, 4, 4, 5, 7, 6, 5, 3, 4, 7, 5, 1, 8, 5, 10, 8, 6, 8, 2, 4, 2, 7, 9, 3, 3, 5, 3, 5, 4, 6, 7, 6, 6, 1, 4)
x3 <- c(7, 9, 5, 5, 6, 9, 5, 5, 6, 2, 6, 7, 2, 10, 6, 9, 9, 3, 8, 2, 4, 4, 9, 6, 2, 3, 7, 4, 3, 6, 4, 7, 6, 5, 1, 4)
x4 <- c(8, 7, 6, 7, 3, 9, 4, 2, 3, 4, 5, 5, 1, 10, 6, 8, 9, 2, 6, 2, 6, 4, 7, 7, 5, 5, 7, 5, 5, 4, 6, 6, 8, 3, 1, 4)
x5 <- c(6, 7, 4, 6, 3, 6, 8, 2, 4, 4, 8, 8, 4, 10, 6, 8, 7, 4, 6, 2, 2, 4, 9, 9, 3, 4, 4, 6, 4, 6, 6, 7, 7, 3, 1, 3)
x6 <- c(8, 6, 4, 7, 5, 9, 8, 3, 5, 7, 6, 6, 5, 10, 4, 10, 9, 4, 6, 3, 6, 2, 8, 8, 5, 3, 7, 7, 4, 5, 6, 7, 5, 3, 2, 4)
datos_orig <- data.frame(x1, x2, x3, x4, x5, x6)

n <- length(x1)
k <- length(datos_orig)
alpha <- 0.05
result <- mvn(datos_orig, subset = NULL, mvnTest = "royston")

# Potencia estadística para H.
R <- cor(datos_orig)
lambda <- 5
mu <- 0.715
ln_n <- log(n)
v_n <- 0.21364 + 0.015124 * ln_n^2 - 0.0018034 * ln_n^3
r_ij <- R[lower.tri(R)]
c_ij <- r_ij^lambda * (1 - (mu / v_n) * (1 - r_ij)^mu)
c_bar <- mean(c_ij)
df_H <- k / (1 + (k - 1) * c_bar)
power_H_5 <- 1 - pchisq(qchisq(1-alpha, df = df_H), df = df_H, ncp = result$multivariateNormality$H, lower.tail = TRUE)
power_H_10 <- 1 - pchisq(qchisq(1-2*alpha, df = df_H), df = df_H, ncp = result$multivariateNormality$H, lower.tail = TRUE)

# Mostrar resultados.
cat("\nTamaño muestral: n =", n, ".\n")
cat("Número de variables en la muestra original: k =", k, ".\n")
cat("Media de las correlaciones entre las variables: m(R) =", round(mean(r_ij), 4), ".\n")
cat("\nPrueba de normalidad H de Royston (1983)\n")
cat("Valor de la estadística de contraste: H =", round(result$multivariateNormality$H, 4), ".\n")
cat("Grados de libertad: gl =", df_H, ".\n")
cat("Valor de probabilidad para la hipótesis nula de normalidade multivariante: p =", result$multivariateNormality$p, ".\n")
cat("\nPotencia para H para un nivel de significancia de", alpha, ": phi =", round(power_H_5, 4), ".\n")
cat("Potencia para H para un nivel de significancia de", 2*alpha, ": phi =", round(power_H_10, 4), ".\n")
cat("Si se rechaza la hipótesis nula, la potencia estadística debería ser superior a 0.5, preferiblemente mayor que 0.8; mientras que, si no se rechaza la hipótesis nula, debería ser inferior a 0.5, preferiblemente menor que 0.2. En otro caso, el resultado es contradictorio o cuestionable.\n")

#--------------------------------------------------------------------------------------------------
# Potencia relativa entre las pruebas H (Royston, 1983) y Q, evaluada mediante la aproximación chi-cuadrado y el enfoque Bootstrap.
#---------------------------------------------------------------------------------------------------
# Para ejecutar correctamente este script, es necesario correr previamente, y de forma secuencial, el script del Apéndice 1 (estadísticas W de Shapiro y Wilk, 1965, estandarizadas por Royston, 1992) o el del Apéndice 2 (estadísticas W’ de Shapiro y Francia, 1972, estandarizadas por Royston, 1993).
# Dichos scripts generan los valores de: power (potencia bajo la aproximación chi-cuadrado) y boot_power (potencia obtenida mediante Bootstrap)
# Después de ello, debe ejecutarse el script del Apéndice 3, que calcula power_H para la prueba H.

# Alternativamente, es posible introducir manualmente los valores obtenidos en una ejecución previa e independiente de los Apéndices 1 o 2. Los valores en azul son los de la muestra analizada de 36 6-tuplas.
# power <- 0.07598624 # Potencia de la aproximación chi-cuadrado de QSW con α = 0.05.
# boot_power <- 0.027 # Potencia según el enfoque Bootstrap de QSW con α = 0.1.
# power <- 0.06032782 # Potencia de la aproximación chi-cuadrado de QSF con α = 0.05.
# boot_power <- 0.018 # Potencia según el enfoque Bootstrap de QSF con α = 0.1.

relative_power1 <- power_H_5 / power
cat("\nPotencia relativa (H / Q_aprox_Chi2) con un nivel de significancia de", alpha, ": Rel_phi =", round(relative_power1, 4), ".\n")
cat("Si se rechaza la hipótesis nula de normalidad multivariante (p < α = 0.05) y los valores de potencia son mayores que 0.5, valores de potencia relativa H / Q_aprox_Chi2 > 1 indican mayor potencia de la prueba H de Royston y < 1 indican mayor potencia de la prueba Q desde la aproximación chi-cuadrado.\n")
cat("Si se mantiene la hipótesis nula de normalidad multivariante (p ≥ α = 0.05)  y los valores de potencia son menores que 0.5, valores de potencia relativa H / Q_aprox_Chi2 > 1 sugieren un mejor rendimiento de la prueba Q desde la aproximación chi-cuadrado y < 1 mejor rendimiento de la prueba H.\n")

relative_power2 <- power_H_10 / boot_power
cat("\nPotencia relativa (H / Q_boot) con un nivel de significancia de", 2*alpha, ": Rel_phi =", round(relative_power2, 4), ".\n")
cat("Si se rechaza la hipótesis nula de normalidad multivariante (p < α = 0.1) y los valores de potencia son mayores que 0.5, valores de potencia relativa H / Q_boot > 1 indican mayor potencia de la prueba H de Royston y < 1 mayor potencia de la prueba Q desde el enfoque Bootstrap.\n")
cat("Si se mantiene la hipótesis nula de normalidad multivariante (p ≥ α = 0.1) y los valores de potencia son menores que 0.5, valores de potencia relativa H / Q_boot > 1 sugieren un mejor rendimiento de la prueba Q desde el enfoque Bootstrap y < 1 mejor rendimiento de la prueba H.\n")



