library(urca)
library(dynlm)
library(lmtest)
library(forecast)
library(tseries)
library(tsoutliers)


################################################################################
# 1. CARGA Y PREPARACIÓN DE DATOS NOAA (CO2 MENSUAL)                           #
################################################################################

# Leemos el archivo ignorando las líneas que empiezan por "#"
datos_co2 <- read.csv("co2_mm_gl.csv", 
                         comment.char = "#", 
                         header = TRUE)

# Filtramos desde 2016 en adelante (para coincidir con la electricidad)
datos_co2_reciente <- subset(datos_co2, year >= 2016)

ts_co2 <- ts(datos_co2$average, 
             start = c(1979, 1), 
             frequency = 12)

plot(ts_co2, main="CO2 Mensual (NOAA) - ppm", col="red")


# Creamos la Serie Temporal (TS) mensual
# Usamos la columna 'average' (dato real medido)
par(mfrow=c(1,1))
ts_co2_mensual <- ts(datos_co2_reciente$average, 
                     start = c(2016, 1), 
                     frequency = 12)



# CARGA Y PREPARACIÓN DE DATOS DE ELECTRICIDAD Italiana (MENSUAL)              #
# Cargamos y seleccionamos columnas
load_electricity <- read.csv("datos_electricidad/datos_electricidad_mensual.csv")

# Nos queremos quedar con los datos de ITALIA, por lo tanto seleccionamos los datos cuya columna 'CountryCode' sea 'IT'
load_it <- subset(load_electricity, CountryCode == "IT")


# Creamos la Serie Temporal
load_it_ts <- ts(load_it[, 3], 
                 start = c(2016, 1), 
                 frequency = 12)

# Graficamos
plot(load_it_ts, 
     main = "Demanda Eléctrica Italia (Mensual)",
     xlab = "Año", ylab = "MWh", col = "blue")
grid()

plot(decompose(load_it_ts))

################################################################################
# 2.1. ESTUDIO DE OUTLIERS E INTERVENCIONES                                    #
################################################################################

# Queremos estudiar los outliers tanto aditivos (AO) como cambios de nivel (LS)
# Empecemos con los outliers aditivos
outliers.load_it.ao <- tso(load_it_ts, types = c("AO"), maxit.iloop = 10)

# Veamos los resultados
print(outliers.load_it.ao)
plot(outliers.load_it.ao)

# RESULTADO DEL ANÁLISIS DE OUTLIERS (ITALIA)
# --------------------------------------------------------------------------
# Outliers:
#   type ind    time    coefhat   tstat
# 1   AO  36 2018:12   22685712   38.73  # ERROR DE DATOS (Data Entry Error): Anomalía inverosímil.
#                                        # Indica un aumento de +22 TWh (casi el 100% de la demanda media).
#                                        # Técnicamente imposible. Se debe corregir obligatoriamente.
# --------

# Vamos a corregirlo con la finalidad de que no enmascare otros posibles outliers
# aditivos que pueda presentar la serie temporal.

tabla_outliers <- outliers.load_it.ao$outliers
ind_tiempo <- tabla_outliers$ind[1]
impacto <- tabla_outliers$coefhat[1]

# Restamos el impacto (el "pico" extra) a la serie original
load_it_ts[ind_tiempo] <- load_it_ts[ind_tiempo] - impacto

plot(load_it_ts, 
     main = "Demanda Eléctrica Italia (Mensual)",
     xlab = "Año", ylab = "MWh", col = "blue")


# Visualización rápida
par(mfrow=c(2,1))
plot(ts_co2_mensual, main="CO2 Mensual (NOAA) - ppm", col="red")
plot(load_it_ts, main="Consumo Eléctrico Italia - MWh", col="blue")


# Ahora tenemos que ver que ambas series sean estacionarias.
# Para ello, vamos a diferenciarlas. Hemos visto que la serie eléctrica necesita
# una diferencia regular y una estacional. 

load_it_sdiff <- diff(load_it_ts, differences=1, lag=12)
load_it_sdiff_diff <- diff(load_it_sdiff, differences=1, lag=1)

# Visualizamos la serie diferenciada
par(mfrow=c(1,1))
plot(load_it_sdiff_diff, 
     main = paste("Serie Estacionaria Final (d=", 1, ", D=", 1, ")", sep=""),
     ylab = "Variación de Demanda (MWh)")
grid()
abline(h=0, col="red", lty=2) 

test_adf <- adf.test(load_it_sdiff_diff)
print(test_adf)

# El pvalor es más pequeño a 0.1 => Rechazamos H0 => La serie es estacionaria

# Veamos lo mismo en la serie del CO2, 

numero_D <- nsdiffs(ts_co2_mensual)
cat("Diferencias Estacionales (D) necesarias:", numero_D, "\n")
co2_sdiff <- diff(ts_co2_mensual, differences=1, lag=12)
numero_d <- ndiffs(co2_sdiff)
cat("Diferencias Regulares (d) necesarias:", numero_d, "\n")
co2_sdiff_diff <- diff(co2_sdiff, differences=1, lag=1)

# Visualizamos la serie diferenciada
par(mfrow=c(1,1))
plot(co2_sdiff_diff, 
     main = paste("Serie Estacionaria Final CO2 (d=", numero_d, ", D=", numero_D, ")", sep=""),
     ylab = "Variación de CO2 (ppm)")
grid()
abline(h=0, col="red", lty=2)

test_adf_co2 <- adf.test(co2_sdiff_diff)
print(test_adf_co2)

# El pvalor es 0.02398 < 0.05 => Rechazamos H0 => La serie es estacionaria

# Podemos entonces empezar con la cointegración de las series.

################################################################################
# COINTEGRACIÓN (procedimiento tipo: Engle–Granger + ECM + Durbin–Watson)
# En este planteamiento:
#   M1.ts = CO2 (serie estacionaria tras d=1 y D=1)
#   M2.ts = Electricidad (serie estacionaria tras d=1 y D=1)
################################################################################

# 1) SINCRONIZACIÓN: nos quedamos solo con los meses comunes entre ambas series
#    (evita desalineaciones temporales y asegura misma longitud)
datos_coint <- ts.intersect(
  CO2  = co2_sdiff_diff,       # CO2 transformado (estacionario)
  Elec = load_it_sdiff_diff    # Electricidad transformada (estacionaria)
)

# 2) ASIGNACIÓN DE VARIABLES: definimos claramente qué es M1 (dependiente) y M2
M1.ts <- datos_coint[, "CO2"]   # Variable dependiente del modelo de largo plazo (CO2)
M2.ts <- datos_coint[, "Elec"]  # Variable explicativa (Electricidad)

cat("Tamaño de la muestra mensual:", length(M1.ts), "meses.\n")

################################################################################
# 3) COMPROBACIÓN (opcional) DE ESTACIONARIEDAD EN LAS SERIES USADAS
#    Como ya las hemos diferenciado (d=1, D=1), deberían ser I(0).
#    Aquí lo confirmamos con ADF (ur.df) sobre las series ya transformadas.
################################################################################

# ADF sin constante ni tendencia (type="none") porque estas series ya oscilan
# alrededor de 0 tras diferenciar (no tienen nivel medio ni tendencia)
summary(ur.df(M1.ts, type = "none", selectlags = "AIC"))
summary(ur.df(M2.ts, type = "none", selectlags = "AIC"))

################################################################################
# 4) REGRESIÓN DINÁMICA
#    Modelo: CO2_t depende de:
#      - su propio retardo (captura inercia/persistencia)
#      - electricidad actual
#      - electricidad con retardo (efectos con retraso)
#
#    CO2_t = α + β1*CO2_{t-1} + β2*Elec_t + β3*Elec_{t-1} + u_t
################################################################################

lr.regkd <- dynlm(M1.ts ~ L(M1.ts) + M2.ts + L(M2.ts))

# Inspección del modelo: significación de retardos y ajuste global
summary(lr.regkd)

# 5) RESIDUOS: u_t = parte de CO2_t no explicada por la ecuación anterior.
#    En el enfoque Engle–Granger, lo importante es comprobar si estos residuos
#    son estacionarios (u_t ~ I(0)).
error <- residuals(lr.regkd)

# Visualización rápida de los residuos (deberían fluctuar alrededor de 0)
plot(error, main = "Residuos de la regresión dinámica", ylab = "u_t")
abline(h = 0, col = "red", lty = 2)
grid()

################################################################################
# 5) TEST ADF SOBRE RESIDUOS (ENGLE–GRANGER)
#    Hipótesis del ADF:
#      H0: u_t tiene raíz unitaria (no estacionario) -> NO hay equilibrio estable
#      H1: u_t es estacionario -> existe relación estable (cointegración)
################################################################################

adf.test(error)

# Interpretación práctica:
# - p-valor pequeño (< 0.05): rechazamos H0 -> residuos estacionarios -> evidencia
#   compatible con cointegración.
# - p-valor grande: no rechazamos H0 -> no hay evidencia de cointegración.

################################################################################
# 6) MODELO DE CORRECCIÓN DEL ERROR (ECM)
#    Si hay cointegración, el ECM describe el ajuste de corto plazo hacia el
#    equilibrio de largo plazo.
#
#    ΔCO2_t = c + γ*u_{t-1} + δ*ΔElec_{t-1} + φ*ΔCO2_{t-1} + e_t
#
#    Donde u_{t-1} = L(error) es el desequilibrio del periodo anterior.
################################################################################

# Construimos un objeto con todas las series necesarias, sincronizadas
data_ecm <- ts.intersect(M1.ts, M2.ts, error)

# Estimamos el ECM:
# - d(M1.ts): variación mensual de CO2 (corto plazo)
# - L(error): término de corrección (desequilibrio pasado)
# - L(d(M2.ts)): variación pasada de electricidad
# - L(d(M1.ts)): inercia en la variación de CO2
ecm <- dynlm(d(M1.ts) ~ L(error) + L(d(M2.ts)) + L(d(M1.ts)), data = data_ecm)

summary(ecm)

# Interpretación clave:
# - Si γ (coef. de L(error)) es NEGATIVO y significativo:
#   hay evidencia de ajuste hacia el equilibrio (corrección del desequilibrio).
# - Si no es significativo: no se detecta mecanismo de ajuste.

################################################################################
# 7) DURBIN–WATSON EN EL ECM
#    Comprueba autocorrelación de primer orden en los residuos del ECM.
#    H0: no autocorrelación
#    Si no rechazamos H0 (p-valor alto), el ECM está “bien comportado” en ese
#    sentido.
################################################################################

dwtest(ecm)

