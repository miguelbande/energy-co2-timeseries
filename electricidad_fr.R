################################################################################
#                TFM - SERIE DEMANDA ELÉCTRICA - FRANCIA                       #
################################################################################

library(haven)
library(readxl)
library(datasets)
library(timsac)
library(forecast)
library(descomponer)
library(lmtest)
library(nortest)
library(tsoutliers)
library(expsmooth)
library(fma)
library(tseries)
library(tsoutliers)


################################################################################
# 1. CARGA DE DATOS Y ANÁLISIS EXPLORATORIO                                    #
################################################################################

# Cargamos y seleccionamos columnas
load_electricity <- read.csv("datos_electricidad/datos_electricidad_mensual.csv")

# Nos queremos quedar con los datos de Francia, por lo tanto seleccionamos los datos cuya columna 'CountryCode' sea 'FR'
load_fr <- subset(load_electricity, CountryCode == "FR")


# Creamos la Serie Temporal
load_fr_ts <- ts(load_fr[, 3], 
                 start = c(2016, 1), 
                 frequency = 12)

# Graficamos
plot(load_fr_ts, 
     main = "Demanda Eléctrica Francia (Mensual)",
     xlab = "Año", ylab = "MWh", col = "blue")
grid()

plot(decompose(load_fr_ts))

################################################################################
# 2.1. ESTUDIO DE OUTLIERS E INTERVENCIONES                                    #
################################################################################

# Queremos estudiar los outliers tanto aditivos (AO) como cambios de nivel (LS)
# Empecemos con los outliers aditivos
outliers.load_fr.ao <- tso(load_fr_ts, types = c("AO"), maxit.iloop = 10)

# Veamos los resultados
print(outliers.load_fr.ao)
plot(outliers.load_fr.ao)

# RESULTADO:
# --------------------------------------------------------------------------
# Outliers:
#   type ind    time   coefhat   tstat
# 1   AO  26 2018:02   4411654   5.084  # OLA DE FRÍO Pico extremo de calefacción eléctrica.
# 2   AO  48 2019:12  -3358731  -3.976  # CLIMA Menor demanda de lo esperado.
# 3   AO  58 2020:10   3707920   4.447  # CLIMA Ola de frío temprana.
# 4   AO  85 2023:01  -2971968  -3.351  # CRISIS ENERGÉTICA: Ahorro masivo voluntario por crisis nuclear/gas.

# Son todos generalmente en meses en los que la temperatura es un factor crucial.
# --------------------------------------------------------------------------

# Ahora veamos los cambios de nivel (LS)
outliers.load_fr.ls <- tso(load_fr_ts, types = c("LS"), maxit.iloop = 10)

# Vemos los resultados
print(outliers.load_fr.ls)
plot(outliers.load_fr.ls)

# RESULTADO DEL ANÁLISIS DE OUTLIERS (FRANCIA) - CAMBIOS DE NIVEL (LS)
# --------------------------------------------------------------------------
# Outliers:
#   type ind    time   coefhat   tstat
# 1   LS  26 2018:02   4722767   4.365  
# 2   LS  28 2018:04  -6265565  -5.846  
# 3   LS  82 2022:10  -2633654  -7.566 
# --------------------------------------------------------------------------


################################################################################
# 2.2. GESTIÓN DE INTERVENCIONES (PREPARACIÓN PARA MODELAR)                    #
################################################################################

# ESTRATEGIA: 
# 1. Vamos a corregir los outliers ADITIVOS.

# 2. Los LS son cambios de nivel -> LO DEJAMOS (No se corrige).
#    Pues las predicciones futuras han de partir de ese escalón más alto.

# Hacemos una copia de la serie original
serie_para_modelar <- load_fr_ts

# Sacamos la tabla de intervenciones detectadas
tabla_outliers <- outliers.load_fr.ao$outliers

# Bucle inteligente: Recorre los outliers y borra SOLO los tipos "AO"
for(i in 1:nrow(tabla_outliers)) {
  
  if(tabla_outliers$type[i] == "AO") {
    # Obtenemos posición e impacto
    ind_tiempo <- tabla_outliers$ind[i]
    impacto <- tabla_outliers$coefhat[i]
    
    # Restamos el impacto para limpiar el dato (rellenar el hueco)
    serie_para_modelar[ind_tiempo] <- serie_para_modelar[ind_tiempo] - impacto
  }
}

# Comprobación visual final
ts.plot(load_fr_ts, serie_para_modelar,
        col = c("red", "blue"),
        lwd = c(1, 2),
        main = "Comparativa: Realidad vs. Serie para Modelar",
        ylab = "MWh")

legend("topright", 
       legend=c("Serie Original (Con Intervenciones)", "Serie Limpia de AO (Mantiene LS)"),
       col=c("red", "blue"), lty=1, lwd=c(1,2))



################################################################################
# 3. TRANSFORMACIONES DE LA SERIE                                              #
################################################################################

# 1. Cálculo de Lambda óptimo para Box-Cox
lambda_optimo <- BoxCox.lambda(serie_para_modelar)

print(paste("El valor de Lambda óptimo es:", round(lambda_optimo, 4)))
# "El valor de Lambda óptimo es: -0.687"

# 2. Aplicamos la transformación
serie_boxcox <- BoxCox(serie_para_modelar, lambda_optimo)

# 3. Gráfico Comparativo: Original vs Transformada

# Gráfico 1: Original
plot(serie_para_modelar, 
     main = "Serie Original (Limpia de AO)",
     ylab = "Mt CO2", col = "black")
grid()

# Gráfico 2: Transformada
plot(serie_boxcox, 
     main = paste("Transformada Box-Cox (lambda=", round(lambda_optimo, 3), ")", sep=""),
     ylab = "Valor Transformado", col = "blue")
grid()

# Volvemos a pantalla normal
par(mfrow = c(1, 1))



################################################################################
# 3. ESTACIONARIEDAD: CÁLCULO DE DIFERENCIAS (D y d)                           #
################################################################################

# IMPORTANTE:
# Al ser datos MENSUALES (freq=12), la estrategia cambia (respecto al CO2).
# 1º Miramos si hay que quitar estacionalidad (Diferencia Estacional D).
# 2º Miramos si, tras quitar la estacionalidad, queda tendencia (Diferencia Regular d).

# 'serie_para_modelar' es la serie de electricidad (original)

# --- PASO 1: DIFERENCIA ESTACIONAL (D) ---
# nsdiffs chequea si la serie necesita diferenciarse estacionalmente (lag 12)
# Test OCSB o SH (Canova-Hansen) son habituales.
numero_D <- nsdiffs(serie_boxcox)
cat("Diferencias Estacionales (D) necesarias:", numero_D, "\n")

# Diferencias Estacionales (D) necesarias: 1 

# Aplicamos la diferencia estacional si hace falta
if (numero_D > 0) {
  # lag=12 significa: Dato(Enero 2017) - Dato(Enero 2016)
  serie_desestacionalizada <- diff(serie_boxcox, lag = 12, differences = numero_D)
  
  plot(serie_desestacionalizada, main="Serie tras Diferencia Estacional")
} else {
  serie_desestacionalizada <- serie_boxcox
  cat("No ha sido necesaria la diferencia estacional.\n")
}

test_adf <- adf.test(serie_desestacionalizada)
print(test_adf)

# El p-valor es 0.08411 > 0.05 => No podemos rechazar H0 => La serie NO es estacionaria aún

# --- PASO 2: DIFERENCIA REGULAR (d) ---
# Ahora cogemos la serie que ha salido del paso 1 y miramos si tiene tendencia
numero_d <- ndiffs(serie_desestacionalizada)

cat("Diferencias Regulares (d) necesarias:", numero_d, "\n")
# Diferencias Regulares (d) necesarias: 0 

# Vamos a forzar una diferencia regular, pues sale que no es necesaria ninguna
numero_d <- 1

# Aplicamos la diferencia regular si hace falta
if (numero_d > 0) {
  # lag=1 significa: Dato(Feb) - Dato(Ene)
  serie_estacionaria <- diff(serie_desestacionalizada, lag = 1, differences = numero_d)
} else {
  serie_estacionaria <- serie_desestacionalizada
}

# --- RESUMEN FINAL ---
cat("Orden total de diferenciación: D =", numero_D, "| d =", numero_d, "\n")


# 3. Graficamos la serie FINAL (Estacionaria)
# Deberías ver "ruido" oscilando alrededor de 0, sin patrones claros ni tendencia.
plot(serie_estacionaria, 
     main = paste("Serie Estacionaria Final (d=", numero_d, ", D=", numero_D, ")", sep=""),
     ylab = "Variación de Demanda (MWh)")
grid()
abline(h=0, col="red", lty=2) 


# 4. Comprobación matemática: Test de Dickey-Fuller Aumentado (ADF)
# H0 (Nula): La serie NO es estacionaria
# H1 (Alternativa): La serie SÍ es estacionaria (p-value < 0.05)
test_adf <- adf.test(serie_estacionaria)
print(test_adf)

# El pvalor es inferior a 0.01 => Rechazamos H0 => La serie es estacionaria


# 5. Correlogramas (ACF y PACF)
# Ajustamos lag.max a 48 (4 años) para ver si hay "ecos" estacionales cada 12 meses
par(mfrow=c(1,2))
acf(serie_estacionaria, lag.max = 48, main="ACF (MA - q / Q)")
pacf(serie_estacionaria, lag.max = 48, main="PACF (AR - p / P)")
par(mfrow=c(1,1))


################################################################################
# 4. FUERZA BRUTA                                                              #
################################################################################


fuerza_bruta_arima <- function(serie, 
                               d = NULL, 
                               D = NULL, 
                               max_p = 3, max_q = 3, 
                               max_P = 2, max_Q = 2) {
  
  # 1. PARAMETROS BÁSICOS
  freq_data <- frequency(serie)
  n <- length(serie) # Longitud de la serie
  es_estacional <- freq_data > 1
  
  cat("------------------------------------------------------\n")
  cat("INICIANDO FUERZA BRUTA \n")
  cat("Datos: n =", n, "| Frecuencia =", freq_data, "\n")
  
  
  # 2. GENERAR GRID
  if (!es_estacional) {
    # ANUAL
    cat("Modo: ANUAL (No Estacional)\n")
    
    # En caso de que no se haya introducido el número de diferencias se precalcula
    if(is.null(d)) {
      d <- ndiffs(serie)
      cat("d calculado:", d, "\n")
    }
    
    combinaciones <- expand.grid(p = 0:max_p, d = d, q = 0:max_q)
    combinaciones <- subset(combinaciones, !(p == 0 & q == 0)) # Filtro p=0,q=0
    
    # Pre-calculamos el lag base (sin contar df del modelo aún)
    # min(10, n/5)
    lag_base <- min(10, n/5)
    
  } else {
    # ESTACIONAL
    cat("Modo: ESTACIONAL\n")
    if(is.null(D)) {
      D <- nsdiffs(serie)
      cat("D calculado:", D, "\n")
    }
    serie_desestacionalizada >- diff(serie, lag = freq_data, differences = D)
    # En caso de que no se haya introducido el número de diferencias se precalcula
    if(is.null(d)) {
      d <- ndiffs(serie)
      cat("d calculado:", d, "\n")
    }
    
    combinaciones <- expand.grid(p = 0:max_p, d = d, q = 0:max_q,
                                 P = 0:max_P, D = D, Q = 0:max_Q)
    combinaciones <- subset(combinaciones, !(P == 0 & Q == 0)) # Filtro estacional nulo
    
    # Pre-calculamos el lag base
    # min(2m, n/5) -> m es freq_data
    lag_base <- min(2 * freq_data, n/5)
  }
  
  resultados <- data.frame()
  total_modelos <- nrow(combinaciones)
  
  # 3. BUCLE DE MODELADO
  for (i in 1:total_modelos) {
    params <- combinaciones[i, ]
    
    tryCatch({
      # A) Ajustar modelo
      if (!es_estacional) {
        modelo <- arima(serie, order = c(params$p, params$d, params$q))
        orden_texto <- paste0("ARIMA(", params$p, ",", params$d, ",", params$q, ")")
        
        # Grados de libertad gastados (p + q)
        n_params <- params$p + params$q
        
      } else {
        modelo <- arima(serie, 
                        order = c(params$p, params$d, params$q),
                        seasonal = list(order = c(params$P, params$D, params$Q), period = freq_data))
        orden_texto <- paste0("ARIMA(", params$p, ",", params$d, ",", params$q, ")(", params$P, ",", params$D, ",", params$Q, ")")
        
        # Grados de libertad gastados (p + q + P + Q)
        n_params <- params$p + params$q + params$P + params$Q
      }
      
      # B) CÁLCULO DINÁMICO DEL LAG (Regla de la función checkresiduals)
      # "Constrained to be at least df + 3"
      lag_final <- max(lag_base, n_params + 3)
      
      # Redondeamos por si acaso n/5 dio decimales
      lag_final <- round(lag_final)
      
      # C) TESTS
      
      # 1. Coeficientes
      test_coef <- coeftest(modelo)
      if(nrow(test_coef) > 0) {
        coef_significativos <- all(test_coef[, 4] < 0.05, na.rm=TRUE)
      } else {
        coef_significativos <- TRUE 
      }
      
      # 2. Ljung-Box (USANDO fitdf)
      # fitdf = n_params le dice al test que hemos estimado parametros
      # Esto ajusta el p-valor para ser estadísticamente correcto
      test_lb <- Box.test(modelo$residuals, lag = lag_final, type = "Ljung-Box", fitdf = n_params)
      
      # 3. Normalidad (tryCatch por si falla con pocos datos)
      test_norm_pval <- tryCatch({ lillie.test(modelo$residuals)$p.value }, error = function(e) 0)
      
      # 4. Homocedasticidad
      res_num <- as.numeric(modelo$residuals)
      tiempo <- 1:length(res_num)
      test_bp <- bptest(res_num ~ tiempo)
      
      # D) GUARDAR
      resultados <- rbind(resultados, data.frame(
        ID = i,
        Modelo = orden_texto,
        AIC = modelo$aic,
        Signif_Coef = coef_significativos,
        Indep_LB = test_lb$p.value > 0.05,
        Norm_Lillie = test_norm_pval > 0.05,
        Homo_BP = test_bp$p.value > 0.05,
        Pval_LB = round(test_lb$p.value, 3),
        Pval_Lill = round(test_norm_pval,3),
        Pval_BP = round(test_bp$p.value, 3),
        Lag_Usado = lag_final # Guardamos qué lag se usó para verificar
      ))
      
    }, error = function(e) { return(NULL) })
    
    if(i %% 20 == 0 || i == total_modelos) cat(sprintf("\rProcesando... %d/%d", i, total_modelos))
  }
  
  cat("\n--- FIN ---\n")
  
  if (nrow(resultados) > 0) {
    mejores <- subset(resultados, Signif_Coef & Indep_LB & Homo_BP)
    return(list(
      ranking_completo = resultados[order(resultados$AIC), ],
      mejores_candidatos = mejores[order(mejores$AIC), ]
    ))
  } else {
    return(NULL)
  }
}


################################################################################
# 6. TORNEO DE VALIDACIÓN: TRAIN (80%) vs TEST (20%) CON BOX-COX               #
################################################################################

# 1. DIVISIÓN DE LA SERIE (80% - 20%)
# -----------------------------------
n_total <- length(serie_para_modelar)
n_train <- floor(0.80 * n_total) 

# Usamos la función subset que respeta las fechas ts
train_set <- subset(serie_para_modelar, end = n_train)
valid_set <- subset(serie_para_modelar, start = n_train + 1)
h_valid   <- length(valid_set)

cat("Datos Totales:", n_total, "meses.\n")
cat("Entrenamiento (80%):", length(train_set), "meses.\n")
cat("Validación (20%):   ", h_valid, "meses.\n")

# --- PASO CLAVE: CÁLCULO DE LAMBDA SOLO EN TRAIN ---
# Calculamos lambda viendo solo el pasado para no hacer trampas
lambda_train <- BoxCox.lambda(train_set)
cat("Lambda calculado sobre Train:", round(lambda_train, 4), "\n")


# 2. ENTRENAMIENTO DEL CANDIDATO 1: SARIMA (Fuerza Bruta en Train)
# ----------------------------------------------------------------
cat("\n--- BUSCANDO MEJOR SARIMA PARA EL CONJUNTO DE TRAIN ---\n")

# IMPORTANTE: A la fuerza bruta le pasamos la serie YA transformada
# para que busque patrones lineales sobre los datos estabilizados.
train_transformed <- BoxCox(train_set, lambda = lambda_train)

res_sarima_train <- fuerza_bruta_arima(train_transformed, 
                                       d = 1,
                                       D = 1, 
                                       max_p = 3, max_q = 3, 
                                       max_P = 3, max_Q = 3) 

# Seleccionamos el mejor del train
candidatos_norm <- subset(res_sarima_train$mejores_candidatos, Norm_Lillie == TRUE)
print(candidatos_norm)

# NO SE ENCUENTRAN CANDIDATOS QUE CUMPLAN NORMALIDAD

if(nrow(candidatos_norm) > 0) {
  modelo_sarima_texto <- candidatos_norm$Modelo[1]
} else {
  modelo_sarima_texto <- res_sarima_train$mejores_candidatos$Modelo[1]
}

cat("Mejor SARIMA encontrado para el Train:", modelo_sarima_texto, "\n")
# Mejor SARIMA encontrado para el Train: ARIMA(1,1,1)(1,1,1) 




################################################################################
# 6. TORNEO DE VALIDACIÓN: SARIMA vs SUAVIZADO (80% / 20%)               #
################################################################################

# 2. SELECCIÓN DE LOS MEJORES SARIMA (Solo viendo el Train)
# -----------------------------------------------------------

top_sarima <- res_sarima_train$mejores_candidatos

df_a_latex_booktabs(top_sarima, 
                    caption = "Modelos SARIMA preseleccionados tras la fuerza bruta sobre el conjunto de entrenamiento de la serie de demanda eléctrica francesa", 
                    label = "tab:fr_sarima_preseleccionados")

# Estos con los que vamos a probar no cumplen normalidad

print(top_sarima[, c("Modelo", "AIC")])


# 3. RONDA DE PREDICCIONES (SARIMAS vs SUAVIZADO)
# -----------------------------------------------

# Creamos un dataframe para guardar los errores
tabla_torneo <- data.frame(
  Modelo = character(),
  RMSE_Validacion = numeric(),
  stringsAsFactors = FALSE
)

# --- A) EVALUAR LOS SARIMAS ---
lista_predicciones <- list() # Para guardar las líneas y graficar luego

for(i in 1:nrow(top_sarima)) {
  
  modelo_txt <- top_sarima$Modelo[i]
  
  # Extraemos números (p,d,q)(P,D,Q)
  nums <- as.numeric(unlist(regmatches(modelo_txt, gregexpr("[0-9]+", modelo_txt))))
  
  # Ajustamos a Train
  fit <- Arima(train_set, 
               order = c(nums[1], nums[2], nums[3]),
               seasonal = list(order = c(nums[4], nums[5], nums[6]), period = 12),
               lambda = lambda_train,
               include.constant = FALSE) # Normalmente FALSE con diferencias
  
  # Predecimos
  fc <- forecast(fit, h = h_valid)
  lista_predicciones[[modelo_txt]] <- fc$mean # Guardamos para el plot
  
  # Calculamos Error
  rmse_val <- accuracy(fc, valid_set)[2, "RMSE"]
  
  # Guardamos resultado
  tabla_torneo[nrow(tabla_torneo) + 1, ] <- list(modelo_txt, rmse_val)
}

# --- B) EVALUAR SUAVIZADO (HOLT-WINTERS) ---
# Usamos ets() que es más robusto y busca automáticamente la mejor forma

fit_ets <- ets(train_set)
fc_ets  <- forecast(fit_ets, h = h_valid)
lista_predicciones[["ETS/Holt-Winters"]] <- fc_ets$mean

rmse_ets <- accuracy(fc_ets, valid_set)[2, "RMSE"]
tabla_torneo[nrow(tabla_torneo) + 1, ] <- list("ETS/Holt-Winters", rmse_ets)


# 4. RESULTADOS Y GANADOR
# -----------------------
# Ordenamos de menor a mayor error
ranking_final <- tabla_torneo[order(tabla_torneo$RMSE_Validacion), ]

cat("\n--- CLASIFICACIÓN FINAL DEL TORNEO (Menor RMSE Gana) ---\n")
print(ranking_final)

# Modelo RMSE_Validacion
# 37    ETS/Holt-Winters         1072035
# 35 ARIMA(2,1,0)(0,1,1)         1133928
# 36 ARIMA(0,1,0)(0,1,1)         1158492
# 30 ARIMA(0,1,2)(0,1,1)         1159948
# 33 ARIMA(1,1,1)(0,1,1)         1168922
# 29 ARIMA(2,1,0)(0,1,2)         1281376
# 34 ARIMA(0,1,0)(0,1,2)         1307472
# 32 ARIMA(1,1,0)(0,1,2)         1308236
# 26 ARIMA(0,1,1)(0,1,2)         1322360
# 15 ARIMA(1,1,1)(0,1,2)         1386636
# 12 ARIMA(0,1,2)(0,1,2)         1389303
# 22 ARIMA(1,1,3)(0,1,2)         1391972
# 5  ARIMA(0,1,2)(1,1,1)         1537867
# 19 ARIMA(1,1,1)(1,1,0)         1552798
# 6  ARIMA(0,1,2)(2,1,0)         1583576
# 16 ARIMA(2,1,0)(1,1,1)         1590115
# 28 ARIMA(0,1,0)(1,1,1)         1605155
# 18 ARIMA(3,1,1)(2,1,1)         1631877
# 10 ARIMA(2,1,0)(2,1,0)         1632663
# 24 ARIMA(1,1,0)(1,1,1)         1632719
# 11 ARIMA(0,1,1)(1,1,1)         1633751
# 20 ARIMA(3,1,1)(3,1,0)         1642991
# 1  ARIMA(1,1,1)(1,1,1)         1645018
# 9  ARIMA(0,1,1)(2,1,0)         1675317
# 27 ARIMA(0,1,0)(2,1,0)         1676918
# 23 ARIMA(1,1,0)(2,1,0)         1688872
# 3  ARIMA(1,1,1)(2,1,0)         1691009
# 17 ARIMA(3,1,1)(1,1,1)         1693961
# 4  ARIMA(0,1,2)(0,1,3)         1694124
# 8  ARIMA(1,1,3)(0,1,3)         1695288
# 31 ARIMA(0,1,0)(1,1,0)         1779451
# 25 ARIMA(0,1,0)(0,1,3)         1807235
# 2  ARIMA(1,1,1)(0,1,3)         1819057
# 13 ARIMA(0,1,1)(0,1,3)         1836697
# 7  ARIMA(2,1,0)(0,1,3)         1855553
# 21 ARIMA(2,1,1)(2,1,1)         1895745
# 14 ARIMA(3,1,1)(0,1,3)         1936440


df_a_latex_booktabs(ranking_final, 
                    caption = "Resultados finales del torneo de validación entre los modelos SARIMA y el modelo de suavizado ETS/Holt-Winters para la serie de demanda eléctrica francesa", 
                    label = "tab:fr_torneo_final")

ganador_nombre <- ranking_final$Modelo[1]
cat("\nEL GANADOR ES:", ganador_nombre, "\n")


# 5. GRÁFICO DE LA BATALLA: ETS vs MEJOR SARIMA
# ---------------------------------------------

# 1. Identificamos dinámicamente el nombre del mejor SARIMA
# Buscamos en el ranking el primer modelo que NO sea el ETS
nombre_mejor_sarima <- ranking_final$Modelo[ranking_final$Modelo != "ETS/Holt-Winters"][1]

# 2. Configuramos el plot
par(mfrow=c(1,1))

# Pintamos la realidad (hacemos zoom en la parte final)
plot(serie_para_modelar, 
     xlim = c(time(valid_set)[1] - 2, time(valid_set)[h_valid]),
     main = paste("Validación: ETS vs Mejor SARIMA (", nombre_mejor_sarima, ")", sep=""),
     ylab = "MWh", lwd = 2, col="black")

# 3. Pintamos el Mejor SARIMA explícitamente - Azul Continuo
lines(lista_predicciones[[nombre_mejor_sarima]], col="blue", lwd=2, lty=1)

# 4. Pintamos el Suavizado (ETS) explícitamente - Rojo Discontinuo
lines(lista_predicciones[["ETS/Holt-Winters"]], col="red", lwd=2, lty=2)


# 5. Leyenda actualizada
legend("topleft", 
       legend = c("Realidad", "ETS/Holt-Winters", nombre_mejor_sarima),
       col = c("black", "red", "blue"),
       lty = c(1, 2, 1), 
       lwd = 2, cex=0.8, bg="white")
grid()




################################################################################
# 7. PREDICCIÓN FINAL A FUTURO (CON EL GANADOR RE-ENTRENADO)                   #
################################################################################

cat("\n--- RE-ENTRENANDO CON EL 100% DE LOS DATOS PARA PREDECIR EL FUTURO ---\n")

h_futuro <- 52 # 

ganador_final <- ganador_nombre

modelo_definitivo <- HoltWinters(serie_para_modelar, seasonal = "additive")
fc_final <- forecast(modelo_definitivo, h = h_futuro)

# Gráfico Final
plot(fc_final, 
     main = paste("Predicción Final a 2030 (Modelo:", ganador_final, ")"),
     ylab = "Demanda (MWh)", col="darkblue", flwd=2)
grid()

# Detalle numérico
print((fc_final$mean))








################################################################################
# 10. TORNEO DE VALIDACIÓN: TRAIN (80%) vs TEST (20%)  SIN TRANSFORMACIÓN       #
################################################################################

# 1. DIVISIÓN DE LA SERIE (80% - 20%)
# -----------------------------------
n_total <- length(serie_para_modelar)
n_train <- floor(0.80 * n_total) # Calculamos el índice del 80%

# Usamos la función subset (del paquete forecast/base) que respeta las fechas ts
train_set <- subset(serie_para_modelar, end = n_train)
valid_set <- subset(serie_para_modelar, start = n_train + 1)
h_valid   <- length(valid_set)

cat("Datos Totales:", n_total, "meses.\n")
cat("Entrenamiento (80%):", length(train_set), "meses.\n")
cat("Validación (20%):   ", h_valid, "meses.\n")

# Caluclamos el lambda en el conjunto de train
lambda_train <- BoxCox.lambda(train_set)
cat("Lambda calculado en Train:", round(lambda_train, 4), "\n")

# VAMOS A EJECUTARLO SIN LA TRANSFORMACIÓN


# 2. ENTRENAMIENTO DEL CANDIDATO 1: SARIMA (Fuerza Bruta en Train)
# ----------------------------------------------------------------
cat("\n--- BUSCANDO MEJOR SARIMA PARA EL CONJUNTO DE TRAIN ---\n")


# Ejecutamos fuerza bruta SOLO sobre el train
res_sarima_train <- fuerza_bruta_arima(train_set, 
                                       d = 1, 
                                       D = 1, 
                                       max_p = 3, max_q = 3, 
                                       max_P = 3, max_Q = 3) 

# Seleccionamos el mejor del train
# Priorizamos Normalidad, si no, AIC
candidatos_norm <- subset(res_sarima_train$mejores_candidatos, Norm_Lillie == TRUE)
print(candidatos_norm)

# TAMPOCO HAY CANDIDATOS QUE CUMPLAN NORMALIDAD

# Nos quedamos con lo que teníamos antes de transformar


