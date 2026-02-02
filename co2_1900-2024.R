################################################################################
#                             TFM - SERIE CO_2                                 #
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
emmisions <- read.csv("annual-co2-emissions-per-country.filtered/annual-co2-emissions-per-country.csv")
emmisions <- emmisions[, 3:4] # Asumiendo que col 3 es Año y col 4 es CO2

# Vamos a eliminar las emisiones anteriores al 1900
emmisions <- emmisions[emmisions[, 1] >= 1900, ]

# Conversión matemática a Mt (Megatoneladas)
emmisions[, 2] <- emmisions[, 2] / 1000000

# Creamos la Serie Temporal
emis_ts <- ts(emmisions[, 2], 
              start = c(min(emmisions[, 1]), 1), 
              frequency = 1)

# Graficamos
plot(emis_ts, 
     main = "Emisiones de CO2 (Serie Temporal)",
     xlab = "Año",
     ylab = "Millones de Toneladas (Mt)",
     col = "blue",
     lwd = 2)


################################################################################
# 2.1. ESTUDIO DE OUTLIERS E INTERVENCIONES                                    #
################################################################################

# Queremos estudiar los outliers tanto aditivos (AO) como cambios de nivel (LS)
# Empecemos con los outliers aditivos
outliers.co2.ao <- tso(emis_ts, types = c("AO"), maxit.iloop = 10)

# Veamos los resultados
print(outliers.co2.ao)
plot(outliers.co2.ao)

# RESULTADO:
# --------------------------------------------------------------------------
# Outliers:
#   type ind time coefhat   tstat
# 1   AO 110 2009   -1185  -5.294  <-- Crisis Financiera 2008
# 2   AO 271 2020   -1817  -8.115  <-- COVID-19

# Es por tanto que no son outliers, son intervenciones (pues conocemos su causa y se debe a efectos externos)

# --------------------------------------------------------------------------

# Ahora veamos los cambios de nivel (LS)
outliers.co2.ls <- tso(emis_ts, types = c("LS"), maxit.iloop = 10)

# Vemos los resultados
print(outliers.co2.ls)
plot(outliers.co2.ls)

# RESULTADO DEL ANÁLISIS DE OUTLIERS (Serie 1900-2024):
# --------------------------------------------------------------------------
# Outliers detectados (tso):
#   type  ind  time   coefhat   tstat   <-- Interpretación Histórica
# 1   LS   46  1945    4746     3.696   <-- Fin WWII: Inicio de la "Gran Aceleración" y reconstrucción industrial.
# 2   LS   82  1981    6353     5.800   <-- Recuperación Post-Crisis Petróleo '79 y auge de nuevos mercados.
# 3   LS   93  1992    5313     5.004   <-- Globalización Post-Guerra Fría y expansión de economías emergentes.
# 4   LS  104  2003    5448     4.993   <-- Boom de China: Entrada en la OMC e industrialización masiva (carbón).
# 5   LS  111  2010    4278     3.472   <-- Rebote Post-Crisis: Estímulos fiscales y recuperación acelerada.

# --------------------------------------------------------------------------

################################################################################
# 2.2. GESTIÓN DE INTERVENCIONES (PREPARACIÓN PARA MODELAR)                    #
################################################################################

# ESTRATEGIA: 
# 1. Los AO (Guerras, COVID) son ruido temporal -> LOS QUITAMOS (Linealizamos).
#    No queremos que nuestro modelo aprenda que cada cierto tiempo puede haber
#    caídas bruscas en las emisiones, ya que se deben a eventos externos.

# 2. Los LS son cambios de nivel -> LO DEJAMOS (No se corrige).
#    Pues las predicciones futuras han de partir de ese escalón más alto.

# Hacemos una copia de la serie original
serie_para_modelar <- emis_ts

# Sacamos la tabla de intervenciones detectadas
tabla_outliers <- outliers.co2.ao$outliers

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
ts.plot(emis_ts, serie_para_modelar,
        col = c("red", "blue"),
        lwd = c(1, 2),
        main = "Comparativa: Realidad vs. Serie para Modelar",
        ylab = "Mt CO2")

legend("topleft", 
       legend=c("Serie Original (Con Intervenciones)", "Serie Limpia de AO (Mantiene LS)"),
       col=c("red", "blue"), lty=1, lwd=c(1,2))

# RESULTADO :
# La línea azul es idéntica a la roja, EXCEPTO en las fechas anteriormente señaladas.
# En 2020 se ve que la línea azul NO cae (ignora el COVID), sigue recta.
# Sin embargo, a partir de 2010, se ve que ambas mantienen el "salto" hacia arriba.
# Esta serie 'serie_para_modelar' es la que se usará en el siguiente paso.


################################################################################
# 3. TRANSFORMACIONES DE LA SERIE                                              #
################################################################################

# NOTA SOBRE DECOMPOSE:
# Al ser datos anuales (frequency=1), NO se puede usar decompose() ni stl(),
# ya que no existe componente estacional que separar.
# Asumimos que la serie es TENDENCIA + RUIDO.

# --- TRANSFORMACIÓN BOX-COX ---

# 1. Calculamos el Lambda óptimo
# Lambda nos dice qué transformación necesita la serie:
# lambda = 1  -> No hacer nada (los datos están bien)
# lambda = 0  -> Aplicar Logaritmo (log)
# lambda = 0.5 -> Raíz cuadrada
lambda_optimo <- BoxCox.lambda(serie_para_modelar)

print(paste("El valor de Lambda óptimo es:", round(lambda_optimo, 4)))

# 2. Aplicamos la transformación
serie_boxcox <- BoxCox(serie_para_modelar, lambda_optimo)

# 3. Gráfico Comparativo: Original vs Transformada
par(mfrow = c(1, 1))

par(mfrow = c(1, 2)) # Dividimos pantalla en 2 columnas

# Gráfico 1: Original
plot(serie_para_modelar, 
     main = "Serie Original (Limpia de AO)",
     ylab = "Mt CO2", col = "black")
grid()

# Gráfico 2: Transformada
plot(serie_boxcox, 
     main = paste("Transformada Box-Cox (lambda=", round(lambda_optimo, 2), ")", sep=""),
     ylab = "Valor Transformado", col = "blue")
grid()

# Volvemos a pantalla normal
par(mfrow = c(1, 1))

################################################################################
# 3. CALCULAMOS LAS DIFERENCIAS, ADF, FAS y FAP                                #
################################################################################

# NOTA IMPORTANTE:
# Al ser datos ANUALES, no existe estacionalidad (lag=12).
# Saltamos nsdiffs y pasamos directamente a diferencias regulares (tendencia).

# 1. Calcular cuántas diferencias regulares necesitamos para quitar la tendencia
# Usamos la serie limpia que creamos antes ('serie_para_modelar' o 'serie_boxcox')
numero.diffreg <- ndiffs(serie_boxcox)

cat('El número de diferencias regulares necesarias es:', numero.diffreg)

# 2. Aplicamos las diferencias
# lag = 1 porque son datos anuales
serie_estacionaria = diff(serie_boxcox, lag = 1, differences = numero.diffreg)

# 3. Graficamos la serie diferenciada
# Deberías ver que ya no crece hacia arriba, sino que oscila alrededor del cero (media constante)
plot(serie_estacionaria, 
     main = "Serie Estacionaria (Diferenciada)",
     ylab = "Variación de Emisiones")
grid()
abline(h=0, col="red", lty=2) # Línea en el cero para referencia


# 4. Comprobación matemática: Test de Dickey-Fuller Aumentado (ADF)
# H0 (Nula): La serie NO es estacionaria
# H1 (Alternativa): La serie SÍ es estacionaria
# Buscamos un p-value < 0.05
test_adf <- adf.test(serie_estacionaria)
print(test_adf)



# Y si solo hacemos 1 diferencia?
serie_estacionaria_1diff = diff(serie_boxcox, lag = 1, differences = 1)
test_adf_1diff <- adf.test(serie_estacionaria_1diff)
print(test_adf_1diff)

# Tanto con 1 diferencia como con dos nos queda que la serie es estacionaria

# 5. Correlogramas (ACF y PACF)
# Podemos intentar decidir un modelo ARIMA (p, d, q)
par(mfrow=c(1,2))
acf(serie_estacionaria, lag.max = 70, main="ACF (MA - q)")
pacf(serie_estacionaria,lag.max = 70, main="PACF (AR - p)")
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
# 6. SELECCIÓN DEL MODELO ENTRE LOS VÁLIDOS                                    #
################################################################################

# 1. Ejecutamos la función 
resultados_co2 <- fuerza_bruta_arima(serie_boxcox, 
                                     d = 2,       # Ya calculamos con ndiffs que era 2
                                     max_p = 4,   # Probamos AR hasta orden 3
                                     max_q = 4)   # Probamos MA hasta orden 3

resultados_co2_d1 <- fuerza_bruta_arima(serie_boxcox, 
                                     d = 1,       # Ya calculamos con ndiffs que era 2
                                     max_p = 4,   # Probamos AR hasta orden 3
                                     max_q = 4)   # Probamos MA hasta orden 3


# 2. Vemos los "finalistas"
# (Estos son los que cumplen Coeficientes Significativos + Independencia + Homocedasticidad)
print("--- MEJORES CANDIDATOS (Ordenados por AIC) ---")
print(resultados_co2$mejores_candidatos)
print(resultados_co2_d1$mejores_candidatos)

# > print(resultados_co2$mejores_candidatos)
#      ID       Modelo      AIC Signif_Coef Indep_LB Norm_Lillie Homo_BP Pval_LB Lag_Usado
# BP16 17 ARIMA(2,2,3) 870.4402        TRUE     TRUE        TRUE    TRUE   0.218        10
# BP4   5 ARIMA(0,2,1) 873.7415        TRUE     TRUE       FALSE    TRUE   0.142        10

# > print(resultados_co2_d1$mejores_candidatos)
#      ID       Modelo      AIC Signif_Coef Indep_LB Norm_Lillie Homo_BP Pval_LB Lag_Usado
# BP14 18 ARIMA(3,1,3) 878.1933        TRUE     TRUE       FALSE    TRUE   0.124        10
# BP5   6 ARIMA(1,1,1) 881.5236        TRUE     TRUE       FALSE    TRUE   0.091        10



# Solo nos queda un modelo que verifique todas las hipótesis estadísticas:
# ARIMA(2,2,3)


################################################################################
# 7. PREDICCIONES                                                              #
################################################################################

# Seleccionamos el mejor modelo (el primero de la lista de mejores candidatos)
# En este caso, el de menor AIC es el que cumple normalidad en los residuos
# condición que no habíamos marcado, pues solamente nos sirve para Intervalos
# de predicción.

mejor_modelo <- resultados_co2$mejores_candidatos$Modelo[1]
cat("El mejor modelo seleccionado es:", mejor_modelo, "\n")

# 1. Recuperamos el lambda que calculamos al principio

# 2. Ajustamos el modelo USANDO LA SERIE ORIGINAL 
# Le pasamos el lambda aquí dentro. R hará la transformación internamente.
modelo_final_auto <- Arima(serie_para_modelar, 
                           order = c(2, 2, 3),
                           lambda = lambda_optimo,
                           include.constant = TRUE)

# 3. Predicción (automáticamente en Mt CO2)
prediccion_real <- forecast(modelo_final_auto, h = 26) # Hasta 2050

# 4. Gráfico final
plot(prediccion_real, 
     main = "Predicción de Emisiones CO2 (2024-2040)",
     ylab = "Millones de Toneladas (Mt CO2)",
     xlab = "Año",
     flwd = 2) # Grosor de la línea de predicción
grid()


################################################################################
# 8. SUAVIZADO.                                                                #
################################################################################

# Tenemos que ver que tipo de suavizado tiene la serie, aplicarlo, predecir y 
# comparar con la modelización ARIMA

# En primer lugar, el suavizado exponencial simple

s_exp_simple <- ses(serie_para_modelar, h = 24)
plot(s_exp_simple)
summary(s_exp_simple)

# Podemos ver que el suavizado exponencial simple no se ajusta bien a los datos.


# En segundo lugar, el suavizado exponencial doble

s_exp_doble <- ets(serie_para_modelar)
pred_sexpdoble <- forecast(s_exp_doble, h=24)
plot(pred_sexpdoble)
summary(s_exp_doble)







################################################################################
# CON SOLO UNA DIFERENCIA FALLA NORMALIDAD EN LOS MODELOS                      #
################################################################################

# Aunque el test de DF nos dice que se verifica la hipótesis de estacionariedad, 
# pues el pvalor es inferior a 0.01, en los modelos resultantes no se verifica
# la hipótesis de normalidad de los residuos.


# 1. Ajustamos el mejor modelo con d=1 (aunque falle normalidad)
#    Fuerza bruta dice que el mejor d=1 es ARIMA(3,1,3)
modelo_d1 <- Arima(serie_para_modelar, 
                   order = c(3, 1, 3), 
                   lambda = lambda_optimo, include.constant = TRUE)

# 2. Ajustamos el modelo con d=2 (el que te da normalidad)
modelo_d2 <- Arima(serie_para_modelar, 
                   order = c(2, 2, 3), # El que tenías puesto
                   lambda = lambda_optimo, include.constant = TRUE)

# 3. Predecimos a muy largo plazo (2050 son 26 años)
pred_d1 <- forecast(modelo_d1, h = 26)
pred_d2 <- forecast(modelo_d2, h = 26)

# 4. Graficamos juntos
par(mfrow=c(1,1))
plot(pred_d2, main = "Peligro: d=1 vs d=2", col="red", fcol="pink", ylim=c(0, 100000))
lines(pred_d1$mean, col="blue", lwd=2)
legend("topleft", legend=c("d=2 (Explosivo)", "d=1 (Tendencial)"), 
       col=c("red", "blue"), lty=1, lwd=2)








################################################################################
# 9. SELECCIÓN ROBUSTA: FUERZA BRUTA + VALIDACIÓN                              #
################################################################################

# 1. DIVISIÓN DE DATOS (ORIGINALES LIMPIOS DE OUTLIERS)
# -----------------------------------------------------
# Usamos 'serie_para_modelar' (que ya tiene corregidos los AO de 2009 y 2020)
train_set <- window(serie_para_modelar, end = 2019)
valid_set <- window(serie_para_modelar, start = 2020)
h_valid <- length(valid_set)

cat("Entrenamiento: 1900-2019 (", length(train_set), "años)\n")
cat("Validación:    2020-2024 (", h_valid, "años)\n")


# 2. CÁLCULO DE LAMBDA (SOLO CON TRAIN) 
# -------------------------------------------------------------
# El algoritmo solo ve el pasado. No sabe cómo se comporta la varianza en 2024.
lambda_train <- BoxCox.lambda(train_set)

print(paste("Lambda calculado solo con Train (1900-2019):", round(lambda_train, 4)))

# Aplicamos este lambda al train para poder hacer la fuerza bruta
train_boxcox <- BoxCox(train_set, lambda = lambda_train)


# 3. EJECUTAMOS FUERZA BRUTA (SOBRE TRAIN TRANSFORMADO)
# -----------------------------------------------------
# Buscamos patrones lineales en los datos históricos transformados

cat("\n--- Buscando candidatos con d=1 ---\n")
# Pasamos train_boxcox que ya está transformado con el lambda purista
res_d1 <- fuerza_bruta_arima(train_boxcox, d=1, max_p=3, max_q=3)

cat("\n--- Buscando candidatos con d=2 ---\n")
res_d2 <- fuerza_bruta_arima(train_boxcox, d=2, max_p=3, max_q=3)

# Unimos los mejores
candidatos_d1 <- head(res_d1$mejores_candidatos, 3) 
candidatos_d2 <- head(res_d2$mejores_candidatos, 3) 
todos_candidatos <- rbind(candidatos_d1, candidatos_d2)

print("--- FINALISTAS (Basado en AIC del pasado) ---")
print(todos_candidatos[, ])


# 4. EL TORNEO: VALIDACIÓN (Usando lambda_train)
# ----------------------------------------------
resultados_validacion <- data.frame()

for(i in 1:nrow(todos_candidatos)) {
  
  modelo_texto <- todos_candidatos$Modelo[i]
  nums <- as.numeric(unlist(regmatches(modelo_texto, gregexpr("[0-9]+", modelo_texto))))
  
  # AJUSTAMOS EL MODELO
  # Usamos train_set (original) + lambda_train
  # R transformará el train con ese lambda, ajustará, y al predecir destransformará.
  fit <- Arima(train_set, 
               order = c(nums[1], nums[2], nums[3]),
               lambda = lambda_train,   # <--- CLAVE: Usamos el lambda del train
               include.constant = TRUE)
  
  # PREDECIMOS 
  # Al comparar con valid_set, vemos qué tal funciona ese lambda + esos coeficientes
  fc <- forecast(fit, h = h_valid)
  
  acc <- accuracy(fc, valid_set)
  rmse_val <- acc[2, "RMSE"] 
  
  resultados_validacion <- rbind(resultados_validacion, data.frame(
    Modelo = modelo_texto,
    AIC_Train = fit$aic,
    RMSE_Validacion = rmse_val
  ))
}

# 5. RESULTADO FINAL
ranking_final <- resultados_validacion[order(resultados_validacion$RMSE_Validacion), ]
print(ranking_final)

ganador_texto <- ranking_final$Modelo[1]

# Los modelos que salen son los siguientes:

#         Modelo AIC_Train RMSE_Validacion
# 1 ARIMA(3,1,3)  918.3458        737.9299
# 2 ARIMA(0,2,1)  913.9786       1099.1218


cat("\nEl modelo ganador es:", ganador_texto, "\n")


# --------------------------------------------------------
# AÑADIMOS EL CANDIDATO DE SUAVIZADO (ETS) AL TORNEO
# --------------------------------------------------------

# 1. Entrenamos ETS con el train_set y el lambda_train
modelo_ets_train <- ets(train_set)

# 2. Predecimos valid_set
fc_ets <- forecast(modelo_ets_train, h = h_valid)

# 3. Calculamos RMSE
acc_ets <- accuracy(fc_ets, valid_set)
rmse_ets <- acc_ets[2, "RMSE"]

# 4. Lo añadimos a la tabla de resultados
resultados_validacion <- rbind(resultados_validacion, data.frame(
  Modelo = "ETS (Suavizado)",
  AIC_Train = modelo_ets_train$aic,
  RMSE_Validacion = rmse_ets
))


# --------------------------------------------------------
# AHORA SÍ: EL RANKING FINAL COMPLETO
# --------------------------------------------------------
ranking_final <- resultados_validacion[order(resultados_validacion$RMSE_Validacion), ]
print("--- CLASIFICACIÓN FINAL DEL TORNEO (Menor RMSE Gana) ---")
print(ranking_final)

ganador_texto <- ranking_final$Modelo[1]
cat("\n El modelo GANADOR ABSOLUTO es:", ganador_texto, "\n")

################################################################################
# VISUALIZACIÓN GRÁFICA DEL TORNEO (VALIDACIÓN)
################################################################################

# 1. Preparamos el lienzo con los datos REALES de validación (Negro grueso)
# -------------------------------------------------------------------------
par(mfrow=c(1,1))

# Calculamos límites para que se vean bien todas las líneas
y_lims <- range(window(serie_para_modelar, start = 2017)) * c(0.96, 1.02) 

plot(window(serie_para_modelar, start = 2017), 
     main = "Torneo de Modelos: Predicciones en el Conjunto de Validación (2020-2024)",
     ylab = "Mt CO2", xlab = "Año",
     ylim = y_lims,
     lwd = 3, col = "black", type = "o", pch = 20) # Realidad en negro con puntos

grid()

# 2. Bucle para graficar cada modelo del ranking
# ----------------------------------------------
colores <- c("blue", "red", "green4", "purple", "orange") # Paleta de colores
ltys <- c(1, 2, 4, 5, 3) # Tipos de línea para diferenciar

# Vectores para la leyenda final
leyenda_txt <- c("Datos Reales")
leyenda_col <- c("black")
leyenda_lty <- c(1)
leyenda_lwd <- c(3)

# Recorremos el ranking (limitado a los top 5 si hubiera muchos)
n_modelos <- min(nrow(ranking_final), 5) 

for(i in 1:n_modelos) {
  
  modelo_nombre <- ranking_final$Modelo[i]
  rmse_val <- round(ranking_final$RMSE_Validacion[i], 2)
  
  # A) SI ES ETS
  if(modelo_nombre == "ETS (Suavizado)") {
    fit <- ets(train_set) # Ajuste sin lambda externo (ets lo gestiona o usas el tuyo)
    fc <- forecast(fit, h = h_valid)
    
    # B) SI ES ARIMA
  } else {
    # Extraemos los órdenes (p,d,q) del texto "ARIMA(3,1,3)"
    nums <- as.numeric(unlist(regmatches(modelo_nombre, gregexpr("[0-9]+", modelo_nombre))))
    
    fit <- Arima(train_set, 
                 order = c(nums[1], nums[2], nums[3]),
                 lambda = lambda_train, 
                 include.constant = TRUE)
    fc <- forecast(fit, h = h_valid)
  }
  
  # Pintamos la línea del modelo
  lines(fc$mean, col = colores[i], lwd = 2, lty = ltys[i])
  
  # Añadimos info a la leyenda
  leyenda_txt <- c(leyenda_txt, paste0(modelo_nombre, " (RMSE: ", rmse_val, ")"))
  leyenda_col <- c(leyenda_col, colores[i])
  leyenda_lty <- c(leyenda_lty, ltys[i])
  leyenda_lwd <- c(leyenda_lwd, 2)
}

# 3. Ponemos la leyenda
legend("bottomleft", 
       legend = leyenda_txt,
       col = leyenda_col,
       lwd = leyenda_lwd,
       lty = leyenda_lty,
       bg = "white", cex = 0.8)



################################################################################
# 10. PREDICCIÓN FINAL (MODELO GANADOR + TODOS LOS DATOS)                      #
################################################################################

cat("--- GENERANDO PREDICCIÓN FINAL A 2040 CON EL GANADOR ---\n")

# 1. PREPARACIÓN DE PARÁMETROS FINALES
# ------------------------------------
# Recalculamos el Lambda usando TODA la serie (1900-2024) para reajustar los coeficientes
lambda_final <- BoxCox.lambda(serie_para_modelar)

# Definimos el horizonte
h_futuro <- 16 


# 2. RE-AJUSTE DEL MODELO SEGÚN QUIÉN HAYA GANADO
# -----------------------------------------------

if (ganador_texto == "ETS (Suavizado)") {
  
  # Si ganó ETS, reajustamos un ETS con toda la serie
  modelo_definitivo <- ets(serie_para_modelar, lambda = lambda_final)
  
} else {
  
  # Si ganó un ARIMA (ej: ARIMA(3,1,3)), extraemos sus números y reajustamos
  nums_win <- as.numeric(unlist(regmatches(ganador_texto, gregexpr("[0-9]+", ganador_texto))))
  
  modelo_definitivo <- Arima(serie_para_modelar, 
                             order = c(nums_win[1], nums_win[2], nums_win[3]),
                             lambda = lambda_final, 
                             include.constant = TRUE)
}


# 3. PREDICCIÓN Y GRÁFICO FINAL TFM
# ---------------------------------
forecast_final <- forecast(modelo_definitivo, h = h_futuro)
print(forecast_final)

# Graficamos bonito
plot(forecast_final, 
     main = paste("Proyección de Emisiones CO2 hasta 2040\nModelo:", ganador_texto),
     ylab = "Millones de Toneladas (Mt CO2)", 
     xlab = "Año", 
     flwd = 2, col="darkblue", shadecols=c("gray80", "gray70"))
grid()

# Zoom a los últimos años para ver el detalle
plot(forecast_final, xlim=c(2000, 2040),
     main = "Zoom: Detalle 2000-2040",
     ylab = "Mt CO2", col="darkblue", flwd=2)
grid()

# 4. DATOS NUMÉRICOS (Para poner en el texto del TFM)
print(head(forecast_final$mean)) # Muestra los valores predichos para los últimos años






################################################################################
# 11. LA DECISIÓN FINAL: COMPARATIVA VISUAL A 2040                             #
################################################################################

# 1. Ajustamos el "Ganador Práctico" (3,1,3) con TODOS los datos
modelo_practico <- Arima(serie_para_modelar, 
                         order = c(3, 1, 3), 
                         lambda = lambda_final, 
                         include.constant = TRUE)

# 2. Ajustamos el "Ganador Teórico" (2,2,3) con TODOS los datos
modelo_teorico <- Arima(serie_para_modelar, 
                        order = c(2, 2, 3), 
                        lambda = lambda_final, 
                        include.constant = TRUE)

# 3. Proyectamos ambos a 2040
fc_practico <- forecast(modelo_practico, h = 16)
fc_teorico  <- forecast(modelo_teorico, h = 16)

# 4. GRÁFICO COMPARATIVO
par(mfrow=c(1,1))
plot(serie_para_modelar, xlim=c(2000, 2040), ylim=c(min(serie_para_modelar), max(fc_teorico$upper, fc_practico$upper)),
     main = "Duelo Final a 2040: Robustez (d=1) vs Ajuste (d=2)",
     ylab = "Mt CO2", lwd=2)

# Línea Azul: El que ganó la validación
lines(fc_practico$mean, col="blue", lwd=3)
lines(fc_practico$lower[,2], col="blue", lty=2) # Intervalo confianza
lines(fc_practico$upper[,2], col="blue", lty=2)

# Línea Roja: El que cumple hipótesis pero tiene d=2
lines(fc_teorico$mean, col="red", lwd=3)
lines(fc_teorico$lower[,2], col="red", lty=2)
lines(fc_teorico$upper[,2], col="red", lty=2)

legend("topleft", 
       legend = c("Histórico", "ARIMA(3,1,3) - Ganador Validación", "ARIMA(2,2,3) - Cumple Hipótesis"),
       col = c("black", "blue", "red"), 
       lty = c(1, 1, 1), lwd = c(2, 3, 3))
grid()


print(tail(fc_teorico$mean))
print(tail(fc_practico$mean))




################################################################################
# CONCLUSIÓN Y SELECCIÓN FINAL DEL MODELO (CO2)                                #
################################################################################

# JUSTIFICACIÓN DE LA DECISIÓN:
# -----------------------------
# Se exploró la posibilidad de utilizar un modelo ARIMA(2,2,3) (d=2), el cual 
# presentaba un excelente ajuste estadístico sobre la serie completa (1900-2024), 
# cumpliendo todas las hipótesis de los residuos.
#
# Sin embargo, al auditar este modelo mediante validación cruzada temporal, se 
# detectó una inestabilidad estructural grave: los coeficientes del modelo no 
# resultaron significativos en el conjunto de entrenamiento (1900-2019). 
# Esto indica que la validez del modelo depende exclusivamente de los últimos 
# datos de la muestra, evidenciando un claro sobreajuste.
#
# En consecuencia, se descarta el modelo con d=2 y se selecciona definitivamente 
# el ARIMA(3,1,3) (d=1), el cual demostró ser estadísticamente robusto en el 
# entrenamiento y minimizó el error de predicción (RMSE) en la fase de validación.

# ------------------------------------------------------------------------------

# Asignación del modelo ganador para la predicción final
ganador_texto <- "ARIMA(3,1,3)"
