# energy-co2-timeseries

Este repositorio contiene el código desarrollado en el Trabajo Fin de Máster (TFM) para el
análisis y modelización de series temporales de **demanda eléctrica** y **emisiones / concentración
de CO₂**, con especial atención al contexto de la transición energética europea.

El trabajo combina técnicas de predicción univariante y análisis multivariante para estudiar
el comportamiento individual de las series y sus relaciones de largo plazo.

---

## Contenido del repositorio

- **preprocessing_load.ipynb**  
  Notebook en Python para el preprocesamiento de los datos de demanda eléctrica (ENTSO-E) y
  la construcción de la serie mensual agregada por país.

- **electricidad_es.R**  
  Modelización y predicción univariante de la demanda eléctrica mensual en España.

- **electricidad_fr.R**  
  Modelización y predicción univariante de la demanda eléctrica mensual en Francia.

- **electricidad_it.R**  
  Modelización y predicción univariante de la demanda eléctrica mensual en Italia.

- **co2_1900-2024.R**  
  Modelización y predicción univariante de la serie histórica anual de emisiones globales de CO₂.

- **coint_co2_es.R**  
  Análisis multivariante de largo plazo entre la demanda eléctrica española y la concentración
  atmosférica de CO₂.

- **coint_co2_fr.R**  
  Análisis multivariante de largo plazo entre la demanda eléctrica francesa y la concentración
  atmosférica de CO₂.

- **coint_co2_it.R**  
  Análisis multivariante de largo plazo entre la demanda eléctrica italiana y la concentración
  atmosférica de CO₂.

---

## Metodología

El flujo metodológico del trabajo se estructura en cuatro bloques:

1. **Preprocesamiento de datos**  
   Limpieza, homogeneización temporal y agregación mensual de los registros horarios de demanda
   eléctrica.

2. **Análisis univariante**  
   Modelización de cada serie mediante modelos ARIMA/SARIMA y métodos de suavizado exponencial,
   incorporando tratamiento de outliers, cambios de nivel y validación cruzada.

3. **Predicción**  
   Generación de predicciones a medio plazo y evaluación de la capacidad de generalización de
   los modelos mediante métricas de error fuera de muestra.

4. **Análisis multivariante**  
   Estudio de relaciones de largo plazo entre demanda eléctrica y concentración atmosférica de
   CO₂ mediante técnicas de cointegración y modelos de corrección del error (ECM).

---


## Datos

El código emplea datos procedentes de fuentes públicas:

- **Demanda eléctrica**: ENTSO-E (registros horarios agregados a frecuencia mensual).
- **Concentración atmosférica de CO₂**: NOAA (serie mensual global).
- **Emisiones globales de CO₂**: IEA (base histórica anual).

Por motivos de tamaño y licencia, los ficheros originales no se incluyen en el repositorio.
El código asume la disponibilidad de los siguientes archivos procesados:

- `datos_electricidad/datos_electricidad_mensual.csv`
- `co2_mm_gl.csv`
- `annual-co2-emissions-per-country.csv`

---

## Reproducibilidad

1. Ejecutar el notebook `preprocessing_load.ipynb` para generar la serie mensual de demanda eléctrica.
2. Ejecutar los scripts `electricidad_*.R` para el análisis y predicción univariante.
3. Ejecutar `co2_1900-2024.R` para la predicción de CO₂ global.
4. Ejecutar los scripts `coint_co2_*.R` para el análisis multivariante.

Se recomienda trabajar con rutas relativas y ejecutar los scripts desde el directorio raíz
del proyecto.

---

## Autor

Miguel Bande Rodríguez


---

## English summary

This repository contains the code developed for a Master's Thesis focused on the
time series analysis and modelling of **electricity demand** and **CO₂ emissions /
atmospheric concentration**, with particular attention to the European energy
transition context.

The project combines univariate forecasting techniques with multivariate time
series analysis to study both the individual behaviour of the series and their
long-run relationships.

### Main components

- Python notebook for preprocessing and aggregating electricity demand data.
- R scripts for univariate modelling and forecasting of electricity demand
  (Spain, France and Italy).
- R script for univariate modelling and forecasting of global CO₂ emissions.
- R scripts for multivariate analysis, including cointegration and error
  correction models (ECM), between electricity demand and atmospheric CO₂
  concentration.

### Data sources

- **Electricity demand**: ENTSO-E (hourly data aggregated to monthly frequency).
- **Atmospheric CO₂ concentration**: NOAA (global monthly series).
- **Global CO₂ emissions**: IEA (annual historical data).

The repository is intended for academic and research purposes.


---

## Licencia

Este repositorio se distribuye con fines académicos y de investigación.
