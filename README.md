# Integrity Risk Income Analysis

## Descripción del Proyecto
Este repositorio alberga un estudio descriptivo y un análisis transversal del **Índice IRIS de SCImago** (Integrity Risk Index), evaluado a través de diferentes grupos de ingreso y países. El proyecto tiene como principal objetivo explorar y comprender cuantitativamente cómo se distribuye y varía el riesgo de integridad científica a nivel global, dependiendo del nivel de ingresos y características socioeconómicas de cada nación.

Los datos originales y métricas provienen de la plataforma oficial: [SCImago IRIS](https://www.scimagoiris.com/).

## Estado Actual y Metodología
Hasta el momento, el proyecto ha completado su fase inicial de recolección y estructuración de los datos:
- **Extracción de datos (Web Scraping / API):** Se desarrolló código en Python (utilizando el entorno de Google Colab) para consultar, extraer y limpiar los datos del índice directamente desde la página web de SCImago IRIS.
- **Almacenamiento y Estructuración:** Los datos brutos extraídos han sido consolidados y guardados localmente en formato CSV (`Scimago_IRIS_Index_Data.csv`). Esta base de datos estructurada será el principal insumo para los futuros análisis estadísticos, modelados y visualizaciones transversales.

## Contenido del Repositorio
- `integrity_risk_extraction_data.py` / `Integrity_Risk_Extraction_Data.ipynb`: Scripts de Python y Jupyter Notebooks utilizados para la extracción inicial y procesamiento de información desde la web.
- `Scimago_IRIS_Index_Data.csv`: Conjunto de datos final resultante y procesado tras la extracción.
- `README.md`: Este archivo, que documenta el propósito, evolución y estructura general del proyecto.
