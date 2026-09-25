# Modelo Cinético del Reformado Seco de Metano (DRM) — DWSIM

Repositorio de apoyo a la tesis doctoral **"Simulación metodológica para implementar un modelo cinético detallado utilizando DWSIM para el reformado seco de metano en un reactor de metal líquido"**, desarrollada por Bruno Agún en Escuela Técnica Superior de Ingenieros Industriales de la Universidad Politécnica de Madrid (UPM), bajo la dirección de Alberto Abánades.

Este repositorio contiene la implementación completa del modelo cinético de tipo Langmuir–Hinshelwood (L-H) para las reacciones de reformado seco de metano (DRM) y de desplazamiento inverso de gas de agua (RWGS) sobre estaño líquido, simulado en el software de código abierto **DWSIM**, junto con las hojas de cálculo de verificación numérica empleadas a lo largo de la tesis.

## Contenido del repositorio

| Archivo | Descripción |
|---|---|
| `DRM_CasoBase_Tesis_BAG.dwxmz` | Archivo de simulación de DWSIM con el caso base del reactor: reactor de flujo pistón (PFR) isotérmico, geometría y condiciones de operación según Plevan et al. (2015), con las dos reacciones cinéticas (DRM y RWGS) configuradas como objetos de reacción heterogénea catalítica (*Heterogeneous Catalyst Reaction*, cinética *Advanced Kinetics*). |
| `DRM_Script` | Script de Python embebido en DWSIM que calcula la velocidad de la reacción de reformado seco ($r_\text{DRM}$), de acuerdo con la expresión cinética L-H derivada en el Capítulo 3 (Sección 3.5.) de la tesis. |
| `RWGS_Script` | Script de Python embebido en DWSIM que calcula la velocidad de la reacción de desplazamiento inverso de gas de agua ($r_\text{RWGS}$), con el mismo denominador de adsorción $\Omega$ que el script anterior derivada en el Capítulo 3 (Sección 3.5.) de la tesis. |
| `Calculos_DRM_BAG.xlsx` | Conjunto de hojas de cálculo con las verificaciones numéricas realizadas a lo largo de la tesis (modificaciones específicas de los metales líquidos; transferibilidad de los parámetros cinéticos con el cambio de dominio y el escalado BEP; el modelado, formación y gestión de carbono; la validación y consistencia termodinámica; la elección de la geometría del reactor con su justificación y el cálculo del holdup de gas). **La primera hoja contiene un índice que indica a qué sección de la tesis corresponde cada hoja del archivo.** |
| `Analsis_Carbon.py` | Archivo que calcula la formación del carbono de forma desacoplada.|

## Resumen del modelo

El reactor se modela como un PFR isotérmico que representa la zona de contacto gas–estaño de una columna de burbujeo, con la geometría del reactor experimental de Plevan et al. (2015) (diámetro interior 35,9 mm, altura de llenado de estaño 600 mm). Las velocidades de reacción se implementan como funciones de Python invocadas por DWSIM en cada paso de integración axial, a partir de la temperatura, la presión y la composición local:

$$
r_\text{DRM} = \frac{k_4\,K_1\,p_{\text{CH}_4}\left(1-\dfrac{p_{\text{CO}}^2\,p_{\text{H}_2}^2}{K_{eq,\text{DRM}}\,p_{\text{CO}_2}\,p_{\text{CH}_4}}\right)}{\Omega^2}
\qquad
r_\text{RWGS} = \frac{k_{11}\,K_2\sqrt{K_3}\,p_{\text{CO}_2}\sqrt{p_{\text{H}_2}}\left(1-\dfrac{p_{\text{CO}}\,p_{\text{H}_2\text{O}}}{K_{eq,\text{RWGS}}\,p_{\text{CO}_2}\,p_{\text{H}_2}}\right)}{\Omega^2}
$$

donde $\Omega = 1 + K_1 p_{\text{CH}_4} + K_2 p_{\text{CO}_2} + \sqrt{K_3\,p_{\text{H}_2}} + K_{-21}\,p_{\text{CO}} + K_{-23}\,p_{\text{H}_2\text{O}}$ es el denominador de adsorción compartido por ambas reacciones, que refleja la competencia de todas las especies presentes por los mismos sitios activos de la superficie del metal líquido.

Los parámetros cinéticos ($A_\text{DRM}=1{,}29\times10^6$, $E_{a,\text{DRM}}=197\,000$ J/mol; $A_\text{RWGS}=0{,}35\times10^6$, $E_{a,\text{RWGS}}=97\,000$ J/mol) y las entalpías de adsorción de las cinco especies se han obtenido transfiriendo la cinética de un catalizador de referencia de Ni/Al₂O₃ al sistema de estaño líquido mediante el escalado de la relación de Brønsted–Evans–Polanyi (BEP), sin ajuste empírico alguno a datos propios del sistema de metal líquido (para el que no existen, hasta la fecha, datos cinéticos experimentales publicados).

### Caso base

| Parámetro | Valor |
|---|---|
| Temperatura | 1173 K (900 °C) |
| Presión | 101 325 Pa (1 atm) |
| Alimentación | CH₄:CO₂ = 1:1 (equimolar) |
| Caudal volumétrico | 100 mL$_n$/min |
| Densidad del Sn líquido | 6400 kg/m³ |
| Fracción de huecos (*gas holdup*) | 2,43 % (correlación de flujo de deriva de Kataoka e Ishii) |
| Carga de catalizador (*Catalyst Loading*) | 6244.7 kg/m³ |
| Paquete termodinámico | Peng-Robinson |

Como resultado del trabajo se obtiene la existencia de una **transición de régimen**: el reactor opera bajo control cinético a baja temperatura y bajo control termodinámico a alta temperatura, con la frontera entre ambos regímenes situada en torno a **1050–1073 K**. Esta frontera constituye una predicción cuantitativa y contrastable experimentalmente del modelo, y su localización precisa es uno de los resultados originales de la tesis.

## Cómo usar este repositorio

1. Descargar e instalar [DWSIM](https://dwsim.org/) (la tesis emplea la versión 9.0.2).
2. Abrir `DRM_CasoBase_Tesis_BAG.dwxmz` directamente en DWSIM: el diagrama de flujo, las corrientes de alimentación y salida, y los dos objetos de reacción (con `DRM_Script` y `RWGS_Script` ya vinculados como cinética avanzada) están completamente configurados para el caso base.
3. Para reproducir el barrido de temperatura de la tesis (723–1273 K, intervalos de 50 K), modificar la temperatura de la corriente de alimentación y resolver de nuevo, o emplear el módulo de análisis de sensibilidad nativo de DWSIM.

## Referencia

Si este código resulta de utilidad para su trabajo, se agradece la cita de los artículos asociados:

> Agún, B., Abánades, A. (2025). *Comprehensive review on dry reforming of methane: Challenges and potential for greenhouse gas mitigation*. International Journal of Hydrogen Energy. (DOI: 10.1016/j.ijhydene.2025.01.160). 

> Agún, B., Abánades, A. (2026). *A methodological simulation framework for implementing a detailed kinetic model using DWSIM for dry reforming of methane in a liquid metal reactor*. International Journal of Hydrogen Energy. (DOI: 10.1016/j.ijhydene.2026.156383).

y/o de la tesis doctoral de la que procede.

## Autor

**Bruno Agún** — Doctorando en Energía Sostenible, Nuclear y Renovable, Universidad Politécnica de Madrid (UPM)

Director de tesis: Alberto Abánades
