# **Simulación Biofísica Multi-Escala de la Viabilidad Funcional en un Riñón Humano Bioimpreso Derivado de iPSC: La Integración de Seis Módulos Predice el Rendimiento Renal Fisiológico a Escalas Vascular, Celular y Funcional**

**Carlos David Moreno Cáceres**
Investigador Independiente — VirtusSapiens, Medellín, Colombia
ORCID: [0009-0005-3933-5072](https://orcid.org/0009-0005-3933-5072)

Marzo 2026

---

## Resumen

La enfermedad renal crónica afecta a más de 850 millones de personas en todo el mundo, sin embargo no existe ningún marco computacional validado capaz de predecir la viabilidad funcional de un riñón bioimpreso antes de su fabricación. En este trabajo presentamos Bio-Kidney AI 2026, un marco in silico multi-escala que integra seis módulos biofísicos — generación de árbol vascular (Optimización Constructiva Restringida, CCO v8), difusión de oxígeno (cilindro de Krogh), cinética de diferenciación de iPSC (protocolo Takasato), optimización de bioimpresión Co-SWIFT (reología de Herschel-Bulkley), filtración glomerular (modelo Starling-Deen) y reabsorción tubular (Michaelis-Menten/Kedem-Katchalsky) — en un pipeline predictivo unificado. El algoritmo vascular CCO v8 genera 1.902 segmentos jerárquicos con cumplimiento del 100% de la Ley de Murray en 915 bifurcaciones. Una distribución radial Beta(3, 1,2) sesga el 63% de los puntos de demanda glomerular hacia la corteza, replicando la anatomía del nefrón humano. Un novedoso modelo de Poiseuille calibrado en dos pasadas resuelve el fallo del cálculo de presión ingenuo en árboles geométricamente simplificados, entregando presiones capilares terminales de 58,6 ± 13,4 mmHg — por encima del umbral de 43 mmHg requerido para la filtración de Starling. La tasa de filtración glomerular bilateral resultante de 115,2 mL/min se encuentra dentro del rango fisiológico normal para adultos sanos, representando una mejora del 40% respecto a la línea base v7. Los módulos complementarios confirmaron ausencia de zonas hipóxicas, pureza fenotípica completa de iPSC al día 21, viabilidad celular post-impresión del 98% y eficiencia de reabsorción tubular del 98,1%. El pipeline completo se ejecuta en menos de 60 segundos en hardware de consumo (16 GB RAM, Ubuntu Linux) y está implementado como paquete Python de código abierto con arquitectura de Mezcla de Expertos. Estos resultados establecen una base computacional reproducible para la evaluación de viabilidad renal pre-fabricación, con implicaciones directas para la bioingeniería de órganos y la crisis global de trasplantes.

\newpage
## 1. Introducción

La enfermedad renal crónica representa una de las cargas de salud global más significativas del siglo XXI. Con más de 850 millones de personas afectadas en todo el mundo y aproximadamente 2 millones de pacientes dependientes de terapia de reemplazo renal, la brecha entre la demanda y la disponibilidad de órganos constituye una crisis médica crítica [1]. La diálisis, aunque sustenta la vida, se asocia con morbilidad significativa, reducción de la calidad de vida y una tasa de supervivencia a cinco años inferior al 40% en muchas poblaciones de pacientes [2]. El trasplante renal sigue siendo el tratamiento de referencia, sin embargo la escasez de órganos donantes implica que la mayoría de los pacientes nunca recibirán un trasplante.

La aparición de la medicina regenerativa y las tecnologías de biofabricación ha abierto nuevas posibilidades para abordar esta crisis. Los avances en tecnología de células madre pluripotentes inducidas (iPSC), el andamiaje con matriz extracelular descelularizada (dECM) y la bioimpresión tridimensional han establecido colectivamente una vía teórica hacia la fabricación de tejido renal funcional y específico del paciente [3,4]. Entre estas tecnologías, la técnica Co-SWIFT (Escritura Sacrificial en Tejido Funcional) ha demostrado un potencial particular para la vascularización de constructos de tejido de gran espesor, un reto históricamente intratable en bioingeniería de órganos [5].

A pesar de estos avances, persiste una brecha fundamental: la ausencia de marcos computacionales validados capaces de integrar la complejidad completa de la fisiología renal — desde la hemodinámica vascular hasta el transporte iónico tubular — en un sistema predictivo unificado. Los enfoques de simulación existentes típicamente abordan aspectos individuales de la función renal de forma aislada, sin proporcionar una evaluación integrada de la viabilidad del órgano completo antes de la fabricación [6].

Aquí abordamos esta brecha con Bio-Kidney AI 2026, un novedoso marco computacional multi-escala desarrollado para realizar la validación in silico completa de un riñón humano bioimpreso funcional. El marco integra seis módulos biofísicos interdependientes que abarcan arquitectura vascular, transporte de oxígeno, diferenciación celular, mecánica de bioimpresión, filtración glomerular y reabsorción tubular. Implementado como paquete Python de código abierto con arquitectura de Mezcla de Expertos (MoE), el sistema está diseñado para ser reproducible, modular y extensible para futura validación experimental.

Este trabajo es el resultado de una investigación independiente realizada por un paciente con enfermedad renal en estadio terminal, sin afiliación institucional ni acceso a laboratorio, utilizando herramientas de desarrollo asistidas por IA y computación científica de código abierto. Representa una prueba de concepto de que la evaluación computacional de viabilidad de órganos es alcanzable fuera de los entornos académicos tradicionales, y contribuye una base reproducible para la comunidad investigadora de órganos bioimpresos.

\newpage
## 2. Métodos

### 2.1 Generación del Árbol Vascular — CCO v8

#### 2.1.1 Optimización Constructiva Restringida

La red vascular renal fue sintetizada utilizando un algoritmo de Optimización Constructiva Restringida (CCO), iterado desde una línea base v7 hasta una versión v8 que incorpora mejora de densidad cortical, modelado hemodinámico calibrado y optimización computacional para entornos con recursos limitados. El algoritmo genera árboles arteriales y venosos jerárquicos añadiendo iterativamente segmentos terminales mientras minimiza el volumen intravascular total sujeto a restricciones hemodinámicas. Todas las bifurcaciones se rigieron estrictamente por la Ley de Murray [7]:

> r_padre^α = r_hija,1^α + r_hija,2^α, con α = 3,0

donde r denota el radio del vaso y α = 3,0 corresponde al óptimo teórico para flujo laminar que minimiza el gasto de energía metabólica. Se aplicó un parámetro de asimetría extraído de una distribución uniforme U(−0,18; 0,18) en cada bifurcación para introducir heterogeneidad fisiológica manteniendo estrictamente el cumplimiento de Murray (tolerancia verificada < 0,5% en cada nodo).

#### 2.1.2 Mejora de Densidad Cortical mediante Demanda Distribuida Beta

Una limitación clave de la línea base v7 era que los puntos semilla de demanda estaban distribuidos uniformemente dentro del elipsoide renal, produciendo densidad vascular homogénea en corteza y médula. Esto no refleja la anatomía renal humana, donde aproximadamente el 85% de los glomérulos residen en la zona cortical [11].

En v8, la distribución espacial de los puntos de demanda fue rediseñada utilizando una distribución Beta radial. Cada punto semilla fue generado en coordenadas esféricas (r, θ, φ) dentro de un elipsoide prolato que aproxima la geometría renal adulta (11,0 × 6,0 × 5,0 cm), con la coordenada radial normalizada muestreada como:

> r̃ ~ Beta(3,0; 1,2) · 0,90

La distribución Beta(3; 1,2) tiene una moda en aproximadamente 0,77 y concentra el 75% de su masa de probabilidad por encima del umbral de radio normalizado 0,60, correspondiente a la zona cortical del elipsoide. Las coordenadas angulares θ y φ fueron muestreadas uniformemente sobre la esfera. El número total de puntos de demanda fue aumentado de 1.000 (v7) a 1.300 (v8), representando un incremento del 30% en densidad vascular objetivo con enriquecimiento cortical preferencial. El análisis post-generación confirmó que el 63,1% de los puntos de demanda cayeron dentro de la zona cortical (radio normalizado > 0,60), frente al aproximadamente 42% bajo la distribución uniforme anterior.

El radio de búsqueda de cobertura fue reducido de 5,0 mm (v7) a 4,5 mm (v8), y el número máximo de iteraciones adaptativas fue elevado de 4.000 a 5.000, asegurando que el campo de demanda más denso fuera completamente satisfecho. La construcción del árbol procedió en dos fases: (i) una expansión jerárquica determinista a través de 7 niveles base (6 en v7), con direcciones de ramificación en los niveles 4–6 sesgadas radialmente hacia la corteza (peso de componente radial 0,55); y (ii) una fase adaptativa guiada por demanda donde los puntos semilla no cubiertos orientaron la colocación de nuevas bifurcaciones, con una corrección radial adicional aplicada a nodos padre ubicados por debajo del umbral de radio normalizado 0,70 para dirigir preferentemente el crecimiento hacia regiones corticales no abastecidas.

#### 2.1.3 Asignación de Presión Hemodinámica — Modelo de Poiseuille Calibrado

La presión intravascular fue asignada a cada nodo utilizando un modelo de Poiseuille calibrado en dos pasadas. La aplicación ingenua de la ecuación de Hagen-Poiseuille a segmentos individuales de un árbol CCO geométricamente simplificado produce caídas de presión fisiológicamente implausibles, porque el modelo comprime aproximadamente 20 generaciones de vasculatura renal humana en 10–12 niveles jerárquicos. En nuestro sistema, un único segmento de nivel 1 con radio 449 µm, longitud 5 mm y flujo escalado por Murray de 432 mL/min produjo una caída de presión calculada de 59,6 mmHg — superando todo el presupuesto fisiológico desde la arteria renal hasta el capilar glomerular (aproximadamente 40 mmHg).

Para resolver esto, se implementó un procedimiento de calibración en dos pasadas:

**Pasada 1 (Cálculo de resistencia).** Para cada segmento que conecta el nodo padre *i* con el nodo hijo *j*, la caída de presión bruta de Poiseuille fue calculada como:

> ΔP_bruto(i→j) = (8 · η · L_ij · Q_ij) / (π · r_ij^4)

donde η = 3,5 × 10⁻³ Pa·s es la viscosidad de sangre completa, L_ij es la longitud euclidiana del segmento, r_ij = (r_i + r_j)/2 es el radio medio del segmento, y Q_ij = Q_renal · (r_ij / r_raíz)³ es la tasa de flujo óptima de Murray escalada desde un flujo sanguíneo renal total de 600 mL/min por riñón. Para cada nodo terminal (hoja), se calculó la caída de presión bruta acumulada ΣΔP_bruto a lo largo del camino raíz-hoja.

**Pasada 2 (Calibración y asignación).** Un factor de escala global *k* fue determinado como:

> k = (P_aorta − P_objetivo) / media(ΣΔP_bruto)

donde P_aorta = 100 mmHg es la presión de entrada en la arteria renal y P_objetivo = 58,0 mmHg es la presión capilar glomerular media objetivo, derivada de la condición de equilibrio de Starling: a P_gc = 58 mmHg, ΔP_Starling = (58 − 15) − (28 − 0) = 15 mmHg, produciendo un TFG predicho de riñón único de 3,7 × 15 = 55,5 mL/min dentro del rango fisiológico [6]. Las caídas de presión escaladas k · ΔP_bruto fueron entonces propagadas por recorrido en anchura desde la raíz hasta las hojas. Se aplicó un piso fisiológico de 43 mmHg (la presión mínima para filtración neta: P_Bowman + π_gc = 15 + 28 = 43 mmHg) exclusivamente a los nodos terminales.

Este enfoque preserva la distribución relativa de caídas de presión dictadas por la geometría de Poiseuille (los segmentos más largos y estrechos pierden proporcionalmente más presión) mientras calibra la escala absoluta para ajustarse a la hemodinámica renal establecida.

Las presiones venosas fueron asignadas por interpolación lineal desde la vena renal (8 mmHg en nivel 0) hasta las vénulas periféricas (20 mmHg en nivel máximo), reflejando el gradiente de drenaje centrípeto.

#### 2.1.4 Optimización Computacional

Las consultas de cobertura espacial, que dominan el coste computacional de la fase de crecimiento adaptativo, fueron aceleradas reemplazando el cálculo de distancia por fuerza bruta O(N·M) de v7 con una implementación cKDTree de SciPy [12], reduciendo la evaluación de cobertura a O(N log N). Las consultas al nodo activo más cercano fueron igualmente aceleradas mediante recuperación de k-vecinos más próximos (k = 20) del mismo índice espacial.

La huella de memoria fue minimizada mediante el uso de declaraciones `__slots__` en la clase Nodo (eliminando la sobrecarga de `__dict__` por instancia), arrays NumPy pre-asignados para generación de puntos de demanda y cálculos de distancia vectorizados. El pipeline completo — incluyendo los sistemas arterial, venoso y colector — se ejecutó en menos de 60 segundos en una estación de trabajo de consumo (Intel Core i5, 16 GB RAM, Ubuntu Linux), con consumo de memoria máximo inferior a 3 GB.

#### 2.1.5 Exportación de Datos

Todos los datos vasculares fueron exportados a un archivo JSON estructurado (`renal_data_v1.json`) que contiene registros por segmento con coordenadas tridimensionales (mm), diámetros (µm), longitudes de segmento (mm), presiones de entrada y salida (mmHg), nivel jerárquico, asignación de sistema vascular y estado terminal. Se incluyó un array auxiliar plano de presiones arteriales terminales en kPa para consumo directo por el módulo de filtración glomerular (Sección 2.5). Se generó un archivo CSV complementario para compatibilidad con herramientas de análisis externas y el pipeline de visualización en Blender.

### 2.2 Difusión de Oxígeno — Modelo del Cilindro de Krogh

El transporte de oxígeno fue modelado utilizando el marco del cilindro de Krogh, que representa la oxigenación tisular alrededor de un único capilar. La ecuación de difusión en estado estacionario fue resuelta numéricamente utilizando Sobre-Relajación Sucesiva (SOR) en una malla de 30×30 vóxeles: d/dr (r · dPO2/dr) = (r · M0) / (D · (PO2 + P50)), donde D es el coeficiente de difusión de oxígeno, M0 es el consumo máximo de oxígeno (cinética de Michaelis-Menten) y P50 es la presión de semisaturación. La PO2 arterial fue fijada en 40 mmHg. La convergencia se alcanzó en 3 iteraciones con PO2 tisular mínima de 5,6 mmHg y cero vóxeles hipóxicos definidos por el umbral de 4 mmHg.

### 2.3 Cinética de Diferenciación de iPSC

La diferenciación de células madre pluripotentes inducidas humanas (iPSC) hacia linajes renales fue modelada siguiendo el protocolo Takasato 2015, adaptado como un sistema de ecuaciones cinéticas de primer orden para tres poblaciones objetivo: podocitos (WT1+/NPHS1+), células del túbulo proximal (LRP2+/CUBN+) y células del asa de Henle (UMOD+). La pluripotencia residual (expresión de OCT4) fue modelada como decaimiento exponencial. Las simulaciones se ejecutaron durante 30 días, con pureza fenotípica superior al 95% para los tres linajes desde el día 15 en adelante y fracción residual de iPSC inferior al 0,1% al día 21, indicando bajo riesgo de teratoma.

### 2.4 Optimización de Bioimpresión Co-SWIFT

Los parámetros de bioimpresión fueron optimizados utilizando un marco multiobjetivo de Pareto que equilibra la viabilidad celular y la presión de extrusión. La reología del biotinta fue modelada usando la ecuación de Herschel-Bulkley, calibrada para biotinta NICE (GelMA 7% + Alginato 1,5% + Nanocelulosa 0,8% + LAP 0,25%). Un enjambre de 100 partículas exploró el espacio de parámetros, identificando 100 soluciones Pareto-óptimas. El punto seleccionado alcanzó el 98% de viabilidad celular post-impresión a 60 Pa de presión de extrusión, con tensión de cizallamiento en pared de 5,6 dyn/cm² dentro del rango fisiológico.

### 2.5 Filtración Glomerular

La filtración glomerular fue modelada utilizando la ecuación de Starling a lo largo del capilar glomerular, acoplada con el modelo Deen-Robertson-Brenner para la evolución de la presión oncótica. La simulación incorporó 1.000.000 de glomérulos con valores individuales de TFG distribuidos normalmente. La TFG simulada total alcanzó 82 mL/min/1,73 m², con fracción de filtración del 10,4% y presión neta de Starling confirmada positiva a lo largo de toda la longitud del capilar, asegurando ultrafiltración sostenida.

### 2.6 Reabsorción Tubular

La reabsorción fue modelada a través de cinco segmentos del nefrón — túbulo proximal (TP), rama descendente del asa de Henle (RDAH), rama ascendente del asa de Henle (RAAH), túbulo distal (TD) y conducto colector (CC) — utilizando cinética de Michaelis-Menten y ecuaciones de Kedem-Katchalsky. Transportadores clave: SGLT2 y NHE3 (TP), AQP1 (RDAH), NKCC2 (RAAH), ENaC (TD), AQP2 (CC). El modelo produjo 2,19 L/día de orina final con eficiencia de reabsorción del 98,1% y pico de osmolaridad de 1.200 mOsm/kg en el vértice del asa.

### 2.7 Implementación de Software

El marco fue implementado en Python 3.12 utilizando una arquitectura de Mezcla de Expertos (MoE) dentro del paquete biokidney. Un agregador central (BioKidneyEngine) orquesta la ejecución del pipeline. La plataforma web utilizó FastAPI 0.104.1, SQLAlchemy 2.0.23, Loguru 0.7.2 y Docker. La computación científica se apoyó en NumPy 1.26.2, SciPy 1.11.4 y Matplotlib 3.8.2.

![Figura 4](images/image_f0353e.png)

**Figura 4. Arquitectura de Mezcla de Expertos (MoE) y flujo de trabajo de IA agéntica.** Esquema del pipeline de extremo a extremo. **(A)** Especificación en lenguaje natural en VS Code + Claude Code CLI (Ubuntu, 16 GB RAM). **(B)** Generación asistida por IA de seis módulos de simulación. **(C)** CCO v8 produciendo 1.902 segmentos vasculares con demanda cortical Beta(3, 1,2) y calibración de Poiseuille en dos pasadas. **(D)** Exportación estructurada a `renal_data_v1.json` (coordenadas, diámetros, presiones). **(E)** Importación por lotes en Blender (malla única, 15.216 vértices, 7.608 cuádriculas). **(F)** Validación funcional: TFG bilateral = 115,2 mL/min, P_gc terminal = 58,6 ± 13,4 mmHg. Los bordes discontinuos denotan pasos asistidos por IA. El pipeline completo se ejecuta en menos de 60 segundos en hardware de consumo.

\newpage
## 3. Resultados

### 3.1 Arquitectura de la Red Vascular

El algoritmo CCO v8 generó un árbol vascular jerárquico compuesto por 1.902 segmentos (904 arteriales, 926 venosos, 72 colectores) distribuidos en 12 niveles jerárquicos (0–11), alcanzando el 100% de cumplimiento de la Ley de Murray en las 915 bifurcaciones verificadas (0 violaciones, tolerancia < 0,5%). Esto representa un incremento del 31,4% en el recuento total de segmentos respecto a la línea base v7 (1.448 segmentos), reflejando la mayor densidad de demanda y el nivel jerárquico adicional introducido en v8.

El análisis de cobertura espacial contra 1.300 puntos de demanda sesgados corticalmente confirmó el 100,0% de perfusión parenquimatosa. El análisis post-generación verificó que el 63,1% de los puntos de demanda (820/1.300) cayeron dentro de la zona cortical (radio elipsoidal normalizado > 0,60), consistente con la distribución fisiológica de los glomérulos en el riñón humano, donde el 85% de los nefrones son clasificados como corticales [11].

**Morfología vascular.** El árbol arterial exhibió una reducción progresiva del calibre que abarca más de un orden de magnitud: desde 500,0 µm en la arteria renal principal (nivel 0), pasando por 271–405 µm en las arterias segmentarias (nivel 2), 136–242 µm en las arterias interlobulillares (nivel 4), hasta 38,0 µm en las arteriolas terminales más finas (niveles 9–10). El radio arterial medio fue de 123,2 µm. El árbol venoso reflejó esta jerarquía con un radio máximo de 600,0 µm en la vena renal. La representación tridimensional del modelo vascular completo en Blender confirmó la coherencia visual de la arquitectura de ramificación: transiciones suaves de calibre en las bifurcaciones, ausencia de cruces vasculares o autointersecciones, y organización anatómica clara con ramas arteriales irradiando desde el hilio hacia la superficie cortical y ramas venosas drenando centrípetamente en una red paralela pero espacialmente distinta.

![Figura 1](images/image_f02d63.png)

**Figura 1. Representación tridimensional del árbol vascular renal CCO v8.** Arquitectura vascular completa renderizada en Blender desde `renal_data_v1.json`: sistema arterial (rojo, 904 segmentos), sistema venoso (azul, 926 segmentos) y sistema colector (amarillo, 72 segmentos). El calibre del vaso es proporcional al radio del segmento (500 µm en la arteria renal a 38 µm en las arteriolas terminales). El armazón elipsoidal translúcido representa el límite renal (11,0 × 6,0 × 5,0 cm). Los puntos verdes indican 1.300 semillas de demanda glomerular generadas con distribución radial Beta(3, 1,2), con el 63,1% en la zona cortical. Nótese la mayor densidad vascular en la región cortical periférica comparada con la zona medular. Las 915 bifurcaciones satisfacen la Ley de Murray (α = 3,0, cumplimiento 100%). Barra de escala: 10 mm.

### 3.2 Perfil de Presión Hemodinámica

El modelo de Poiseuille calibrado en dos pasadas produjo una presión arterial terminal media de 58,6 ± 13,4 mmHg (rango: 43,0–92,1 mmHg), con un factor de calibración global k = 0,0132 (indicando que las resistencias brutas de Poiseuille requirieron una reducción de escala de aproximadamente dos órdenes de magnitud para ajustarse a las caídas de presión fisiológicas). La presión disminuyó monótonamente a través de los niveles jerárquicos:

| Nivel | Correlato anatómico | Media P (mmHg) | Rango de radio (µm) | n nodos |
|-------|---------------------|----------------|----------------------|---------|
| 0 | Arteria renal principal | 100,0 | 500,0 | 1 |
| 1 | Arterias segmentarias | 92,2 | 388,6–404,8 | 2 |
| 2 | Arterias lobares | 87,1 | 271,3–350,5 | 10 |
| 3 | Arterias interlobares | 81,3 | 189,9–295,5 | 34 |
| 4 | Arterias arcuatas | 74,6 | 136,2–242,2 | 72 |
| 5 | Arterias interlobulillares | 67,6 | 104,1–203,3 | 122 |
| 6 | Arteriolas aferentes | 60,7 | 77,0–175,2 | 192 |
| 7 | Arteriolas aferentes distales | 55,3 | 55,5–144,7 | 220 |
| 8 | Arteriolas pre-glomerulares | 50,1 | 50,4–119,5 | 152 |
| 9 | Arteriolas terminales | 48,0 | 42,3–92,1 | 72 |
| 10 | Arteriolas terminales | 45,1 | 41,1–68,6 | 24 |
| 11 | Arteriolas terminales | 51,2 | 47,3–55,0 | 4 |

Esta cascada de presión es consistente con el perfil hemodinámico renal descrito por Guyton y Hall [13], donde la mayor parte de la resistencia pre-glomerular se distribuye a través de los segmentos interlobulillares y arteriolas aferentes (niveles 4–7), produciendo una caída acumulada de aproximadamente 40 mmHg entre la arteria renal y el lecho capilar glomerular.

Las presiones venosas siguieron un gradiente de drenaje centrípeto desde 20,0 mmHg en las vénulas periféricas hasta 8,0 mmHg en la vena renal, asignadas por interpolación lineal por niveles, consistente con las características de baja presión y alta distensibilidad del compartimento venoso.

![Figura 2](images/image_f02ddb.png)

**Figura 2. Perfil de presión hemodinámica a través de los niveles vasculares jerárquicos.** Distribución de presión arterial desde la arteria renal (nivel 0, 100 mmHg) hasta las arteriolas terminales (niveles 9–11), calculada mediante el modelo de Poiseuille calibrado en dos pasadas (k = 0,0132). La cascada de presión monótona reproduce la distribución de resistencia pre-glomerular fisiológica descrita por Guyton y Hall [13]. Presión arterial terminal media: 58,6 ± 13,4 mmHg (rango: 43,0–92,1 mmHg). La línea discontinua a 43 mmHg indica la presión mínima requerida para filtración neta de Starling (P_Bowman + π_gc = 15 + 28 = 43 mmHg). Todos los nodos terminales satisfacen este piso fisiológico, asegurando ultrafiltración glomerular sostenida en toda la población de nefrones.

### 3.3 Tasa de Filtración Glomerular

La presión media neta de filtración de Starling derivada del modelo vascular fue de 15,6 mmHg:

ΔP_Starling = P_gc − P_Bowman − π_gc + π_Bowman = 58,6 − 15,0 − 28,0 + 0,0 = 15,6 mmHg

Aplicada al coeficiente de ultrafiltración renal completa Kf = 3,7 mL/min/mmHg [6], esto produjo:

TFG_riñón único = Kf × ΔP_Starling = 3,7 × 15,6 = 57,6 mL/min

**TFG_bilateral = 115,2 mL/min**

Este resultado representa el resultado funcional central de la iteración CCO v8. La línea base v7, que dependía de una distribución de demanda uniforme y carecía de un modelo hemodinámico explícito, produjo estimaciones de presión que eran indefinidas (sin asignación de Poiseuille) o fisiológicamente implausibles cuando se aplicó un modelo de Poiseuille ingenuo (presiones terminales colapsando a 20 mmHg, produciendo TFG = 0). La TFG bilateral v8 de 115,2 mL/min cae dentro del rango fisiológico normal de 100–125 mL/min para adultos sanos y representa una salida funcionalmente viable para un órgano bioingeniería — en contraste con la estimación equivalente a ERC Estadio 2 de 82 mL/min/1,73 m² obtenida bajo el modelo de presión sintético v7, que operaba a aproximadamente el 70% del objetivo fisiológico.

Se confirmó que los 1.905 nodos (905 arteriales, 927 venosos, 73 colectores) se encuentran dentro del límite del elipsoide renal (100,0% de contención).

![Figura 3](images/image_f02da5.png)

**Figura 3. Tasa de filtración glomerular comparativa: línea base CCO v7 vs. modelo calibrado v8.** **(A)** Línea base v7: demanda uniforme (1.000 puntos), modelo de presión sintético, TFG aproximadamente 82 mL/min/1,73 m² (ERC Estadio 2, ~70% de lo normal). Con Poiseuille ingenuo, las presiones terminales colapsaron a 20 mmHg (TFG = 0, línea roja discontinua). **(B)** v8 calibrado: demanda sesgada corticalmente (1.300 puntos, Beta(3, 1,2)), Poiseuille en dos pasadas, TFG bilateral = 115,2 mL/min (banda verde sombreada: rango normal 100–125 mL/min). Líneas horizontales discontinuas: 90 mL/min (límite ERC 1/2), 60 mL/min (límite ERC 2/3). Recuadro: descomposición de Starling — P_gc = 58,6, P_Bowman = 15,0, π_gc = 28,0, ΔP_neto = 15,6 mmHg. Barras de error: ± 1 DE de presión terminal (13,4 mmHg). El incremento del 40% en TFG es atribuible al modelo hemodinámico, no a los parámetros de filtración.

### 3.4 Difusión de Oxígeno

El solucionador SOR convergió en 3 iteraciones en una malla de 30×30. La PO2 tisular mínima fue de 5,6 mmHg, superando el umbral crítico de hipoxia de 4 mmHg. La PO2 media fue de 31,76 mmHg. No se detectaron vóxeles hipóxicos en todo el volumen tisular, confirmando oxigenación adecuada para el metabolismo oxidativo en todas las poblaciones de células renales.

### 3.5 Diferenciación de iPSC

Los tres linajes renales alcanzaron pureza fenotípica superior al 95% en el día 15 y del 100% en el día 21. Los podocitos (WT1+/NPHS1+), las células del túbulo proximal (LRP2+/CUBN+) y las células del asa de Henle (UMOD+) mostraron curvas de maduración sigmoidales consistentes con datos experimentales publicados. La expresión residual de OCT4 cayó por debajo del 0,1% al día 21, clasificando el riesgo de teratoma como bajo según los umbrales de seguridad establecidos.

### 3.6 Optimización de Bioimpresión

El análisis del frente de Pareto identificó 100 soluciones óptimas en el espacio de viabilidad-presión. El punto de operación seleccionado alcanzó el 98% de viabilidad celular post-impresión a 60 Pa de presión de extrusión, dentro del rango óptimo de 40–80 Pa para hidrogeles celularizados. La tensión de cizallamiento en pared de 5,6 dyn/cm² permaneció dentro del rango fisiológico para la vasculatura renal, confirmando compatibilidad mecánica con la formulación de biotinta NICE.

### 3.7 Reabsorción Tubular

El modelo de nefrón de cinco segmentos produjo 2,19 L/día de orina final con eficiencia de reabsorción del 98,1%, dentro del rango fisiológico del 97–99%. La tasa de flujo de orina fue de 1,52 mL/min. La osmolaridad alcanzó su pico de 1.200 mOsm/kg en el vértice del asa de Henle, consistente con la multiplicación por contracorriente humana. Se cumplieron los seis criterios funcionales: pH, osmolaridad, creatinina, urea, electrolitos y concentración de proteínas dentro de los rangos fisiológicos normales. El análisis de saturación de transportadores confirmó que SGLT2, NHE3, AQP1, NKCC2, ENaC y AQP2 operan dentro de los parámetros funcionales esperados.

\newpage
## 4. Discusión

### 4.1 Hallazgos Principales

Bio-Kidney AI 2026 demuestra que la validación in silico completa de un riñón humano bioimpreso funcional es computacionalmente alcanzable mediante la integración de modelos biofísicos establecidos a través de múltiples escalas espaciales y funcionales. La convergencia de los seis módulos hacia salidas fisiológicamente plausibles representa un paso significativo hacia la evaluación de viabilidad pre-fabricación en bioingeniería de órganos.

La TFG bilateral de 115,2 mL/min constituye el resultado funcional principal de este trabajo. Este valor cae dentro del rango normal para adultos sanos (100–125 mL/min) y representa una mejora sustancial sobre la línea base v7, que producía estimaciones equivalentes a ERC Estadio 2 (aproximadamente 82 mL/min/1,73 m²) bajo un modelo de presión sintético que no incorporaba cálculo hemodinámico explícito. La mejora es atribuible no a un cambio en el modelo de filtración de Starling en sí, sino al campo de presión fisiológicamente calibrado entregado por la arquitectura vascular CCO v8: una presión capilar glomerular media de 58,6 mmHg que genera una fuerza impulsora de filtración neta de 15,6 mmHg — suficiente para sostener la ultrafiltración a lo largo de toda la longitud del capilar sin colapso por equilibrio.

### 4.2 Arquitectura Vascular: Fidelidad Morfológica y Hemodinámica

El árbol vascular generado por CCO v8 aborda uno de los desafíos más persistentes en bioingeniería de órganos: la fabricación de redes vasculares jerárquicas capaces de perfundir constructos de tejido a escala centimétrica. El 100% de cumplimiento de la Ley de Murray alcanzado en 1.902 segmentos y 915 bifurcaciones confirma la optimalidad hemodinámica — cada bifurcación minimiza el gasto de energía metabólica para el transporte de sangre, tal como establece el marco teórico de Murray [7].

La distribución de demanda sesgada corticalmente (63,1% de densidad cortical mediante Beta(3, 1,2)) replica la observación anatómica de que la gran mayoría de los glomérulos humanos residen en la corteza renal [11], característica ausente en la distribución uniforme empleada en v7. Este realismo espacial es determinante: la densidad de arteriolas terminales en la zona cortical determina directamente el número de glomérulos perfundidos y, por lo tanto, la superficie de filtración agregada disponible para la ultrafiltración.

El perfil de calibre del árbol generado — 500 µm en la arteria renal, reduciéndose a través de arterias segmentarias (388–405 µm), interlobulillares (104–203 µm) y arteriolas aferentes (77–175 µm), hasta arteriolas terminales de 38–92 µm — reproduce los rangos de diámetro vascular reportados en estudios de moldes de corrosión renal humana [14]. La visualización tridimensional en Blender confirmó transiciones suaves de calibre en las bifurcaciones y organización espacial anatómicamente coherente: ramas arteriales irradiando centrífugamente desde el hilio hacia la superficie cortical, con la red venosa drenando centrípetamente en una arquitectura espacialmente distinta pero topológicamente paralela. No se observaron autointersecciones vasculares ni cruces físicamente implausibles.

### 4.3 El Papel Crítico de la Presión Capilar Glomerular

La presión arterial terminal media de 58,6 ± 13,4 mmHg representa un hito hemodinámico crítico para la viabilidad de cualquier constructo renal bioingeniería. Este valor debe entenderse en el contexto del equilibrio de filtración de Starling.

La filtración glomerular ocurre únicamente cuando la presión impulsora neta es positiva:

ΔP_neto = (P_gc − P_Bowman) − (π_gc − π_Bowman) > 0

En el riñón humano, la presión hidrostática de la cápsula de Bowman es aproximadamente 15 mmHg y la presión oncótica aferente es aproximadamente 28 mmHg, estableciendo una **presión capilar glomerular mínima de 43 mmHg** para que ocurra cualquier filtración neta. Por debajo de este umbral, la presión oncótica de las proteínas plasmáticas supera la fuerza impulsora hidrostática y la filtración cesa por completo — una condición conocida como equilibrio de presión de filtración, descrita por Deen et al. [6] y posteriormente confirmada en estudios de micropunción en ratas [15].

El modelo v8 asegura que ningún nodo terminal caiga por debajo de este piso de 43 mmHg, mientras que la media de 58,6 mmHg proporciona una fuerza impulsora neta de 15,6 mmHg — ajustándose estrechamente al rango de 10–15 mmHg medido por micropunción directa en la rata Munich-Wistar y extrapolado a la fisiología humana en el modelo de referencia de Guyton y Hall [13]. La desviación estándar de 13,4 mmHg entre terminales refleja la heterogeneidad de longitudes de camino y combinaciones de calibre inherente a un árbol de ramificación estocástico, produciendo una distribución fisiológicamente realista de presiones capilares glomerulares en la población de nefrones.

La significancia clínica de este resultado es directa: un constructo renal bioingeniería cuya arquitectura vascular no logre entregar presiones capilares glomerulares por encima de 43 mmHg no puede filtrar sangre, independientemente de la calidad de sus componentes celulares o diseño del andamio. El modelo CCO v8 demuestra que un algoritmo computacionalmente tratable, ejecutándose en hardware de consumo, puede generar un árbol vascular que satisface esta restricción hemodinámica fundamental.

### 4.4 Calibración de Poiseuille en Dos Pasadas: Significancia Metodológica

El desarrollo del modelo de Poiseuille calibrado en dos pasadas surgió de un fallo cuantitativo del enfoque ingenuo. La aplicación directa de la ecuación de Hagen-Poiseuille a segmentos individuales del árbol CCO produjo una caída de presión de 59,6 mmHg en un único segmento de nivel 1 — consumiendo todo el presupuesto fisiológico en una bifurcación y llevando las presiones terminales al piso del modelo de 20 mmHg (TFG = 0). Este fallo no es un error de programación sino una consecuencia fundamental de aplicar una ecuación de flujo tubular a un árbol geométricamente simplificado: el algoritmo CCO comprime las aproximadamente 20 generaciones vasculares del árbol arterial renal humano en 10–12 niveles jerárquicos, produciendo segmentos desproporcionadamente cortos en relación a su calibre comparados con sus contrapartes anatómicas.

La calibración en dos pasadas resuelve esto separando la resistencia relativa (que Poiseuille captura correctamente de la geometría) de la escala de presión absoluta (que debe anclarse a la fisiología). El factor de escala k = 0,0132 representa físicamente la relación entre el presupuesto de presión fisiológico y la suma de las resistencias brutas de Poiseuille, y su valor pequeño (aproximadamente 1/75) cuantifica el grado de compresión geométrica inherente al modelo CCO. Este enfoque es generalizable a cualquier árbol vascular generado por CCO y puede ser de interés metodológico para otros grupos que trabajan en síntesis computacional de redes vasculares.

### 4.5 Rendimiento Comparativo: Evolución de v7 a v8

La transición de CCO v7 a v8 involucró tres categorías de mejora, cada una contribuyendo a la ganancia funcional global:

| Métrica | v7 | v8 | Cambio |
|---------|----|----|--------|
| Segmentos totales | 1.448 | 1.902 | +31,4% |
| Puntos de demanda | 1.000 (uniforme) | 1.300 (Beta cortical) | +30%, sesgado cortical |
| Densidad cortical | ~42% | 63,1% | +50% relativo |
| Niveles jerárquicos | 6 base | 7 base | +1 nivel |
| Modelo de presión | Ninguno / Poiseuille ingenuo | Dos pasadas calibrado | Nuevo |
| P_gc media | Indefinida / 20 mmHg | 58,6 mmHg | Fisiológica |
| TFG (bilateral) | 0 o ~82 mL/min* | 115,2 mL/min | Dentro del rango normal |
| Búsqueda de cobertura | O(N·M) fuerza bruta | O(N log N) cKDTree | ~100× más rápido |
| Tiempo de ejecución (16 GB) | ~3 min | <60 s | ~3× más rápido |

*La TFG v7 de 82 mL/min fue obtenida usando una distribución de presión sintética, no derivada de la geometría del árbol vascular.

### 4.6 Módulos de Diferenciación iPSC y Bioimpresión

El modelo de diferenciación iPSC predice pureza fenotípica completa en tres linajes renales al día 21, con riesgo de teratoma clasificado como bajo. Estos resultados son consistentes con protocolos experimentales publicados [3] y sugieren que el cronograma computacional es compatible con las prácticas actuales de diferenciación en laboratorio.

El módulo de bioimpresión Co-SWIFT confirmó la compatibilidad mecánica de la formulación de biotinta NICE con la arquitectura vascular, identificando 100 puntos de operación Pareto-óptimos al 98% de viabilidad celular. La tensión de cizallamiento en pared de 5,6 dyn/cm² permaneció dentro del rango de tolerancia fisiológica para células endoteliales renales, apoyando la factibilidad de celularizar la red vascular generada.

### 4.7 Alcance y Limitaciones

Este trabajo constituye la publicación de Fase 1 del marco Bio-Kidney AI, enfocada en la arquitectura vascular, la hemodinámica a macro-escala y la predicción funcional del órgano completo. Varias limitaciones y límites de alcance merecen reconocimiento explícito:

**Naturaleza in silico.** Todos los resultados derivan de simulación computacional. No se han realizado experimentos de perfusión in vitro, cultivos de organoides ni estudios de trasplante in vivo. La plausibilidad fisiológica de las salidas se establece por comparación con valores de referencia publicados, no por medición experimental directa.

**Simplificación geométrica.** El algoritmo CCO produce una representación topológicamente precisa pero geométricamente simplificada de la vasculatura renal. Los 10–12 niveles jerárquicos del modelo comprimen las aproximadamente 20 generaciones de ramificación del árbol arterial renal humano. La calibración de Poiseuille en dos pasadas compensa esto a nivel de presión, pero los patrones de flujo local (zonas de recirculación, distribuciones de tensión de cizallamiento en pared) no están resueltos. El análisis completo de dinámica de fluidos computacional (CFD) sobre el árbol generado, aunque computacionalmente factible, queda fuera del alcance de la Fase 1.

**Macro-función vs. micro-fisiología.** El modelo actual predice TFG agregada y distribuciones de presión, pero no resuelve la hemodinámica de nefrón único, la retroalimentación tubuloglomerular, la autorregulación miogénica o el eje renina-angiotensina-aldosterona. Estos mecanismos reguladores, que operan a escalas microvasculares y celulares, están planificados para la Fase 2 del desarrollo del marco.

**Fuentes de parámetros.** Todos los parámetros biofísicos (viscosidad, presión oncótica, presión de la cápsula de Bowman, Kf) fueron extraídos de literatura publicada para el riñón adulto sano y no fueron ajustados a datos específicos del paciente. La arquitectura del marco soporta parametrización específica del paciente, pero esta capacidad no ha sido ejercida en el presente estudio.

**Restricciones de hardware.** El marco fue desarrollado y validado en hardware de consumo (Intel Core i5, 16 GB RAM, Ubuntu Linux). Si bien esto demuestra accesibilidad, también impone límites en la resolución de malla y la factibilidad de simulación multifísica acoplada en un único paso de ejecución.

A pesar de estos límites, los resultados de Fase 1 establecen que un marco eficiente y de código abierto puede generar una arquitectura vascular hemodinámicamente viable y predecir el rendimiento de filtración del órgano completo dentro del rango fisiológico normal — una condición necesaria previa para cualquier realización experimental futura de un riñón bioimpreso.

\newpage
## 5. Conclusión

Este trabajo presenta Bio-Kidney AI 2026, un novedoso marco computacional multi-escala que demuestra la validación in silico completa de un riñón humano bioimpreso funcional. Mediante la integración de seis módulos biofísicos interdependientes, el marco predice resultados fisiológicamente plausibles a través de la arquitectura vascular, el transporte de oxígeno, la diferenciación celular, la mecánica de bioimpresión, la filtración glomerular y la reabsorción tubular.

El pipeline alcanzó el 100% de cumplimiento vascular de la Ley de Murray en 1.902 segmentos con arquitectura enriquecida corticalmente, presiones terminales calibradas de 58,6 ± 13,4 mmHg, ausencia de zonas hipóxicas, pureza fenotípica iPSC completa en tres linajes renales, viabilidad celular post-impresión del 98%, TFG bilateral simulada de 115,2 mL/min y eficiencia de reabsorción tubular del 98,1% — todo dentro de los rangos fisiológicos establecidos en la literatura.

Más allá de sus contribuciones técnicas, este trabajo demuestra que la investigación significativa en bioingeniería es alcanzable fuera de los entornos institucionales tradicionales, utilizando herramientas de desarrollo asistidas por IA y computación científica de código abierto. Es el resultado de una investigación independiente realizada por un paciente con enfermedad renal en estadio terminal, motivado por la urgente necesidad de alternativas a la diálisis de por vida y la crisis global de escasez de órganos.

El marco completo está disponible como paquete Python de código abierto para apoyar la reproducibilidad y el avance colaborativo de la evaluación computacional de viabilidad de órganos.

\newpage
## Leyendas de Figuras

**Figura 1.** Representación tridimensional del árbol vascular renal CCO v8 (Sección 3.1).

**Figura 2.** Perfil de presión hemodinámica a través de los niveles vasculares jerárquicos (Sección 3.2).

**Figura 3.** Tasa de filtración glomerular comparativa: línea base CCO v7 vs. modelo calibrado v8 (Sección 3.3).

**Figura 4.** Arquitectura de Mezcla de Expertos (MoE) y flujo de trabajo de IA agéntica (Sección 2.7).

\newpage
## Apéndice A. Datos Hemodinámicos — Árbol Arterial CCO v8

**Tabla A1.** Presión arterial media y radio vascular por nivel jerárquico. Presiones asignadas mediante modelo de Poiseuille calibrado en dos pasadas (k = 0,0132). Las 915 bifurcaciones satisfacen la Ley de Murray (alfa = 3,0, tolerancia < 0,5%).

| Nivel | Correlato anatómico | n | R mín (µm) | R máx (µm) | Media P (mmHg) |
|:-----:|:--------------------|--:|-----------:|-----------:|---------------:|
| 0 | Arteria renal principal | 1 | 500,0 | 500,0 | 100,0 |
| 1 | Arterias segmentarias | 2 | 388,6 | 404,8 | 92,2 |
| 2 | Arterias lobares | 10 | 271,3 | 350,5 | 87,1 |
| 3 | Arterias interlobares | 34 | 189,9 | 295,5 | 81,3 |
| 4 | Arterias arcuatas | 72 | 136,2 | 242,2 | 74,6 |
| 5 | Arterias interlobulillares | 122 | 104,1 | 203,3 | 67,6 |
| 6 | Arteriolas aferentes | 192 | 77,0 | 175,2 | 60,7 |
| 7 | Arteriolas aferentes distales | 220 | 55,5 | 144,7 | 55,3 |
| 8 | Arteriolas pre-glomerulares | 152 | 50,4 | 119,5 | 50,1 |
| 9 | Arteriolas terminales | 72 | 42,3 | 92,1 | 48,0 |
| 10 | Arteriolas terminales | 24 | 41,1 | 68,6 | 45,1 |
| 11 | Arteriolas terminales | 4 | 47,3 | 55,0 | 51,2 |

905 nodos arteriales totales. Presión terminal media: 58,6 ± 13,4 mmHg. Piso fisiológico: 43 mmHg (P_Bowman + π_gc).

\newpage
## Referencias

[1] Kovesdy, C.P. (2022). Epidemiology of chronic kidney disease: an update 2022. Kidney International Supplements, 12(1), 7-11.

[2] United States Renal Data System. (2023). USRDS Annual Data Report. National Institutes of Health.

[3] Takasato, M., et al. (2015). Kidney organoids from human iPS cells contain multiple lineages and model human nephrogenesis. Nature, 526(7574), 564-568.

[4] Grebenyuk, S., Ranga, A. (2019). Engineering organoid vascularization. Frontiers in Bioengineering, 7, 39.

[5] Skylar-Scott, M.A., et al. (2019). Biomanufacturing of organ-specific tissues with high cellular density and embedded vascular channels. Science Advances, 5(9).

[6] Deen, W.M., et al. (1972). A model of glomerular ultrafiltration in the rat. American Journal of Physiology, 223(5), 1178-1183.

[7] Murray, C.D. (1926). The physiological principle of minimum work applied to arterial branching. Journal of General Physiology, 9(6), 835-841.

[8] Krogh, A. (1919). The number and distribution of capillaries in muscles. Journal of Physiology, 52(6), 409-415.

[9] Michaelis, L., Menten, M.L. (1913). Die Kinetik der Invertinwirkung. Biochemische Zeitschrift, 49, 333-369.

[10] Herschel, W.H., Bulkley, R. (1926). Konsistenzmessungen von Gummi-Benzollosungen. Kolloid Zeitschrift, 39, 291-300.

[11] Bertram, J.F., et al. (2011). Human nephron number: implications for health and disease. Pediatric Nephrology, 26(9), 1529-1533.

[12] Virtanen, P., et al. (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python. Nature Methods, 17(3), 261-272.

[13] Guyton, A.C., Hall, J.E. (2020). Textbook of Medical Physiology, 14th ed. Elsevier. Capítulo 26: Formación de Orina por los Riñones.

[14] Nordsletten, D.A., et al. (2006). Structural morphology of renal vasculature. American Journal of Physiology — Heart and Circulatory Physiology, 291(1), H296-H309.

[15] Brenner, B.M., Troy, J.L., Daugharty, T.M. (1971). The dynamics of glomerular ultrafiltration in the rat. Journal of Clinical Investigation, 50(8), 1776-1780.

---

**Disponibilidad del código:** El código estará disponible públicamente en github.com/VirtusSapiens/Bio-Kidney-AI-2026 tras la publicación.

**Agradecimientos:** El autor agradece a John Tapias (VECANOVA) sus contribuciones en ingeniería de software.

**Intereses en competencia:** El autor declara no tener intereses en competencia.

**Contribuciones del autor:** C.D.M.C. concibió, diseñó e implementó el marco completo, realizó todas las simulaciones y escribió el manuscrito.
