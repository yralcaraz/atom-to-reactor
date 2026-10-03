#!/usr/bin/env python3
"""
Comprehensive synchronization of Linear project (TFM - Uppsala)
with the complete findings and architectural deliverables of improvement_plan.md.

Classification: CONFIDENTIAL (see CLASSIFICATION.md)
Source: Y. Alcaraz Galván (project management)
"""
import os
import requests
import json
from pathlib import Path

ENV_FILE = Path(__file__).resolve().parent.parent / ".env"
API_KEY = os.environ.get("LINEAR_API_KEY")

if not API_KEY and ENV_FILE.exists():
    for line in ENV_FILE.read_text().splitlines():
        if line.startswith("LINEAR_API_KEY="):
            API_KEY = line.split("=", 1)[1].strip().strip('"').strip("'")
            break

if not API_KEY:
    raise ValueError("LINEAR_API_KEY not found in environment or .env file.")

HEADERS = {
    "Authorization": API_KEY,
    "Content-Type": "application/json"
}

def gql(query, variables=None):
    res = requests.post("https://api.linear.app/graphql", headers=HEADERS, json={"query": query, "variables": variables or {}})
    data = res.json()
    if "errors" in data:
        raise RuntimeError(f"GraphQL Error: {data['errors']}")
    return data.get("data", {})

def run_sync():
    print("🚀 Sincronizando y reformulando issues en Linear...")

    # Fetch team, states, labels, milestones, and existing issues
    q = gql('''
    query {
      teams {
        nodes {
          id
          name
          key
          states { nodes { id name type } }
          labels { nodes { id name } }
        }
      }
      project(id: "865d4428-acb3-49ad-83d3-b517df19f04b") {
        id
        name
        state
        projectMilestones { nodes { id name } }
        issues(first: 100) {
          nodes {
            id
            identifier
            title
            state { id name }
            projectMilestone { id name }
          }
        }
      }
    }
    ''')

    team = q["teams"]["nodes"][0]
    team_id = team["id"]
    project = q["project"]
    project_id = project["id"]

    states = {s["name"]: s["id"] for s in team["states"]["nodes"]}
    done_id = states["Done"]
    todo_id = states.get("Todo", states.get("Backlog"))
    backlog_id = states["Backlog"]

    milestones = {m["name"]: m["id"] for m in project["projectMilestones"]["nodes"]}
    labels = {l["name"].lower(): l["id"] for l in team["labels"]["nodes"]}
    existing_issues = {i["identifier"]: i for i in project["issues"]["nodes"]}

    print(f"✓ Conectado a equipo {team['name']}, Proyecto {project['name']}")

    # 1. Update project state
    if project["state"] != "started":
        gql('''
        mutation UpdateProject($id: String!, $input: ProjectUpdateInput!) {
          projectUpdate(id: $id, input: $input) { success }
        }
        ''', {"id": project_id, "input": {"state": "started"}})

    # 2. Mark legacy/solved input tasks as Done with precise notes
    legacy_done = {
        "UNI-70": ("Decidir con Peter: microcinética determinista vs kMC",
                   "**Resuelto:** Se adoptó microcinética determinista de ecuaciones diferenciales rígidas (ODE) integradas con solver implícito Radau IIA de orden 5 con Jacobiano analítico y preservación estricta de invariantes de masa (Si, P)."),
        "UNI-74": ("Confirmar alcance DFT con Peter (reutilización de datos)",
                   "**Resuelto:** Se acordó la reutilización del dataset B3LYP-D3/def2-TZVP de Peter como contrato de entrada cerrado (27 especies químicas) en lugar de recalcular estructuras desde cero."),
        "UNI-75": ("Definir metodología computacional DFT",
                   "**Resuelto:** Metodología consolidada en el contrato de entrada: B3LYP-D3/def2-TZVP para energías electrónicas y correcciones MACE-OMol para energías libres de solvatación."),
        "UNI-76": ("Optimizaciones de geometría (reactivos, intermedios, productos)",
                   "**Resuelto:** Geometrías optimizadas provistas en el dataset de Peter para las 27 especies del sistema."),
        "UNI-77": ("Estados de transición y barreras de activación",
                   "**Resuelto:** En lugar de buscar estados de transición DFT costosos para todas las familias, se adopta el formalismo de escalado lineal de energía libre BEP (Bell-Evans-Polanyi) y Marcus acoplado a la ecuación de Eyring y relaciones de Wegscheider."),
        "UNI-78": ("Compilar dataset de energías de reacción",
                   "**Resuelto:** Dataset cargado en `multiscale_microkinetics.ipynb` con las 27 especies en Hartree y conversión a eV."),
        "UNI-79": ("Seleccionar modelo de solvatación",
                   "**Resuelto:** Modelo de solvatación adoptado: correcciones híbridas MACE-OMol MD con salto de estado estándar (1 atm gas -> 1 M solución, +1.89 kcal/mol)."),
        "UNI-80": ("Aplicar correcciones de solvatación a energías DFT",
                   "**Resuelto:** Implementado en `calculate_solution_gibbs.py` sumando ΔE_solv y el término de concentración."),
        "UNI-81": ("Comprobación cruzada de termodinámica corregida",
                   "**Resuelto:** Verificado en Bloques 3 y 4 del notebook."),
        "UNI-72": ("Definir pasos elementales e intermedios",
                   "**Resuelto:** 9 reacciones elementales (R1 a R9) y 27 intermedios definidos y balanceados estequiométricamente."),
        "UNI-82": ("Formular ecuaciones cinéticas",
                   "**Resuelto:** Leyes de acción de masas reversibles con k_r = k_f / K_eq y velocidades r_j = k_{f,j} ∏ C_r - k_{r,j} ∏ C_p."),
        "UNI-85": ("Definir modelo de reactor (tanque homogéneo)",
                   "**Resuelto:** Reactor por lotes (batch) isotermo perfectamente agitado gobernado por dC/dt = S · r(C, T)."),
        "UNI-86": ("Acoplar microcinética al modelo de reactor",
                   "**Resuelto:** Acoplado en `simulate_tank_reactor.py` usando scipy `solve_ivp` con método Radau IIA.")
    }

    for ident, (title, desc) in legacy_done.items():
        if ident in existing_issues:
            iss = existing_issues[ident]
            gql('''
            mutation UpdateIssue($id: String!, $input: IssueUpdateInput!) {
              issueUpdate(id: $id, input: $input) { success }
            }
            ''', {
                "id": iss["id"],
                "input": {
                    "title": f"[CONTRATO RESUELTO] {title}",
                    "description": desc,
                    "stateId": done_id
                }
            })
            print(f"✓ {ident} actualizado a [CONTRATO RESUELTO] y Done.")

    # 3. Reformulate legacy backlog issues into concrete architectural deliverables
    reformulables = {
        "UNI-83": {
            "title": "[SOLVER] [B7] Implementar Jacobiano analítico y aserciones de conservación de Si/P en solve_ivp",
            "milestone": "Reactor simulations complete",
            "priority": 2,
            "labels": [labels.get("reactor"), labels.get("kinetics")],
            "stateId": todo_id,
            "description": """### Componente: Solver Numérico Radau IIA (Master Architecture Section 4)
Para acelerar la integración de ecuaciones diferenciales rígidas y garantizar invariantes físicas:
- [ ] Implementar el cálculo del **Jacobiano analítico**:
  $$J = S \\cdot \\frac{\\partial r}{\\partial C}$$
  en lugar de calcularlo mediante diferencias finitas aproximadas en `simulate_tank_reactor.py`.
- [ ] Añadir aserciones matemáticas estrictas de conservación de masa tras la integración:
  $$u_{Si}^T \\cdot S = 0, \\quad u_P^T \\cdot S = 0$$
- [ ] Validar que la tolerancia `rtol=1e-8, atol=1e-12` converja sin oscilaciones numéricas."""
        },
        "UNI-87": {
            "title": "[SIM] [B7] Simular escenarios comparativos: Modelo 1 vs 2 vs 3 bajo condiciones reales",
            "milestone": "Reactor simulations complete",
            "priority": 2,
            "labels": [labels.get("reactor"), labels.get("kinetics")],
            "stateId": todo_id,
            "description": """### Componente: Simulación de Escenarios Operativos (Master Architecture Section 4)
Ejecutar simulaciones sistemáticas en el reactor de tanque para evaluar la respuesta del electrolito:
- [ ] **Escenario 1 (Condiciones basales de batería):** T = 298.15 K, [H2O]_0 = 50 ppm, [TMSPA]_0 = 50 mM en EC tamponado (4.5 M).
- [ ] **Escenario 2 (Temperatura extrema):** T = 333.15 K (60 °C, test de degradación acelerada).
- [ ] **Escenario 3 (Sobrecarga de humedad):** [H2O]_0 = 200 ppm (contaminación severa).
- [ ] Generar curvas comparativas de tiempo de vida media de TMSPA y acumulación de H3PO4 / siloxilo."""
        },
        "UNI-88": {
            "title": "[DATA] [B8] Ingesta y estandarización de series temporales de RMN operando del laboratorio",
            "milestone": "Validation & parameter fitting complete",
            "priority": 2,
            "labels": [labels.get("nmr")],
            "stateId": todo_id,
            "description": """### Componente: Ingesta de Datos Experimentales (Master Architecture Section 5)
Preparar el pipeline para recibir las series temporales de RMN operando de 29Si y 31P del laboratorio de Erik:
- [ ] Definir contrato de datos tabular (CSV/JSON) para almacenar pares (tiempo, área de pico integrado, desplazamiento ppm, error experimental).
- [ ] Script de normalización frente a patrón interno o concentración inicial conocida.
- [ ] Rutina de interpolación temporal para emparejar la cuadrícula de simulación con los puntos de adquisición de RMN."""
        },
        "UNI-89": {
            "title": "[VERIF] [B8] Confrontación de perfiles simulados vs integrales experimentales de RMN",
            "milestone": "Validation & parameter fitting complete",
            "priority": 2,
            "labels": [labels.get("nmr"), labels.get("kinetics")],
            "stateId": todo_id,
            "description": """### Componente: Discriminación de Modelos (Master Architecture Section 5)
Comparar cuantitativamente las predicciones de los modelos frente a las observaciones de RMN:
- [ ] Graficar en paneles compartidos las curvas simuladas C_i(t) junto a los puntos experimentales con barras de error.
- [ ] Cuantificar la incapacidad del **Modelo 1** (hidrólisis pura) para reproducir el periodo de inducción o la regeneración de agua.
- [ ] Evaluar si el **Modelo 2** (condensación R4) o **Modelo 3** (transferencia de sililo) reproduce la tasa de emergencia del pico de siloxilo (HMDSO, ~7 ppm)."""
        },
        "UNI-90": {
            "title": "[CALIBRACIÓN] [B9] Optimización inversa de mínimos cuadrados (chi-cuadrado) sobre barreras E0",
            "milestone": "Validation & parameter fitting complete",
            "priority": 2,
            "labels": [labels.get("nmr"), labels.get("kinetics")],
            "stateId": todo_id,
            "description": """### Componente: Calibración Inversa de Parámetros (Master Architecture Section 5)
Ajuste fino de las barreras intrínsecas BEP (E0) de cada familia para minimizar la discrepancia con el laboratorio:
- [ ] Formular función objetivo chi-cuadrado ponderada:
  $$\\chi^2(E_0) = \\sum_{k} \\frac{(I_{sim}(t_k) - I_{exp}(t_k))^2}{\\sigma_k^2}$$
- [ ] Implementar optimizador con restricciones físicas (e.g. Scipy `minimize` con método L-BFGS-B o Nelder-Mead).
- [ ] Calcular intervalos de confianza de los parámetros ajustados mediante la matriz Hessiana del ajuste."""
        }
    }

    for ident, data in reformulables.items():
        if ident in existing_issues:
            iss = existing_issues[ident]
            m_id = milestones.get(data["milestone"])
            filtered_labels = [l for l in data["labels"] if l]
            gql('''
            mutation UpdateIssue($id: String!, $input: IssueUpdateInput!) {
              issueUpdate(id: $id, input: $input) { success }
            }
            ''', {
                "id": iss["id"],
                "input": {
                    "title": data["title"],
                    "description": data["description"],
                    "priority": data["priority"],
                    "stateId": data["stateId"],
                    "projectMilestoneId": m_id,
                    "labelIds": filtered_labels
                }
            })
            print(f"✓ {ident} reformulado como: {data['title']}")

    # 4. Standardize titles and descriptions for UNI-100 to UNI-107
    updates_100_107 = {
        "UNI-100": {
            "title": "[CRÍTICO] [B5] Calibrar barrera E0 de solvent_attack a >= 1.30 eV (Finding 8)",
            "priority": 1,
            "labels": [labels.get("critical"), labels.get("kinetics")],
            "description": """### Contexto Físico & Severidad Crítica (Finding 8)
El diccionario `family_bep_parameters` asigna un valor uniforme de $E_0 = 0.80\\text{ eV}$ a todas las familias de reacción, incluida `solvent_attack` ($R_8, R_9$).

**Evidencia experimental:** Gogoi et al. 2024 demuestra que el TMSOH **no reacciona** con carbonato de etileno (EC) a temperatura ambiente. Con $E_0 = 0.80\\text{ eV}$, $R_8$ y $R_9$ avanzan $\\sim 10^5$ veces más rápido de lo real, causando:
1. Desgasificación masiva y ficticia de $\\text{CO}_2$ en los Bloques 7 y 8.
2. Consumo prematuro de TMSOH antes de que pueda intervenir en las reacciones clave de transferencia de sililo ($R_5\\text{--}R_7$).
3. Bloqueo de la cinética observada por Peter en el laboratorio.

### Tareas
- [ ] Modificar `family_bep_parameters['solvent_attack']['E0_eV'] = 1.30` en `calculate_rate_constants.py`.
- [ ] Actualizar la celda del Bloque 5 en `multiscale_microkinetics.ipynb`.
- [ ] Re-ejecutar el reactor del Bloque 7 y comprobar que TMSOH sobrevive y participa en la química real."""
        },
        "UNI-101": {
            "title": "[BUG] [B7] Corregir transitorio inicial en integración ODE a t0 = 1e-3 s (Finding 10)",
            "priority": 2,
            "labels": [labels.get("bug"), labels.get("reactor")],
            "description": """### Defecto Numérico en Reactor (Finding 10)
En `simulate_tank_reactor.py`:
`t_eval = np.logspace(0, np.log10(t_end_s), n_points)`
fija el inicio de evaluación en $10^0 = 1.0\\text{ s}$.

Como la constante de hidrólisis inicial es rápida ($k_f \\sim 0.2\\text{ s}^{-1}$), en el primer segundo transcurre hasta un $20\\%$ de la conversión total de agua. Al omitir el rango $t < 1\\text{ s}$, las curvas cinéticas pierden la fase inicial de choque.

### Tareas
- [ ] Modificar el mallado temporal a `np.logspace(-3, np.log10(t_end_s), n_points)` en `simulate_tank_reactor.py`.
- [ ] Actualizar gráficos de evolución de especies para capturar la escala sub-segundo ($10^{-3}\\text{ s}$ a $1\\text{ s}$)."""
        },
        "UNI-102": {
            "title": "[VERIF] [B4] Implementar y verificar los 3 ciclos independientes de Wegscheider (Finding 5)",
            "priority": 2,
            "labels": [labels.get("thermodynamics"), labels.get("kinetics")],
            "description": """### Integridad Termodinámica Cíclica (Finding 5)
La red de 9 reacciones contiene **3 ciclos de Hess independientes**:
1. $R_1 + R_4 - R_5 = 0$ (Actualmente probado en el notebook)
2. $R_2 + R_4 - R_6 = 0$ (No verificado explícitamente en el código)
3. $R_3 + R_4 - R_7 = 0$ (No verificado explícitamente en el código)

### Tareas
- [ ] Añadir en `calculate_reaction_thermo.py` la verificación numérica de los 3 ciclos:
  $$\\Delta G(R_2) + \\Delta G(R_4) - \\Delta G(R_6) = 0$$
  $$\\Delta G(R_3) + \\Delta G(R_4) - \\Delta G(R_7) = 0$$
- [ ] Añadir aserciones equivalentes para las razones cinéticas de velocidades $k_f / k_r$."""
        },
        "UNI-103": {
            "title": "[DATA] [B1] Sustituir frecuencias e inercias placeholder por datos reales de Peter para qRRHO (Finding 1)",
            "priority": 3,
            "labels": [labels.get("thermodynamics")],
            "description": """### Defecto en Modo Termodinámico qRRHO (Finding 1)
Actualmente, las 27 especies químicas comparten **exactamente los mismos momentos de inercia dummy** (`(150, 200, 300)`) y el mismo espectro ficticio de 8 frecuencias (`[45, 85, 150, 300, 650, 1050, 1260, 2960] cm⁻¹`).

Esto causa que el modo `qRRHO` produzca entropías y $G_{sol}(T)$ físicamente erróneas (razón por la cual el notebook se apoya en el modo benchmark directo).

### Tareas
- [ ] Extraer del archivo de salida DFT de Peter los modos normales de vibración reales para cada especie.
- [ ] Extraer los momentos principales de inercia reales a partir de las coordenadas nucleares optimizadas.
- [ ] Poblar la base de datos de especies con los valores reales y validar el modo `qRRHO`."""
        },
        "UNI-104": {
            "title": "[MODEL] [B4/B7] Formalizar matrices estequiométricas modulares para Modelos KIN 1 a 4",
            "priority": 2,
            "labels": [labels.get("kinetics")],
            "description": """### Espacio de Hipótesis de Redes Cinéticas (Master Architecture Section 2)
Implementar la formulación modular que permita alternar entre 4 topologías de red en el reactor:
- **Modelo 1 [Hidrólisis secuencial]:** $R_1\\text{--}R_3$ (sumidero irreversible de agua).
- **Modelo 2 [Hidrólisis + Condensación]:** Añade $R_4$ ($2\\text{ TMSOH} \\rightleftharpoons \\text{siloxyl} + \\text{H}_2\\text{O}$, reciclaje catalítico).
- **Modelo 3 [Transferencia de sililo + Disolvente]:** Añade $R_5\\text{--}R_7$ y $R_8\\text{--}R_9$ (reacción con EC).
- **Modelo 4 [Red extendida]:** Competencia de aditivos (VC) y evolución de volátiles.

### Tareas
- [ ] Formalizar las matrices estequiométricas $\\mathbf{S}_1, \\mathbf{S}_2, \\mathbf{S}_3, \\mathbf{S}_4$ en un módulo desacoplado.
- [ ] Validar conservación estequiométrica en cada una de ellas."""
        },
        "UNI-105": {
            "title": "[SWEEP] [B9] Implementar barridos de sensibilidad multi-paramétricos por familia (Finding 14)",
            "priority": 3,
            "labels": [labels.get("reactor"), labels.get("kinetics")],
            "description": """### Sensibilidad Paramétrica Desacoplada (Finding 14)
El Bloque 9 varía un escalar global $E_0$ idéntico para todas las familias simultáneamente.
Para calibrar el modelo frente al laboratorio, se requiere variar las barreras de forma desacoplada:
- $\\Delta E_0$ para `hydrolysis`
- $\\Delta E_0$ para `condensation`
- $\\Delta E_0$ para `transfer`
- $\\Delta E_0$ para `solvent_attack`

### Tareas
- [ ] Refactorizar la función de barrido en el Bloque 9 para aceptar variaciones individuales o matriciales por familia.
- [ ] Generar mapa de calor (heatmap) de sensibilidad de la vida media de TMSPA."""
        },
        "UNI-106": {
            "title": "[RMN] [B8] Tratar degeneración espectral BMSPA/TMSOEG a -18.01 ppm en RMN virtual (Finding 13)",
            "priority": 2,
            "labels": [labels.get("nmr")],
            "description": """### Degeneración Espectroscópica en ²⁹Si RMN (Finding 13)
BMSPA y TMSOEG tienen asignado exactamente el mismo desplazamiento químico en $^{29}\\text{Si}$ ($-18.01\\text{ ppm}$).
En el espectro experimental de 1D, estos picos se solapan perfectamente y no son deconvolucionables.

### Tareas
- [ ] Implementar en `simulate_virtual_nmr.py` el espectro virtual sintético de $^{31}\\text{P}$ (donde BMSPA sí se distingue de TMSOEG al no tener fósforo este último).
- [ ] Documentar en la teoría cómo la combinación $^{29}\\text{Si} + ^{31}\\text{P}$ resuelve la ambigüedad espectral."""
        }
    }

    for ident, data in updates_100_107.items():
        if ident in existing_issues:
            iss = existing_issues[ident]
            filtered_labels = [l for l in data["labels"] if l]
            gql('''
            mutation UpdateIssue($id: String!, $input: IssueUpdateInput!) {
              issueUpdate(id: $id, input: $input) { success }
            }
            ''', {
                "id": iss["id"],
                "input": {
                    "title": data["title"],
                    "description": data["description"],
                    "priority": data["priority"],
                    "labelIds": filtered_labels
                }
            })
            print(f"✓ {ident} estandarizado: {data['title']}")

    # 5. Create missing findings from improvement_plan.md
    missing_issues = [
        {
            "title": "[THEORY] [B1] Corregir heurística de números de simetría molecular sigma_rot (Finding 2)",
            "milestone": "Microkinetic model complete",
            "priority": 4, # Low
            "labels": [labels.get("thermodynamics")],
            "description": """### Inexactitud en Números de Simetría (Finding 2)
En el Bloque 1 de `multiscale_microkinetics.ipynb`:
`'sigma_rot': 1 if 'TMS' in name else 2`
asigna $\\sigma_{rot} = 1$ a cualquier molécula con TMS y $\\sigma_{rot} = 2$ al resto.

**Errores identificados:**
- $\\text{CH}_4$ tiene simetría tetraédrica $T_d$ ($\\sigma = 12$), pero recibe $\\sigma = 2$.
- $\\text{H}_2\\text{O}$ ($C_{2v}$, $\\sigma = 2$) y $\\text{CO}_2$ (lineal con centro de inversión, $\\sigma = 2$) coinciden por casualidad.

### Tareas
- [ ] Sustituir la heurística por una tabla explícita de números de simetría grupo-puntual por especie.
- [ ] Validar que la entropía rotacional en modo qRRHO refleje los valores correctos."""
        },
        {
            "title": "[THEORY] [B3] Documentar e incorporar contribución entrópica en solvatación ΔS_solv (Finding 4)",
            "milestone": "Solvation corrections complete",
            "priority": 3, # Medium
            "labels": [labels.get("thermodynamics")],
            "description": """### Suposición Físico-Química en Solvatación (Finding 4)
En `calculate_solution_gibbs.py`, la energía libre de solvatación se calcula como:
$$G_{sol} = G_{gas}^\\circ + \\Delta E_{solv} + \\Delta G^{\\circ \\to *}$$
El término $\\Delta E_{solv}$ (obtenido de MACE-OMol) se trata como puramente energético con entropía de solvatación nula ($\\Delta S_{solv} = 0$).

En disolución de carbonato de etileno, la reordenación de la primera capa de solvatación conlleva costos entrópicos significativos.

### Tareas
- [ ] Documentar explícitamente en `theory-multiscale_microkinetics.md` el impacto de omitir $\\Delta S_{solv}$.
- [ ] Evaluar correcciones empíricas de entropía de cavitación para líquidos polares."""
        },
        {
            "title": "[DOC] [B4] Añadir advertencia de datos placeholder en comparativa qRRHO del Bloque 4 (Finding 6)",
            "milestone": "Microkinetic model complete",
            "priority": 4, # Low
            "labels": [labels.get("thermodynamics")],
            "description": """### Advertencia de Claridad en Notebook (Finding 6)
La tabla del Bloque 4 muestra lado a lado los valores termodinámicos de modo `benchmark` vs `qRRHO`.
Dado que qRRHO utiliza actualmente frecuencias dummy (Finding 1), las columnas de $\\Delta(\\Delta G)$ y $K_{eq}$ en qRRHO no tienen sentido físico real.

### Tareas
- [ ] Añadir un banner `> [!WARNING]` en el markdown previo a la tabla en `multiscale_microkinetics.ipynb` indicando que los datos de qRRHO son ilustrativos hasta que se cargue el Hessiano real."""
        },
        {
            "title": "[DOC] [B6] Clarificar la naturaleza de ajuste Arrhenius R²=1.0000 en benchmark (Finding 9)",
            "milestone": "Microkinetic model complete",
            "priority": 4, # Low
            "labels": [labels.get("kinetics")],
            "description": """### Advertencia Metodológica (Finding 9)
En la regresión de Arrhenius modificada:
$$\\ln k = \\ln A + \\beta \\ln(T/T_0) - \\frac{E_a}{RT}$$
El modelo en modo benchmark obtiene $R^2 = 1.0000$ de forma exacta por construcción, debido a que proviene de la ecuación de Eyring con barreras independientes de la temperatura.

### Tareas
- [ ] Clarificar en el notebook y en la teoría que este ajuste perfecto es una consecuencia algebraica y no una validación física independiente.
- [ ] Explicar que la curvatura genuina no-Arrhenius solo aparecerá al acoplar la dependencia térmica de $\\Delta G^\\ddagger(T)$ con qRRHO real."""
        },
        {
            "title": "[REFACTOR] [B7] Heredar selected_mode global en simulate_tank_reactor (Finding 11)",
            "milestone": "Reactor simulations complete",
            "priority": 4, # Low
            "labels": [labels.get("reactor")],
            "description": """### Consistencia de Parámetros en Solver (Finding 11)
La función `simulate_tank_reactor` tiene por defecto `mode='wb97mv'`.
Si un usuario modifica la variable global `selected_mode = 'qRRHO'` en las celdas superiores, el simulador del reactor continuará silenciosamente usando el modo benchmark salvo que se le pase explícitamente el argumento.

### Tareas
- [ ] Modificar la llamada en el Bloque 7 para que pase explícitamente `mode=selected_mode`.
- [ ] Añadir una comprobación en el simulador para advertir qué modo termodinámico está conduciendo la integración."""
        }
    ]

    existing_titles = {i["title"].strip() for i in existing_issues.values()}

    for item in missing_issues:
        if item["title"].strip() in existing_titles:
            print(f"↷ Omitiendo (ya existe): {item['title']}")
            continue

        m_id = milestones.get(item["milestone"])
        filtered_labels = [l for l in item["labels"] if l]

        res = gql('''
        mutation CreateIssue($input: IssueCreateInput!) {
          issueCreate(input: $input) {
            success
            issue {
              id
              identifier
              title
            }
          }
        }
        ''', {
            "input": {
              "teamId": team_id,
              "projectId": project_id,
              "projectMilestoneId": m_id,
              "title": item["title"],
              "description": item["description"],
              "priority": item["priority"],
              "stateId": todo_id,
              "labelIds": filtered_labels
            }
        })
        new_i = res["issueCreate"]["issue"]
        print(f"✓ Creada issue [{new_i['identifier']}]: {new_i['title']}")

    print("\n🎉 Todas las issues han sido reformuladas y alineadas con el implementation_plan.")

if __name__ == "__main__":
    run_sync()
