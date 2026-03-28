# BioKidney AI 2026 — Guia de Despliegue Local

Plataforma web de validacion in-silico para rinon bioimpreso.
6 modulos: Vascular CCO, Filtracion Glomerular, Reabsorcion Tubular, Oxigenacion, iPSC, Bioprinting.

---

## Requisitos previos

| Componente      | Version minima | Verificar con          |
|-----------------|---------------|------------------------|
| Python          | 3.10+         | `python3 --version`    |
| pip             | 22+           | `pip --version`        |
| Docker          | 24+           | `docker --version`     |
| Docker Compose  | 2.20+         | `docker compose version` |
| Git             | 2.30+         | `git --version`        |

> Docker y Docker Compose solo son necesarios para el despliegue con contenedores (opcion B).

---

## Opcion A: Despliegue local con Python (desarrollo)

### 1. Clonar el repositorio

```bash
git clone <url-del-repositorio> BioKidney-AI
cd BioKidney-AI
```

Si ya tienes el proyecto:

```bash
cd ~/Escritorio/BioKidney-AI
```

### 2. Crear entorno virtual

```bash
python3 -m venv .venv
```

### 3. Activar el entorno

**Linux / macOS:**
```bash
source .venv/bin/activate
```

**Windows (PowerShell):**
```powershell
.venv\Scripts\Activate.ps1
```

**Windows (CMD):**
```cmd
.venv\Scripts\activate.bat
```

### 4. Instalar dependencias

```bash
pip install -r requirements.txt
```

Verificar que se instalaron correctamente:

```bash
python -c "import fastapi, uvicorn, numpy, scipy, matplotlib, sqlalchemy, loguru; print('OK')"
```

### 5. Crear directorios de datos

```bash
mkdir -p web_app/database web_app/logs
```

### 6. Iniciar el servidor

```bash
PYTHONPATH=. python web_app/backend/main.py
```

En Windows:

```cmd
set PYTHONPATH=.
python web_app/backend/main.py
```

### 7. Abrir la aplicacion

Abrir en el navegador:

```
http://localhost:8000
```

### 8. Verificar que funciona

```bash
curl http://localhost:8000/api/system/health
```

Respuesta esperada:

```json
{
  "status": "Healthy",
  "version": "2026.2.0",
  "engine": "BioKidney-Aggregator-MoE-v2",
  "database": "SQLite-Active",
  "modules": ["vascular", "filtration", "reabsorption", "oxygen", "ipsc", "bioprinting"]
}
```

### 9. Ejecutar una simulacion de prueba

```bash
curl -X POST http://localhost:8000/api/simulate/filtration \
  -H "Content-Type: application/json" \
  -d '{"params": {"pgc": 60}}'
```

### 10. Detener el servidor

Presionar `Ctrl+C` en la terminal donde se inicio.

---

## Opcion B: Despliegue con Docker Compose (produccion)

### 1. Construir la imagen

```bash
cd ~/Escritorio/BioKidney-AI
docker compose build
```

Esto crea la imagen `biokidney-ai:latest` (~400 MB).

### 2. Iniciar el contenedor

```bash
docker compose up -d
```

El flag `-d` ejecuta en segundo plano (detached).

### 3. Verificar que esta corriendo

```bash
docker compose ps
```

Salida esperada:

```
NAME              STATUS          PORTS
biokidney-ai     Up (healthy)    0.0.0.0:8000->8000/tcp
```

### 4. Ver logs en tiempo real

```bash
docker compose logs -f
```

### 5. Abrir la aplicacion

```
http://localhost:8000
```

### 6. Detener el servicio

```bash
docker compose down
```

### 7. Detener y eliminar datos persistidos

```bash
docker compose down -v
```

> Esto elimina los volumenes (base de datos e historial de simulaciones).

---

## Estructura de archivos relevantes

```
BioKidney-AI/
├── biokidney/                  # Core del framework (MoE)
│   ├── core/config.py          # Parametros fisiologicos
│   ├── experts/                # Expertos: vascular, fluids, cellular
│   └── aggregator.py           # Orquestador del pipeline
├── web_app/
│   ├── backend/
│   │   ├── main.py             # FastAPI — endpoints
│   │   ├── services/
│   │   │   └── simulation_service.py  # Logica de simulacion
│   │   ├── database/
│   │   │   └── models.py       # SQLAlchemy ORM
│   │   └── utils/
│   │       └── logger.py       # Loguru logging
│   ├── frontend/
│   │   └── dashboard.html      # SPA (Tailwind + Chart.js)
│   ├── database/               # SQLite (auto-creado)
│   └── logs/                   # Logs de auditoria
├── Dockerfile
├── docker-compose.yml
├── .dockerignore
└── requirements.txt
```

---

## Endpoints de la API

| Metodo | Ruta                          | Descripcion                       |
|--------|-------------------------------|-----------------------------------|
| GET    | `/`                           | Dashboard SPA                     |
| POST   | `/api/simulate`               | Pipeline completo (6 modulos)     |
| POST   | `/api/simulate/vascular`      | Vascular CCO + Murray             |
| POST   | `/api/simulate/filtration`    | Filtracion Glomerular (Starling)  |
| POST   | `/api/simulate/reabsorption`  | Reabsorcion Tubular (5 segmentos) |
| POST   | `/api/simulate/oxygen`        | Difusion de Oxigeno (Krogh)       |
| POST   | `/api/simulate/ipsc`          | Diferenciacion iPSC (3 protocolos)|
| POST   | `/api/simulate/bioprinting`   | Co-SWIFT Bioprinting (Pareto)     |
| GET    | `/api/history`                | Historial de simulaciones         |
| GET    | `/api/system/health`          | Estado del sistema                |

Todos los endpoints POST aceptan:

```json
{
  "params": {
    "clave": "valor"
  }
}
```

---

## Parametros configurables por modulo

### Vascular CCO
| Parametro        | Default  | Rango       |
|-----------------|----------|-------------|
| `n_seeds`       | 1000     | 100 - 5000  |
| `murray_exponent`| 3.0     | 2.0 - 4.0   |

### Filtracion Glomerular
| Parametro        | Default  | Rango       |
|-----------------|----------|-------------|
| `pgc`           | 60 mmHg  | 40 - 80     |
| `pi_gc`         | 28 mmHg  | 15 - 40     |
| `n_glomerulos`  | 1000000  | 100k - 2M   |

### Reabsorcion Tubular
| Parametro            | Default | Rango     |
|---------------------|---------|-----------|
| `gfr_input`         | 82      | 30 - 150  |
| `aldosterone_factor`| 1.0     | 0 - 3     |
| `adh_factor`        | 1.0     | 0 - 3     |

### Oxigenacion
| Parametro    | Default | Rango       |
|-------------|---------|-------------|
| `grid_size` | 40      | 15 - 60    |
| `p_art_o2`  | 40 mmHg | 20 - 100   |

### iPSC
| Parametro | Default | Rango   |
|----------|---------|---------|
| `t_end`  | 30 dias | 15 - 60 |

### Bioprinting
| Parametro      | Default | Rango     |
|---------------|---------|-----------|
| `n_particles` | 100     | 20 - 500  |

---

## Solucion de problemas

### El servidor no inicia

**Error: `ModuleNotFoundError: No module named 'web_app'`**

Asegurate de establecer `PYTHONPATH`:

```bash
PYTHONPATH=. python web_app/backend/main.py
```

**Error: `unable to open database file`**

Crear el directorio de datos:

```bash
mkdir -p web_app/database web_app/logs
```

### Puerto 8000 ocupado

```bash
# Ver que proceso usa el puerto
lsof -i :8000

# Detenerlo
kill $(lsof -t -i:8000)
```

O cambiar el puerto editando `web_app/backend/main.py` linea final:

```python
uvicorn.run(app, host="0.0.0.0", port=9000)
```

### Docker: imagen no construye

Verificar que `.dockerignore` existe y que no se esta copiando `.venv/`:

```bash
cat .dockerignore
```

### Docker: contenedor no arranca

```bash
docker compose logs biokidney-platform
```

---

## Variables de entorno

| Variable       | Default                                       | Descripcion              |
|---------------|-----------------------------------------------|--------------------------|
| `PYTHONPATH`  | `.` (local) o `/app` (Docker)                 | Ruta de modulos Python   |
| `DATABASE_URL`| `sqlite:///web_app/database/biokidney.db`     | Conexion a base de datos |

---

## Resultados esperados del pipeline completo

| Modulo        | Metrica principal            | Valor esperado    | Estado  |
|--------------|------------------------------|-------------------|---------|
| Vascular     | Murray compliance            | 100%              | OPTIMO  |
| Filtracion   | TFG                          | ~65 mL/min        | FUNCIONAL |
| Reabsorcion  | Tasa reabsorcion             | ~98.7%            | FUNCIONAL |
| Oxigeno      | Zonas hipoxicas              | 0%                | OPTIMO  |
| iPSC         | Protocolos con pureza >95%   | 3/3               | OPTIMO  |
| Bioprinting  | Viabilidad celular           | ~98%              | OPTIMO  |

---

*BioKidney AI 2026 — VirtusSapiens — Carlos David Moreno Caceres*
