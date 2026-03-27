# Stage 1: Build/Run
FROM python:3.10-slim

WORKDIR /app

# Instalar dependencias del sistema
RUN apt-get update && apt-get install -y \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Copiar requerimientos e instalar
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copiar el core del proyecto y la web app
COPY biokidney/ ./biokidney/
COPY web_app/ ./web_app/

# Crear directorios para BD y Logs
RUN mkdir -p web_app/database web_app/logs

# Exponer puerto FastAPI
EXPOSE 8000

# Variable de entorno para Python Path
ENV PYTHONPATH=/app

# Comando de inicio
CMD ["python", "web_app/backend/main.py"]
