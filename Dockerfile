FROM python:3.12-slim

WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential curl \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY biokidney/ ./biokidney/
COPY web_app/ ./web_app/

RUN mkdir -p web_app/database web_app/logs

ENV PYTHONPATH=/app
ENV DATABASE_URL=sqlite:////app/web_app/database/biokidney.db

EXPOSE 8000

HEALTHCHECK --interval=30s --timeout=5s --start-period=10s --retries=3 \
    CMD curl -f http://localhost:8000/api/system/health || exit 1

CMD ["python", "web_app/backend/main.py"]
