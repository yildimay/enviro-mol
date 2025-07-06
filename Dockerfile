FROM python:3.11-slim

# 1. install deps
WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt \
    && pip install fastapi uvicorn[standard] gunicorn

# 2. add code
COPY green_core ./green_core
COPY backend ./backend
COPY frontend ./frontend

# 3. simple entry: serve static + API
ENV PORT=8000
CMD gunicorn backend.main:app \
    --workers 2 --worker-class uvicorn.workers.UvicornWorker \
    --bind 0.0.0.0:$PORT \
    --access-logfile - \
    --log-level info \
    --forwarded-allow-ips="*"