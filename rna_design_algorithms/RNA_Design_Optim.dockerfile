FROM python:3.8

WORKDIR /app
COPY environment.yml /app/
COPY requirements.txt /app/
RUN pip install --upgrade pip
COPY . /app/