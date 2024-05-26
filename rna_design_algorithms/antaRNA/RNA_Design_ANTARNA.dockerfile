FROM python:2.7

WORKDIR /app
COPY anta_rna.yml /app/
RUN pip install --upgrade pip
COPY . /app/

CMD ["python", "main.py"]