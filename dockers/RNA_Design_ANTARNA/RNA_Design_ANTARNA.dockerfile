FROM python:2.7

WORKDIR /app
COPY anta_rna.yml /app/
RUN pip install --upgrade pip
RUN pip install -r requirements.txt
COPY . /app/

CMD ["python", "main.py"]