FROM python:3.8

WORKDIR /app
COPY environment.yml /app/
COPY requirements.txt /app/
RUN pip install --upgrade pip
RUN pip install -r requirements.txt
COPY . /app/

CMD ["python", "main.py"]