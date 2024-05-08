FROM python:3.9-slim

WORKDIR /app
# python requirements (including ViennaRNA's python3 binding)
COPY requirements.txt .
RUN pip install -r requirements.txt

WORKDIR /work
COPY . .
CMD ["/bin/bash"]
