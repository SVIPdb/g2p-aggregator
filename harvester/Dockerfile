FROM python:2

WORKDIR /app

RUN apt-get update && apt-get install -y postgresql-client jq

COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

CMD [ "/bin/bash", "entrypoint.sh" ]
