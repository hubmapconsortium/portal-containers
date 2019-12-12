FROM python:3.8.0-alpine

COPY . .
RUN pip install -r ./requirements.txt

CMD [ "python", "main.py", "--input_dir", "/input", "--output_dir", "/output" ]
