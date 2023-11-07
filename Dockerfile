FROM python:3.10-slim

RUN apt-get update -y && apt-get install -y git libz-dev build-essential default-jdk wget

# Set the working directory in the container
WORKDIR /usr/src/app
COPY . /usr/src/app
RUN git clone https://bitbucket.org/genomicepidemiology/kma.git
RUN cd kma && make && mv kma /usr/src/app/bin && chmod +x /usr/src/app/bin/kma

ENV PATH=$PATH:/usr/bin:/usr/src/app/bin/

# Install any needed packages specified in requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

