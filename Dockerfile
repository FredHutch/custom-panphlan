FROM ubuntu:16.04
MAINTAINER Sam Minot sminot@fredhutch.org

# Install pre-requisites
RUN apt update && \
    apt install -y python3 python3-pip

# Small utility to facilitate execution via sciluigi
RUN pip3 install bucket_command_wrapper==0.3.0 pandas

ADD scripts/get_patric_data.py /usr/local/bin/
