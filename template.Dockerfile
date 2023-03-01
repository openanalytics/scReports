#include packamon.disclaimer

#include packamon.from

RUN apt-get update && apt-get install -y \
  libfontconfig1-dev \
  libxt-dev \
  python3 \
  python3-pip

#include packamon.system-dependencies

RUN apt-get update && apt-get install -y \
    patch
    
RUN R -e "options(warn=2)"

#include packamon.r-repos

#include packamon.r-dependencies

#include packamon.local-r-dependencies

#include packamon.runtime-settings

RUN pip install --no-cache-dir "anndata==0.7.8"
RUN pip install --no-cache-dir "scanpy==1.8.2"
