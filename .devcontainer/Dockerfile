FROM rocker/verse:4.5.3

RUN usermod -l ed rstudio \
    && usermod -d /home/ed -m ed \
    && groupmod -n ed rstudio \
    && mkdir -p /home/ed/.claude /renv/cache /renv/library \
    && chown -R ed:ed /home/ed/.claude /renv/cache /renv/library

RUN apt-get update \
    && apt-get install -y python3-pip ca-certificates curl gnupg \
        libudunits2-dev libgdal-dev libgeos-dev libproj-dev \
    && pip3 install radian --break-system-packages \
    && mkdir -p /etc/apt/keyrings \
    && curl -fsSL https://deb.nodesource.com/gpgkey/nodesource-repo.gpg.key | gpg --dearmor -o /etc/apt/keyrings/nodesource.gpg \
    && echo "deb [signed-by=/etc/apt/keyrings/nodesource.gpg] https://deb.nodesource.com/node_18.x nodistro main" > /etc/apt/sources.list.d/nodesource.list \
    && apt-get update \
    && apt-get install -y nodejs \
    && npm install -g @anthropic-ai/claude-code \
    && rm -rf /var/lib/apt/lists/*
