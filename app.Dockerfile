# Load base image

FROM cgt/base:latest

# Copy app files

COPY app.R /srv/shiny-server/
COPY R /srv/shiny-server/R
COPY data /srv/shiny-server/data

# Select port and add permission (LEFT OFF FOR NOW - USING googleComputeEngineR)

# EXPOSE 3838
# RUN sudo chown -R shiny:shiny /srv/shiny-server

# Run app

# CMD ["/usr/bin/shiny-server.sh"]
