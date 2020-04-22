# Load base image

FROM cgt/base:latest

# Copy app files

COPY app.R /srv/shiny-server/
COPY R /srv/shiny-server/R
COPY data /srv/shiny-server/data
COPY www /srv/shiny-server/www

# Permissions
RUN chmod -R +r /srv/shiny-server/


