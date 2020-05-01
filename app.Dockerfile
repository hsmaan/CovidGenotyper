# Load base image

FROM cgt/base:latest

# Copy app files

COPY app.R /srv/shiny-server/
COPY R /srv/shiny-server/R
COPY data /srv/shiny-server/data
COPY www /srv/shiny-server/www

# Permissions
RUN chmod -R +r /srv/shiny-server/

# Expose http port
RUN sed -i 's/3838/80/g' /etc/shiny-server/shiny-server.conf
EXPOSE 80

# Run command
CMD ["/usr/bin/shiny-server.sh"]
