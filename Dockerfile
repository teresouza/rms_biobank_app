FROM rocker/shiny:3.5.2

# System libraries of general use and some required for the R packages
RUN apt-get update && apt-get install -y \
  sudo \
  libssl-dev \
  libxml2-dev \
  libxml2


# Add shiny user
#  RUN groupadd  shiny \
#  && useradd --gid shiny --shell /bin/bash --create-home shiny

# System library dependency for the genetic interaction tool app
RUN R -e "install.packages(c('devtools'))"
RUN R -e "install.packages('shiny', repo='https://mran.microsoft.com/snapshot/2021-05-03/')"
RUN R -e "install.packages('shinyjs', repo='https://mran.microsoft.com/snapshot/2021-05-03/')"
RUN R -e "install.packages('shinyWidgets', repo='https://mran.microsoft.com/snapshot/2021-05-03/')"
RUN R -e "install.packages('shinycssloaders', repo='https://mran.microsoft.com/snapshot/2021-05-03/')"
RUN R -e "install.packages('ggplot2', repo='https://mran.microsoft.com/snapshot/2021-05-03/')"
RUN R -e "install.packages('tidyr', repo='https://mran.microsoft.com/snapshot/2021-05-03/')"
RUN R -e "install.packages('dplyr', repo='https://mran.microsoft.com/snapshot/2021-05-03/')"
RUN R -e "install.packages('tibble', repo='https://mran.microsoft.com/snapshot/2021-05-03/')"
RUN R -e "install.packages('data.table', repo='https://mran.microsoft.com/snapshot/2021-05-03/')"
RUN R -e "install.packages('DT', repo='https://mran.microsoft.com/snapshot/2021-05-03/')"
RUN R -e "install.packages('plotly', repo='https://mran.microsoft.com/snapshot/2021-05-03/')"
RUN R -e "install.packages('ggthemes', repo='https://mran.microsoft.com/snapshot/2021-05-03/')"
RUN R -e "install.packages('ggsci', repo='https://mran.microsoft.com/snapshot/2021-05-03/')"


# Copy the app to the image
WORKDIR /srv/shiny-server
COPY global.R .
COPY app.R .
COPY www www/
COPY data data/

# Setup from R
RUN Rscript global.R

# Allow permission
RUN chown -R shiny:shiny /srv/shiny-server
RUN chown -R shiny:shiny /var/lib/shiny-server

# Document which port can be published
EXPOSE 3838

# Run app
CMD ["/usr/bin/shiny-server.sh"]
