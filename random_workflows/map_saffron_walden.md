# making a map of Saffron Walden


- using this as inspiration: http://joshuamccrain.com/tutorials/maps/streets_tutorial.html
```R


library(remotes)
remotes::install_github("ropensci/osmdata")
library(tidyverse)
library(osmdata) # package for working with streets
library(showtext) # for custom fonts
library(ggmap)
library(rvest)



big_streets <- getbb("Saffron Walden United Kingdom")%>%
  opq()%>%
  add_osm_feature(key = "highway",
                  value = c("motorway", "primary", "motorway_link", "primary_link")) %>%
  osmdata_sf()

med_streets <- getbb("Saffron Walden United Kingdom")%>%
  opq()%>%
  add_osm_feature(key = "highway",
                  value = c("secondary", "tertiary", "secondary_link", "tertiary_link")) %>%
  osmdata_sf()


small_streets <- getbb("Saffron Walden United Kingdom")%>%
  opq()%>%
  add_osm_feature(key = "highway",
                  value = c("residential", "living_street",
                            "unclassified",
                            "service"
                  )) %>%
  osmdata_sf()

bridleway <- getbb("Saffron Walden United Kingdom")%>%
    opq()%>%
    add_osm_feature(key = "highway",
                    value = c("bridleway","footway","cycleway","track")) %>%
    osmdata_sf()

river <- getbb("Saffron Walden United Kingdom")%>%
  opq()%>%
  add_osm_feature(key = "waterway", value = c("river","stream","ditch")) %>%
  osmdata_sf()

railway <- getbb("Saffron Walden United Kingdom")%>%
  opq()%>%
  add_osm_feature(key = "railway", value=c("rail")) %>%
  osmdata_sf()
railway2 <- getbb("Saffron Walden United Kingdom")%>%
    opq()%>%
    add_osm_feature(key = "railway", value=c("platform")) %>%
    osmdata_sf()



# add fonts
font_add_google(name = "Lato", family = "lato")
showtext_auto()


# make the plot
PLOT <- ggplot() +
     geom_sf(data = river$osm_lines,
          inherit.aes = FALSE,
          color = "steelblue",
          size = 2,
          alpha = .5) +
     geom_sf(data = railway$osm_lines,
          inherit.aes = FALSE,
          color = "black",
          size = 1,
          linetype="dotted",
          alpha = .5) +
     geom_sf(data = med_streets$osm_lines,
          inherit.aes = FALSE,
          color = "black",
          size = .3,
          alpha = .6) +
     geom_sf(data = small_streets$osm_lines,
          inherit.aes = FALSE,
          color = "#666666",
          size = .2,
          alpha = .5) +
     geom_sf(data = big_streets$osm_lines,
          inherit.aes = FALSE,
          color = "black",
          size = .5,
          alpha = .7) +
     geom_sf(data = bridleway$osm_lines,
          inherit.aes = FALSE,
          color = "brown",
          size = .4,
          alpha = .4)

PLOT <-  PLOT + geom_point(aes(x = 0.256, y = 52.024), shape = 19, col = "red", fill = "red", size = 10, alpha = 0.6) +
          geom_rect(aes(xmin = 0.2011256, ymin = 51.9815690, xmax = 0.2811256, ymax = 51.988), fill="white")


PLOT <- PLOT + coord_sf(xlim = c(0.2011256, 0.2811256),
           ylim = c(51.9815690, 52.0615690),
           expand = FALSE) +
     theme_void() +
     labs(title = "\nSAFFRON WALDEN", subtitle = "52.024°N / 0.256°E")

PLOT <- PLOT + theme(plot.title = element_text(size = 100, family = "lato", face="bold", hjust=0.5), plot.subtitle = element_text(family = "lato", size = 40, hjust=0.5, margin=margin(2, 0, 5, 0)))


# save it
ggsave(plot=PLOT, "Desktop/saffron_walden.pdf",height=1189,width=841,units="mm", useDingbats=FALSE)
