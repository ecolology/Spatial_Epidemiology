
library(readr) 
library(tidyverse) # cleaning, wrangling 
library(sf) # spatial manipulation
library(tmap) # map visualisations
library(leaflet) # leaflet maps 
library(terra) # spatial data analysis with raster  
#install.packages("remotes") 
#remotes::install_github("wfmackey/absmapsdata") # tons of data in package!
# remotes::install_github("runapp-aus/strayr") # loads the data with strayr::read_absmap("sa42021")
library(strayr)
library(geodata)


# Open data ------------------------------------------------------
library(ozmaps)
aus <- ozmaps::ozmap('abs_ced') # Map of Australia  
head(aus)

library(Census2016)
census <- Census2016::Census2016_wide_by_SA2_year %>% 
  filter(year == 2016)
head(census)

library(absmapsdata)
# absmapsdata::absmapsdata_file_list
# aus2016 <- absmapsdata::sa22016
aus2016 <- strayr::read_absmap(area = "sa2", year = 2016)
head(aus2016)

ft <- readr::read_csv("raw_data/Farm_transparency.csv") %>%
  janitor::clean_names()
head(ft)


library(skimr)
skim(ft)

ggplot(ft, aes(x=lng, y=lat, color= type))+ 
  geom_point()+  
  theme_bw()

# this is a tibble, so if we use tmap, this is not going to work
# tmap::tm_shape(ft) +
#   tm_dots()

# Now we transfom the tibble to a spatial data using sf
ft_s <-  ft %>%
  sf::st_as_sf(coords = c('lng', 'lat'), crs = 4326) # WGS 84

st_crs(ft_s) #EPSG:4326 

ft_s <- st_transform(ft_s, 4283) #GDA94 
st_crs(ft_s) #EPSG:4326 


tm_shape(ft_s) + 
  tm_dots(col = 'type', size = 0.1, title = "farms in Australia",  palette = "RdBu" )


tmap_mode('view')  # use 'plot' to turn off the interactive view with Leaflet!

tm_shape(ft_s)+  
  tm_dots(col = 'type', size = 0.1, 
          title = "Farms of Australia", 
          palette = "RdBu" )

tmap_mode('plot')  # turn off the interactive view of tmap

# Now lets focus on QLD
names(aus2016)
unique(aus2016$sa2_code_2016)

qld <- aus2016 %>%
  filter(state_name_2016=="Queensland")

qld
plot(qld)

# we do not need all of those variables, so we will plot only SA2_code_2016
qld %>% 
  dplyr::select(sa2_code_2016) %>%
  plot()

# tm_shape(qld)+
#   tm_polygons()


qld %>%
  ggplot() +
  geom_sf(aes(geometry = geometry, fill = sa2_code_2016),
          show.legend = FALSE)

# the shape contains empty units
# Check for empty units
st_is_empty(qld) # The shape qld contains empty units. 

qld <- qld %>% filter (!st_is_empty(.))

# tm_shape(qld)+
#   tm_polygons()
# 
# tm_shape(qld) + 
#   tm_polygons(col = "grey10", alpha = 0.1)+
#   tm_shape(ft_s)+ 
#   tm_dots(col = 'type',      
#           size = 0.1,  
#           title = "Farms in QLD",        
#           palette = "RdBu") +  
#   tm_layout(legend.position = c('right', 'top'),   
#             legend.outside = TRUE)+ 
#   tm_compass(position = c('right', 'top'))+ 
#   tm_grid(lines = F)+   
#   tm_scale_bar(text.size = 1, position = c('left', 'bottom'))


# Now lets focus on only piggeries in QLD and calculate how many piggeries we have
# per 100,000 people
names(ft_s)
unique(ft_s$type)

farms <- ft_s %>% 
  filter(type %in% c('Pigs', 'Dairy', 'Broiler (Meat) Chickens',
                                   'Cattle (Beef)', 'Horse Racing') &
        state=='QLD')
head(farms)

# add the census to a spatial layer
# census + qld

names(qld)
names(census)

census_qld <- qld %>% 
  dplyr::select(sa2_code_2016) %>%
  mutate(sa2_code_2016 = as.integer(sa2_code_2016)) %>%
  left_join(census, by=c('sa2_code_2016'='sa2_code'))
head(census_qld)


# now we need the number of piggeries per SA2
#change the CRS
farms_sa2 <- farms %>% 
  st_transform(st_crs(census_qld)) %>% #GDA94 
  st_join(census_qld) 
 
  
census_qld2 <- census_qld %>% 
  mutate(farms = lengths(st_intersects(census_qld, farms_sa2)),
         incidence = farms/persons) # add the new column of piggeries incidence

# map the incidence
census_qld2 %>%
  ggplot() +
  geom_sf(aes(fill = farms)) +
  scale_fill_viridis_c() +
  geom_sf(data = farms_sa2, color = "grey", size=0.2, alpha = 0.2)



# tm_shape(census_qld) + 
#   tm_polygons() +
#   tm_bubbles(size = "farms")


# tmap_style("natural") #natural, cobalt

# tm_shape(census_qld) +
#   tm_polygons("farms", title = "Farms", style = "cont") +
#   tm_layout(legend.outside = TRUE) +
#   tm_scale_bar(position = c("left", "bottom")) + # add scale bar to the top right
#   tm_compass(type = "arrow",
#              position = c("right", "top"))


# -------------------------------------------------------------------------
#Leaflet

# color palette 
pal <- 
  colorBin(
    palette = "YlOrRd",
    domain = census_qld$farms) #change to your data here

# pop up message
labels <- 
  sprintf(
    "<strong>%s</strong><br/>%g",
  census_qld$sa2_name,  census_qld$farms) %>% #change to your data here
  lapply(htmltools::HTML)




library(shiny)
library(leaflet)
shinyApp(
  ui <- navbarPage("Leaflet", id="nav", 
                   # a tab for the map 
                   tabPanel(
                     "Interactive map",
                     shinycssloaders::withSpinner(leafletOutput(
                       outputId = "mymap", 
                       width = "900px", 
                       height = "500px"))),
                   # A tab to explore the data in table format
                   tabPanel("Explore the data",
                            DT::dataTableOutput("table"))
  ),
  
  server <- function(input, output) {
    
    # map panel 
    output$mymap <- renderLeaflet({
      
      # passing the shp df to leaflet
      leaflet(census_qld) %>% #change to your data here
        # zooming in on Brisbane
        setView(153.0260, -27.4705, 8) %>% # long/lat
        # adding tiles, without labels to minimize clutter
        addProviderTiles("CartoDB.PositronNoLabels") %>%
        # parameters for the polygons
        addPolygons(
          fillColor = ~pal(farms), 
          weight = 1,
          opacity = 1,
          color = "white",
          fillOpacity = 0.7,
          highlight = highlightOptions(
            weight = 2,
            color = "#666",
            fillOpacity = 0.7,
            bringToFront = TRUE),
          label = labels,
          labelOptions = labelOptions(
            style = list("font-weight" = "normal"),
            textsize = "15px",
            direction = "auto")) %>%
        # legend
        addLegend(pal = pal,
                  values = census_qld$farms, #change to your data here
                  position = "bottomright",
                  title = "Farms", #change to your data here
                  opacity = 0.8,
                  na.label = "No data")
    })
    
    # data panel
    output$table <- DT::renderDataTable({
      DT::datatable(pop_shp %>% st_drop_geometry(), rownames = F,  filter = 'top', #change to your data here
                    extensions = c('Buttons', 'FixedHeader', 'Scroller'),
                    options = list(pageLength = 15, lengthChange = F,
                                   fixedHeader = TRUE,
                                   dom = 'lfBrtip',
                                   list('copy', 'print', list(
                                     extend = 'collection',
                                     buttons = c('csv', 'excel', 'pdf'),
                                     text = 'Download'
                                   ))
                    ))
    })
    
  }
  
  ,
  
  options = list(height = 700)
)



# -------------------------------------------------------------------------

library(terra) 
library(geodata)


australia <- geodata::worldclim_country("Australia", var="tmin", path=tempdir())


#extract temperature
tmin_qld <- terra::crop(australia, qld)

tm_shape(tmin_qld)+
  tm_raster()


points <- ft %>% select(lng, lat,type, state) %>% 
  filter(state =='QLD' & 
           type %in% c('Pigs', 'Dairy', 'Broiler (Meat) Chickens',
                   'Cattle (Beef)', 'Horse Racing')) %>% 
  sf::st_as_sf(coords = c('lng', 'lat'), crs = 4326) # WGS 84
  
points$tmin <- terra::extract(tmin_qld, points)

###
# Load the ggplot2 package
library(ggplot2)

points$type <- as.factor(points$type)

ggplot()+
  geom_boxplot(aes(x=type, y=tmin$AUS_wc2.1_30s_tmin_1), data=points)


# Perform linear regression
model1 <- lm(points$tmin$AUS_wc2.1_30s_tmin_1 ~  points$type)
summary(model1)




