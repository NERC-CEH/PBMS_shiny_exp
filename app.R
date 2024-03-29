# PBMS basic app 24/2/2020

# ShinyHelper pakcagehttps://cwthom94.shinyapps.io/shinyhelper-demo/

# try this sometime https://github.com/cwthom/earl2019-shiny

#fontawesome: crow, dove, kiwi-bird,earlybirds,feather
# in the future

# leafelet markers color: a bit tricky https://stackoverflow.com/questions/32940617/change-color-of-leaflet-marker

# image in leaflet popup markers: https://stackoverflow.com/questions/36433899/image-in-r-leaflet-marker-popups

# overlay image exact location (say old map or EM map)

# good package: https://cran.r-project.org/web/packages/googleway/vignettes/googleway-vignette.html
# save a mapshot: https://r-spatial.github.io/mapview/reference/mapshot.html

## to fix 3/2/2021:
# 1. randomize 10 KM gridref (skip)
# 2. fix chem data specified by names
# 2b. fix region by species or year color (as.factor )
# 3. economic regions not number but names
# 4. bird details: see Lee's comments: comment out PCB concs., add post mortem year, cause of death(NA= to be examined, need live OCDB connection, need to test this)
# 5. Finish mapsave button
# 6. Fix logo disappearing (https://github.com/r-spatial/leafem/issues/22) --> fix later, backup first
# 7. Fix randomnes for clustering (why change when seed is fixed, prob data change)
# 8. Use htmltools::htmlEscape in leaflet popup for security
# ? how to give a feel of where most birds are coming from?


# DONE. added help tool tip https://stackoverflow.com/questions/46648471/popover-tooltip-for-a-text-in-shiny-app-using-shinybs /MT: 20210219
# turn to mobile app with ShinyMobile?
# plot a ggplot occurence heat map? >> https://ourcodingclub.github.io/tutorials/seecc_1/
# geom_density and geom_contour >> doable in R but much easier with folium
# https://gis.stackexchange.com/questions/168886/r-how-to-build-heatmap-with-the-leaflet-package


packrat::on(project = '/data/PBMS/')
#shinyWidgets::shinyWidgetsGallery()

set.seed(100)

library(shiny)
library(leaflet)
#library(leafem) # imported by mapview
library(readxl)
library(dplyr)
library(shinyWidgets)
library(shinycssloaders)
library(rgdal)
library(raster)
library(shinythemes)
library(stringr)
library(lubridate)
library(mapview)   # more advanced leaflet, may need manually untar Phantom_JS, need copy bin/phantomjs to add to a location in $PATH, save leaflet map： https://stackoverflow.com/questions/44259716/how-to-save-a-leaflet-map-in-shiny
library(plotly)
library(readr)
library(sf)
library(rnrfa) # for osg_parse
library(rnaturalearth) # dist to coast: https://dominicroye.github.io/en/2019/calculating-the-distance-to-the-sea-in-r/
#library(htmltools) # important to use it for security on leaflet popup

Sys.getenv("PATH")
#mapshot(leaflet() %>% addTiles(),file='test.png')



chemcode_lookup = read_csv('/data/PBMS/PBMS_CHEM_CODE_DESC.csv')

lookup_dist_coast = function(PBMSdata,latlong){
  
  uk <- ne_countries(scale = 'medium', country = "united kingdom", returnclass = "sf")
  uk <- st_transform(uk, 3055) # transform to UTM
  #grid <- st_make_grid(uk, cellsize = 5000, what = "centers")
  #grid <- st_intersection(grid, uk)   
  uk_coastline <- st_cast(uk, "MULTILINESTRING")
  #dist <- st_distance(uk_coastline, grid)
  
  points = st_as_sf( PBMSdata %>% select(Long,Lat) , 
                     coords = c("Long", "Lat"), 
                     crs = latlong)  %>% 
    st_transform(3055)
  dist = st_distance(uk_coastline, points) %>% as.numeric
  return(dist)
}

lookup_EU = function(X,Y){
  coords <- c(X,Y)  # temporary override MT 191122
  proj <-  "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +towgs84=446.448,-125.157,542.06,0.1502,0.247,0.8421,-20.4894 +units=m +no_defs"
  coords_sp <- SpatialPoints(matrix(coords, nrow = length(coords)/2, ncol = 2), proj4string= CRS(proj))
  regnum <- regrast[cellFromXY(regrast , coords_sp)]
  
  EU_names <- data.frame(num=1:11, names=c("Wales","London" , "North West","Scotland","West Midlands","Yorkshire and the Humber",
                                           "South East", "East of England", "East Midlands", "North East", "South West" ))
  EU_names <- setNames(EU_names$names,EU_names$num)
  regnum <- EU_names[regnum]
  return(regnum)
}
# previewColors(colorFactor("Paired", domain = NULL), EU_names)

#distinct(chem,CHEM_CHEMICALCODE)

chem = read_csv('/data/PBMS/BOP_ANALYTICAL_SAMPLE.csv') %>% 
  filter(!is.na(CHEM_CONCENTRATION)) %>% 
  mutate(TISSUE = ifelse(TISSUE=='liver','Liver',TISSUE)) %>% 
  select(CARCASS_ID, TISSUE ,BATCH,CHEM_CHEMICALCODE ,CHEM_CONCENTRATION) %>% 
  mutate(temp_chem_name = ifelse(is.na(as.numeric(CHEM_CHEMICALCODE)),CHEM_CHEMICALCODE,NA)) %>% 
  mutate(CHEM_CHEMICALCODE = ifelse(is.na(as.numeric(CHEM_CHEMICALCODE)),
                                    case_when( temp_chem_name == 'A-HCH' ~ 1001,
                                               temp_chem_name == 'Brod' ~ 1004,
                                               temp_chem_name == 'Brodifacoum' ~ 1004,
                                               temp_chem_name == 'Brom' ~ 1005,
                                               temp_chem_name == 'Bromadiolone' ~ 1005,
                                               temp_chem_name == 'DDE' ~ 1010,
                                               temp_chem_name == 'DDT' ~ 1011,
                                               temp_chem_name == 'Dife' ~ 1012,
                                               temp_chem_name == 'Difenacoum' ~ 1012,
                                               temp_chem_name == 'Difenacoum-Corrected' ~ 1013,
                                               temp_chem_name == 'Floc' ~ 1016,
                                               temp_chem_name == 'Flocumafen' ~ 1016,
                                               temp_chem_name == 'Flocumafen-Corrected' ~ 1017,
                                               temp_chem_name == 'G-HCH (BHC)' ~ 1018,
                                               temp_chem_name == 'HCB' ~ 1019,
                                               temp_chem_name == 'HEOD' ~ 1020,
                                               temp_chem_name == 'Hepox' ~ 1021,
                                               temp_chem_name == 'Matched PCB' ~ 1028,
                                               temp_chem_name == 'mercury' ~ 1023,
                                               temp_chem_name == 'Mercury' ~ 1023,
                                               temp_chem_name == 'Sum Congener' ~ 1025,
                                               temp_chem_name == 'TDE' ~ 1026,
                                               temp_chem_name == 'Total PCBs.' ~ 1028,
                                               TRUE ~ NA_real_ ),
                                    CHEM_CHEMICALCODE))%>% 
  filter(!is.na(as.numeric(CHEM_CHEMICALCODE))) %>%           ## need to correct for those put as chem names
  left_join(chemcode_lookup %>% mutate(CHEM_CHEMICALCODE = as.character(CHEM_CHEMICALCODE))) %>% 
  filter(CHEM_CHEMICALCODE %in% c(128,138,153,179,180)) %>%    ## five common PCB congenors in CEH 2003 paper
  group_by(CARCASS_ID) %>% 
  summarise(PCB = sum(CHEM_CONCENTRATION ))



cc = read_csv('/data/PBMS/BOP_CARCASS_DETAILS.csv') %>% 
  mutate(PM_DATE = make_date(PM_YEAR,PM_MONTH,PM_DAY)) %>% 
  select(CARCASS_ID,ORIGINAL_SPECIES,COLLECTION_YEAR,LOCATION_10KM_GRIDREF,LOCATION_GRIDREF,PM_DATE,PM_YEAR) %>% 
  rename(`Bird species` = ORIGINAL_SPECIES, `Year submitted to PBMS` = COLLECTION_YEAR) %>% 
  left_join(chem)

GRIDREF_rand = function(GRIDREF_10KM){
  a = runif(1, min=0, max=10000) %>% floor() %>%  str_pad( 6, pad = "0")
  return(paste0(substr(GRIDREF_10KM, 1, 3) , substr(a, 2, 4), 
                substr(GRIDREF_10KM, 4, 4) , substr(a, 1, 1), substr(a, 5, 2) ))
}

cc = cc %>%   # randomize 10KM grid
  mutate(LOCATION_GRIDREF = if_else(str_length(LOCATION_GRIDREF) == 8 ,
                                   LOCATION_GRIDREF, 
                                   GRIDREF_rand(LOCATION_10KM_GRIDREF)  ))  

osg_parse2 <- function(grid_refs) { 
  # error handling and randomize if only 10km is avail.
  out <- tryCatch(
    {
      # if (nchar == 4) {
      #   #randomize
      #   grid_refs = grid_refs
      # } else{
      #   grid_refs = grid_refs
      # }
      osg_parse(grid_refs)
    },
    error=function(cond) {
      # message("Here's the original error message:")
      # message(cond)
      # message("")
      # print(grid_refs)
      # Choose a return value in case of error
      return(list(easting = NA, northing = NA))
    },
    warning=function(cond) {
      # message(paste("URL caused a warning"))
      # print(grid_refs)
      return(list(easting = NA, northing = NA))
    },
    finally={
      #message("Some other message at the end")
    }
  )
  
return(out)
}
osg_parse2('TL252100')
#cc2a = lapply(cc2[1:2000],FUN = osg_parse2)
#cc2a = lapply(cc$LOCATION_GRIDREF,FUN = osg_parse2)
XY_from_osg = sapply((cc) %>% select(LOCATION_GRIDREF) %>% pull(),FUN = osg_parse2) %>% t() %>% as.data.frame()

cc = cc %>% mutate(X = XY_from_osg$easting  %>% unlist() %>% as.vector(),
              Y = XY_from_osg$northing %>% unlist() %>% as.vector()) # IT WORKS!

uk <- ne_countries(scale = 'medium', country = "united kingdom", returnclass = "sf")
uk <- st_transform(uk, 3055) # transform to UTM
#grid <- st_make_grid(uk, cellsize = 5000, what = "centers")
#grid <- st_intersection(grid, uk)   
uk_coastline <- st_cast(uk, "MULTILINESTRING")
#dist <- st_distance(uk_coastline, grid)

pbms_logo <- "pbms-logo.png"

cbPalette <- c("#0072B2", "#E69F00",   "#009E73","#CC79A7",  "#D55E00","#F0E442", "#999999", "#56B4E9",
               "#0072B2","#56B4E9",  "#009E73","#CC79A7","#E69F00",    "#D55E00","#F0E442", "#999999",
               "#0072B2","#56B4E9",  "#009E73","#CC79A7","#E69F00",    "#D55E00","#F0E442", "#999999",
               "#0072B2","#56B4E9",  "#009E73","#CC79A7","#E69F00",    "#D55E00","#F0E442", "#999999")


# PBMSdata <- read_excel("Data for Michael Tso.xlsx") %>% 
#   mutate(`Bird species` = str_to_title(`Bird species`, locale = "en"))  # upper case 

PBMSdata = cc %>% filter(!is.na(X))

ukgrid <- "+init=epsg:27700"
latlong <- "+init=epsg:4326"
coords <- cbind(Easting = as.numeric(as.character(PBMSdata$X)),
                Northing = as.numeric(as.character(PBMSdata$Y)))
dat_SP <- SpatialPointsDataFrame(coords,
                                 data = PBMSdata %>% select(X,Y),
                                 proj4string = CRS(ukgrid))
dat_SP_OSGB <- spTransform(dat_SP, CRS(latlong))
PBMSdata$Long <- coordinates(dat_SP_OSGB)[, 1]
PBMSdata$Lat  <- coordinates(dat_SP_OSGB)[, 2]

# convert to lat long first
# lat/long to sea: https://www.doogal.co.uk/DistanceToSea.php
# elevation data from Google: http://more.stevemorse.org/latlonbatch2.html?direction=altitude

species_choices <- PBMSdata %>% distinct(`Bird species`) %>% pull()
min_yr = min(PBMSdata$`Year submitted to PBMS`,na.rm = T)
max_yr = max(PBMSdata$`Year submitted to PBMS`,na.rm = T)

### add colum of distance to coast and altitude (run 3rd party)
#a = read_csv('www/distance_to_coast.csv')


PBMSdata = PBMSdata %>% 
  mutate(region = as.factor(1),
         #dist_to_coast2 = as.numeric(a$`Distance KMs`), 
         dist_to_coast = lookup_dist_coast(PBMSdata,latlong),
         elevation = raster::extract(raster::getData('alt', country = "GB") , PBMSdata %>% select(Long,Lat), method = "bilinear") ) %>% 
  filter(!is.na(elevation))

### add column of EU economic regions
regrast <- raster("/data/PBMS/www/regionrast")
#levels(regrast)[[1]]$NAME <- str_replace(as.character(levels(regrast)[[1]]$NAME)," Euro Region","") # alter region names

PBMSdata = PBMSdata %>% dplyr::mutate(EU_region = as.factor(lookup_EU(X,Y))) %>% 
      filter(!is.na(EU_region))

find_clusters = function(a){  # automatically employ elbows method (last point to have reduction of misfit by >25%)
  # Nbclust: https://www.datanovia.com/en/lessons/determining-the-optimal-number-of-clusters-3-must-know-methods/ 
  #f_of_x <- splinefun(x,sin(pi/x)) 
  #f_of_x(x,deriv = 1) # take first derivative of function
  a2 = a[complete.cases(a),]
  clusters = rep(NA, length(a))
  k_candidate = 1:10 # no gaps
  misfit_k = rep(NA, length(k_candidate))
  i=1
  for (k in k_candidate) {
    misfit_k[i] = kmeans(scale(a2),k)$tot.withinss
    i = i+1
  }
  ratio = misfit_k[2:length(misfit_k)]/misfit_k[1:length(misfit_k)-1] 
  k_opt = which(ratio <0.75) %>% max() + 1 # add error handling, set default to 5
  clusters[complete.cases(a)] = kmeans(scale(a2),k_opt)$cluster
  return(as.factor(clusters))
}
#find_clusters(USArrests)

# 
# PBMSdata = left_join(PBMSdata,cbind(a,b) %>% select(), by = "Long") %>% rename(cluster = b)
# 
# leaflet(data=PBMSdata) %>% addTiles() %>%
#   addCircleMarkers(~Long,~Lat,
#                    #clusterOptions = markerClusterOptions(),
#                    color = ~factpal(as.factor(clusters)),
#                    popup = ~paste("<b><font color='#0281a4'>UK Predatory Bird Monitoring Scheme</font></b>",
#                                   "<br/><b>Year:</b>",`Year submitted to PBMS`,
#                                   "<br/><b>Species:</b>",`Bird species`,
#                                   "<p><font color='#c2c5cc'> &copy;" ,year(Sys.Date()),"UK Centre for Ecology and Hydrology </font></p>") ,
#                    label = ~as.character(`Year submitted to PBMS`)) %>%
#   addLegend('bottomright', pal = factpal, values = PBMSdata$cluster,
#             title = 'Bird locations<br>Region:',
#             opacity = 1)


######## UI #########

ui <- fillPage(
  #title = "UK Predatory Bird Monitoring Scheme",
  #theme = shinytheme("simplex"),
  #div(id="logo", style="position:fixed;top:0px;left:0px;z-index:11000;",
  #       "pbms-logo.png">'),
  tags$style(".fa-info-circle {color:#E87722}"), # change color of font awesome
  leafletOutput("mymap", width = "100%", height = "100%"), 
  absolutePanel(
    id = "controls", class = "panel panel-default", fixed = TRUE,
    draggable = TRUE, top = 55, left = "8%", #125 
    right = "auto", bottom = "auto",
    width = 0, height = 0,
    dropdownButton(
      inputId = "mydropdown",
      tags$h4("UK Predatory Bird Monitoring Scheme"),
      tabsetPanel(
        tabPanel("Home",
                 p(),
      p(icon('info'), 'The map shows PBMS bird submissions based on user input selection.'),
      p(icon('info'), 'Hover on the ',icon('info-circle'), 'button next to input boxes for help text.'),
      p(icon('info'), 'Open/minimize the control panel by clicking the red button.'),
      p(icon('info'), 'Move the control panel by click and drag on its white space when it is open.'),
      p(icon('info'), 'Hover and click on bird for individual bird details.'),
      
      
      tags$div(title="Select the year range. Use the play button to animate changes every N years.",
          sliderInput(inputId = 'years', label = tagList("Year range:",icon('info-circle')),
                  min = min_yr, max = max_yr, value = c(1994,2002),animate = T,sep="")),
      tags$div(title="Select the species to show on map. Click on box for options. Backspace to delete choice.",
          selectInput(inputId = 'species', label = tagList("Species:",icon('info-circle')), 
                  multiple =TRUE,
                  selectize = TRUE,
                  selected = c('Barn Owl','Sparrowhawk','Kestrel'), 
                  choices = species_choices)),
      #textInput(inputId = 'id_filter',label = 'Filter by Carcass ID:', value = ''),
      shinyWidgets::searchInput(
        inputId = "id_filter", 
        label = "Filter by Carcass ID:", 
        placeholder = "Enter a single PBMS Carcass ID", 
        btnSearch = icon("search"), 
        btnReset = icon("remove"), 
        width = "70%"
      ),
      tags$div(title="It will take a few seconds.",
               checkboxInput(inputId = 'overlay_built', tagList("Overlay built areas (England and Wales)",icon('info-circle')))),
      p(),
      tags$div(title="You can download a static image of the interactive map.",
        downloadButton("dl",label = "Download map")),
      
      tags$div(title='The type of marker to represesnt each bird on map.',
               radioButtons(inputId = 'marker_option', label = tagList('Marker option:',icon('info-circle')), 
                  selected = 'circles', 
                  choices = c('pins','circles','PCB'), inline = T))
     
        ),
      div(class="flexcontainer", 
          
          # action button
          actionButton(inputId="pdf", label="Go to User Manual", class="btn-success" , 
                       onclick = "window.open('PBMS explorer App User Manual.pdf')")
      )
      ,
      
        tabPanel("Regions",
                 # selectInput(inputId = 'loc_cluster_var', label = 'Location clustering variables:', 
                 #             selected = c('Latitude','Longitude'), 
                 #             mulitple = TRUE,
                 #             choices = c('Latitude','Longitude', 'Altitude', 'Distance to Sea')) ,
                 p(),
                 p(icon('info'), 'Seclect options to color circle markers and to group variables for analysis.'),
                 tags$div(title="This controls the groupings in the circle markers and analysis tab.",
                     radioButtons(inputId = 'region_option', label = tagList('Regionalization option:',icon('info-circle')), 
                              selected = 'default', 
                              choices = c('default','UK Economic regions','clusters','species','year'), inline = T)),
                 tags$div(title="The variables used for clustering if clusters is selected above. Click on box for options. Backspace to delete choice.",
                 selectInput(inputId = 'cluster_vars', label = tagList('Clustering variables',icon('info-circle')), 
                             multiple =TRUE,
                             selectize = TRUE,
                             selected = c("Long","Lat","elevation","dist_to_coast"), 
                             choices = c("Longitude" = "Long","Latitude" = "Lat","Elevation" = "elevation",
                                         "Distance to coast" = "dist_to_coast"))),
                 plotlyOutput("region_boxplot")
                 ),
        tabPanel("Analysis",
                 p(),
                 p(icon('info'), 'These plots show the changes in numbers of PBMS records over the years, grouped by the option selected in the "Regions" tab.'),
                 # selectInput(inputId = 'species2', label = 'Species', 
                 #             selectize = TRUE,
                 #             selected = species_choices, 
                 #             choices = species_choices),
                 # sliderInput("a_val", "prediction interval opaqueness (%)", 50, min=0, max=100, step = 1),
                 # checkboxInput("show_pred_int", label="show prediction intervals",value = TRUE),
                 # checkboxInput("remove_negative", label="crop negative prediction intervals",value = TRUE),
                 plotlyOutput("prediction") %>% withSpinner(color="#0dc5c1"),
                 plotlyOutput("prediction_PCB", height = "600px") %>% withSpinner(color="#0dc5c1")
                 )
      ),
      circle = TRUE, status = "danger", icon = icon("sliders"), inline = F, width = "550px",
      tooltip = tooltipOptions(title = "Click to see inputs !")
    )    
  ),
  p()
)

server <- function(input, output, session) {
  
  PBMSdata_filter<- reactive({
    print(input$id_filter)
    PBMSdata %>% 
      dplyr::filter(`Bird species` %in% input$species) %>% 
      dplyr::filter(between(`Year submitted to PBMS`, input$years[1], input$years[2])) 
  })
  
  PBMSdata_highlight<- reactive({
     PBMSdata_filter () %>% 
        {if (!input$id_filter=="") dplyr::filter(.,  CARCASS_ID == as.numeric(input$id_filter)) else .}
  })  
  
  PBMSdata2<- reactive({
    PBMSdata_filter() %>% 
      dplyr::mutate(clusters = find_clusters(PBMSdata_filter() %>% select(!!input$cluster_vars))) %>%  ### here
      dplyr::mutate(`Bird species` = as.factor(`Bird species`),
                    `Year submitted to PBMS` = as.factor(`Year submitted to PBMS`)) %>% 
      #dplyr::select(-region) %>% 
      #dplyr::mutate(region=as.character(region)) %>%  # 210204: need to un-factor then re-factor to avoid NA
      # dplyr::mutate(region = case_when(input$region_option == "UK Economic regions" ~ EU_region,  # case_when doesn't work for non-integers
      #                                  input$region_option == "clusters" ~ clusters,
      #                                  input$region_option == "species" ~ `Bird species`,       
      #                                  input$region_option == "year" ~ `Year submitted to PBMS`,
      #                                  #a == 0 | a == 1 | a == 4 | a == 3 |  c == 4 ~ 3,
      #                                  TRUE ~ as.factor(1))) %>% 
      #mutate(region = `Bird species`)
    
    #%>%  
      {if (input$region_option == "clusters") dplyr::mutate(., region= clusters) else .} %>% 
      {if (input$region_option == "UK Economic regions") dplyr::mutate(., region= EU_region) else .} %>% 
      {if (input$region_option == "year") dplyr::mutate(., region= `Year submitted to PBMS`) else .} %>% 
    {if (input$region_option == "species") dplyr::mutate(., region= `Bird species`) else .} 
    
    #  mutate(region = ifelse(input$region_option == "species", `Bird species`, region)) %>% 
     # mutate(region = ifelse(input$region_option == "year", `Year submitted to PBMS`, region))
  })
  
  #factpal <- colorFactor("RdYlBu", 1:2, n = 2)
  #factpal <- colorFactor(topo.colors(5), 1:5, n = 5)
  #factpal <- colorFactor("viridis", 1:5, n = 5)
  #factpal <- colorFactor(cbPalette[1:5], 1:5, n = 5)
  
  factpal <- colorFactor("Paired",domain=NULL)
  
  
  
  map <- reactiveValues(dat = 0)
  
  # leaflet map, do this way to can pass reactive output to mapshot: https://stackoverflow.com/questions/44259716/how-to-save-a-leaflet-map-in-shiny                 
  # Create foundational leaflet map
  # and store it as a reactive expression
  foundational.map <- reactive({
    print(head(PBMSdata2() ))
    print(input$region_option)
    nlevels <- length(unique(PBMSdata2()$region))
    # if (input$region_option== 'year'| input$region_option== 'species'){
    #   factpal <- colorFactor("viridis",levels = sort(unique(PBMSdata2()$region),decreasing=TRUE)) # new factors not quite working
    # } else {
    #   factpal <- colorFactor("Paired",domain=NULL)
    # }
    
    m = leaflet(data=PBMSdata2() ) %>% # addTiles() %>% 
      addProviderTiles(providers$Stamen.TonerLite,
                       options = providerTileOptions(noWrap = TRUE)) %>% 
          leafem::addLogo(pbms_logo, url = "https://pbms.ceh.ac.uk",
                          position = "bottomleft",
                          offset.x = 10, offset.y = 10,
                          width = 320, height = 214)
    if (input$marker_option == "circles") {
        m = m %>% addCircleMarkers(~Long,~Lat,
                 color = ~factpal(region),
                 popup = ~(paste("<b><font color='#0281a4'>UK Predatory Bird Monitoring Scheme</font></b>",
                                "<br/><b>Submission year:</b>",`Year submitted to PBMS`,
                              #  "<br/><b>Post mortem date:</b>",as.character(`PM_DATE`),
                                "<br/><b>PBMS ID:</b>",`CARCASS_ID`, 
                                "<br/><b>Species:</b>",`Bird species`,
                              #  "<br/><b>Common PCB conc.:</b>", PCB, 
                                "<p><font color='#c2c5cc'> &copy;" ,year(Sys.Date()),"UK Centre for Ecology and Hydrology </font></p>")) ,
                 label = ~(as.character(`Year submitted to PBMS`))) %>% 
        addLegend('bottomright', pal = factpal, values = PBMSdata2()$region,
                title = paste0('Bird locations<br>',input$region_option,':'),
                opacity = 1)
    } else if (input$marker_option == "PCB") {
      pal <- colorNumeric(
        palette = "YlGnBu",
        domain = PBMSdata2()$PCB
      )
      m = m %>% addCircleMarkers(~Long,~Lat,
                                 color = ~pal(PCB),
                                 popup = ~paste("<b><font color='#0281a4'>UK Predatory Bird Monitoring Scheme</font></b>",
                                                "<br/><b>Year:</b>",`Year submitted to PBMS`, 
                                            #  "<br/><b>Post mortem date:</b>",as.character(`PM_DATE`),
                                                "<br/><b>PBMS ID:</b>",`CARCASS_ID`, 
                                                "<br/><b>Species:</b>",`Bird species`,
                                            #    "<br/><b>Common PCB conc.:</b>", PCB, 
                                                "<p><font color='#c2c5cc'> &copy;" ,year(Sys.Date()),"UK Centre for Ecology and Hydrology </font></p>") ,
                                 #radius = ~ifelse(type == "ship", 6, 10),
                                 radius = ~(PCB)^2 ,
                                 label = ~as.character(`Year submitted to PBMS`)) %>% 
        addLegend('bottomright', pal = pal, values = ~PCB,
                  title = paste0('PCB conc.<br>','(mg/kg) :'),
                  opacity = 1)
    } else {
      m = m %>% addMarkers(~Long,~Lat,
                       clusterOptions = markerClusterOptions(),
                       #color = ~factpal(region),
                       popup = ~paste("<b><font color='#0281a4'>UK Predatory Bird Monitoring Scheme</font></b>",
                                      "<br/><b>Year:</b>",`Year submitted to PBMS`, 
                                      "<br/><b>PBMS ID:</b>",`CARCASS_ID`, 
                                      "<br/><b>Species:</b>",`Bird species`,
                                      "<br/><b>Common PCB conc.:</b>", PCB, 
                                      "<p><font color='#c2c5cc'> &copy;" ,year(Sys.Date()),"UK Centre for Ecology and Hydrology </font></p>") ,
                       label = ~as.character(`Year submitted to PBMS`)) 
    }
    
    if (input$overlay_built == 1) {
      ukbuilt <- rgdal::readOGR("www/uk_built_areas_2011/63efb944-5faa-48ac-843b-1b8bd0019bac202044-1-1pvcf1x.xo1ch.shp")
      
      m = m %>% addPolygons(data=ukbuilt, color = "#444444", weight = 1, smoothFactor = 0.5,
                            opacity = 1.0, fillOpacity = 0.5)
    }
    
    
    #leaflet(df.20) %>% addTiles() %>%
    #  addAwesomeMarkers(~long, ~lat, icon=icons, label=~as.character(mag))
    if (input$id_filter != '') {
      m = m %>% addAwesomeMarkers(data=PBMSdata_highlight(), ~Long,~Lat, icon = awesomeIcons(icon = 'ios-close',
                                                                                              iconColor = 'black',
                                                                                              library = 'ion',
                                                                                              markerColor = "pink"),
                                  popup = ~paste("<b><font color='#0281a4'>UK Predatory Bird Monitoring Scheme</font></b>",
                                                 "<br/><b>Year:</b>",`Year submitted to PBMS`, 
                                                 "<br/><b>PBMS ID:</b>",`CARCASS_ID`, 
                                                 "<br/><b>Species:</b>",`Bird species`,
                                                 "<br/><b>Common PCB conc.:</b>", PCB, 
                                                 "<p><font color='#c2c5cc'> &copy;" ,year(Sys.Date()),"UK Centre for Ecology and Hydrology </font></p>") ,
                                  label = ~as.character(`Year submitted to PBMS`))
    }
    map$dat <- m
      
    m 
  })
  # observe({
  #   leafletProxy("mymap") %>%
  #     leafem::addLogo(pbms_logo, url = "https://pbms.ceh.ac.uk",
  #               position = "bottomleft",
  #               offset.x = 10, offset.y = 10,
  #               width = 320, height = 214)
  # })
  
  output$mymap <- renderLeaflet({  
    withProgress(message = 'Loading plot', {
    foundational.map()
    })
  })
  
  ######## Map shot stuff #########
  # store the current user-created version of the Leaflet map for download in 
  # a reactive expression
  user.created.map <- reactive({
    
    # call the foundational Leaflet map
    foundational.map() %>%
      
      # store the view based on UI
      setView( lng = input$map_center$lng
               ,  lat = input$map_center$lat
               , zoom = input$map_zoom
      )
    
  }) # end of creating user.created.map()
  # create the output file name
  # and specify how the download button will take
  # a screenshot - using the mapview::mapshot() function
  # and save as a PDF
  # output$dl <- downloadHandler(
  #   filename = paste0( Sys.Date()
  #                      , "_customLeafletmap"
  #                      , ".png"
  #   )
  #   
  #   , content = function(file) {
  #     mapshot( x = user.created.map()
  #              , file = file
  #              , cliprect = "viewport" # the clipping rectangle matches the height & width from the viewing port
  #              , selfcontained = FALSE # when this was not specified, the function for produced a PDF of two pages: one of the leaflet map, the other a blank page.
  #     )
  #   } # end of content() function
  # ) # end of downloadHandler() function

  output$dl <- downloadHandler(
    filename = "map.png",
    
    content = function(file) {
      mapshot(map$dat, file = file)
    }
  )
    
  #### location clustering ###
  clusterdata = reactive({       # perform clustering
    
    cluster_data = tidydata() %>%
      dplyr::select(c("DATE",input$MA_choices )) #("DATE","DRYTMP","SOLAR","WSPEED","RAIN")
    cluster_data = cluster_data[ complete.cases(cluster_data) ,]
    
    clusters = kmeans( scale( cluster_data %>% dplyr::select(-DATE) ), input$n_cluster)  ##
    cluster_data = cluster_data %>%
      dplyr::mutate(states = as.factor(clusters$cluster))
    
    centers = clusters$centers
    data = tidydata()%>%
      dplyr::left_join(cluster_data %>% dplyr::select(DATE,states),by="DATE")
    

    
    out = list(data = data,centers = centers , stats = stats, lag_stat=lag_stat )
    return( out  )
  })
  
  
  pred_data_aa <- reactive({ 
    aa = PBMSdata2() %>%
      dplyr::group_by(`Year submitted to PBMS`,region) %>% 
      dplyr::summarise(VALUE = n()) #%>% 
      # dplyr::mutate(counts.mean = mean(VALUE), 
      #               counts.sd = sd(VALUE),
      #               counts.max = mean(VALUE)+as.numeric(input$n_sd)*sd(VALUE),
      #               counts.min = mean(VALUE)-as.numeric(input$n_sd)*sd(VALUE))
    # if(input$remove_negative == TRUE) {
    #   aa = aa %>% dplyr::mutate(counts.min = replace(counts.min,which(counts.min <0),0)) # input$remove_negative      
    # }
    return(aa)
  })
  
  output$prediction <- renderPlotly({
    aa = pred_data_aa() %>% tidyr::drop_na()
    p1 = ggplot( aa ,aes(x = `Year submitted to PBMS`, y = VALUE,colour=region, group = region) )
    
    # if(input$show_pred_int == TRUE) {
    #   p1 = p1 + geom_ribbon(aes(x = DATE, 
    #                             ymin=counts.min,
    #                             ymax=counts.max
    #   ),fill="steelblue2",color="steelblue2",alpha=input$a_val/100)
    # }
    p1 = p1 + geom_point() + geom_line() +
      scale_colour_manual(values=cbPalette)+
      #scale_x_date() +
      ggtitle(str_c("Number of sumbissions")) +
      ylab("Counts") + xlab('Year') +
      theme(legend.position = "s",
            axis.title.y = element_text(size = 10),
            plot.title = element_text(lineheight = 0.8, face = "bold")) 
    ggplotly(p1)  %>% 
      plotly::config(toImageButtonOptions = list(format="png"))
  })
  
  output$prediction_PCB <- renderPlotly({
    aa = PBMSdata2() 
    p1 = ggplot( aa ,aes(x = `Year submitted to PBMS`, y = PCB, fill = region) ) +
      geom_boxplot() +  facet_wrap(~region,ncol = 1) + 
      geom_hline(yintercept =  (aa %>% pull(PCB) %>% mean(na.rm=T)) , size=1.5, color="red")
    ggplotly(p1)  
  })
  
  ##### input controls ####
  
  # observeEvent(input$pdf, {
  #   # Absolute path to a pdf, use file.path()
  #   file.show("www/PBMS explorer App User Manual.docx")
  # })
  # 
  observeEvent(list(), {   ### show dropdown content during startup
    toggleDropdownButton(inputId = "mydropdown")
  }, ignoreInit = FALSE)
  
}

shinyApp(ui, server)
# 
# ukgrid <- "+init=epsg:27700"
# latlong <- "+init=epsg:4326"
# 
# coords <- cbind(Easting = as.numeric(as.character(PBMSdata$X)),
#                 Northing = as.numeric(as.character(PBMSdata$Y)))
# 
# dat_SP <- SpatialPointsDataFrame(coords,
#                                  data = PBMSdata,
#                                  proj4string = CRS(ukgrid))
# 
# dat_SP_OSGB <- spTransform(dat_SP, CRS(latlong))
# PBMSdata$Long <- coordinates(dat_SP_OSGB)[, 1]
# PBMSdata$Lat  <- coordinates(dat_SP_OSGB)[, 2]
# m = leaflet() %>% addTiles() %>%
#       addMarkers(data = PBMSdata, ~Long,~Lat,
#              #clusterOptions = markerClusterOptions(),
#              popup = ~paste("<b><font color='#0281a4'>UK Predatory Bird Monitoring Scheme</font></b>",
#                             "<br/><b>Year:</b>",`Year submitted to PBMS`, 
#                             "<br/><b>Species:</b>",`Bird species`,
#                             "<p><font color='#c2c5cc'> &copy;" ,year(Sys.Date()),"UK Centre for Ecology and Hydrology </font></p>") ,
#              label = ~as.character(`Year submitted to PBMS`)) 
# if (a == 1) {
#   m = m %>% addPolygons(data=ukbuilt, color = "#444444", weight = 1, smoothFactor = 0.5,
#                         opacity = 1.0, fillOpacity = 0.5)
# }
# m
# 
