# PBMS basic app 24/2/2020

#fontawesome: crow, dove, kiwi-bird,earlybirds,feather
# in the future

# leafelet markers color: a bit tricky https://stackoverflow.com/questions/32940617/change-color-of-leaflet-marker

# good package: https://cran.r-project.org/web/packages/googleway/vignettes/googleway-vignette.html
# save a mapshot: https://r-spatial.github.io/mapview/reference/mapshot.html

## to fix 30/6/2020:
# 1. randomize 10 KM gridref
# 2. fix chem data specified by names
# 2b. fix region by species or year color (as.factor )
# 3. state tag plots for PCB
# 4. some PCB map?

packrat::on(project = '/data/PBMS/')
#shinyWidgets::shinyWidgetsGallery()

set.seed(100)

library(shiny)
library(leaflet)
library(readxl)
library(dplyr)
library(shinyWidgets)
library(shinycssloaders)
library(rgdal)
library(raster)
library(shinythemes)
library(stringr)
library(lubridate)
library(mapview)   # more advanced leaflet
library(plotly)
library(readr)
library(sf)
library(rnrfa) # for osg_parse
library(rnaturalearth) # dist to coast: https://dominicroye.github.io/en/2019/calculating-the-distance-to-the-sea-in-r/

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
  return(regnum)
}


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
  select(CARCASS_ID,ORIGINAL_SPECIES,COLLECTION_YEAR,LOCATION_10KM_GRIDREF,LOCATION_GRIDREF) %>% 
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

ui <- fillPage(
  #title = "UK Predatory Bird Monitoring Scheme",
  #theme = shinytheme("simplex"),
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
      sliderInput(inputId = 'years', label = 'Year range', 
                  min = min_yr, max = max_yr, value = c(1994,2002),animate = T),
      selectInput(inputId = 'species', label = 'Species', 
                  multiple =TRUE,
                  selectize = TRUE,
                  selected = c('Barn Owl','Heron','Sparrowhawk','Kestrel'), 
                  choices = species_choices),
      #textInput(inputId = 'id_filter',label = 'Filter by Carcass ID:', value = ''),
      shinyWidgets::searchInput(
        inputId = "id_filter", 
        label = "Filter by Carcass ID:", 
        placeholder = "This is a placeholder", 
        btnSearch = icon("search"), 
        btnReset = icon("remove"), 
        width = "100%"
      ),
      checkboxInput(inputId = 'overlay_built', 'Overlay built areas (England and Wales)'),
      p("Save plot"),
      radioButtons(inputId = 'marker_option', label = 'Marker option:', 
                  selected = 'circles', 
                  choices = c('pins','circles','PCB'), inline = T)
        ),
        tabPanel("Regions",
                 # selectInput(inputId = 'loc_cluster_var', label = 'Location clustering variables:', 
                 #             selected = c('Latitude','Longitude'), 
                 #             mulitple = TRUE,
                 #             choices = c('Latitude','Longitude', 'Altitude', 'Distance to Sea')) ,
                 radioButtons(inputId = 'region_option', label = 'Regionalization option:', 
                              selected = 'default', 
                              choices = c('default','UK Economic regions','clusters','species','year'), inline = T),
                 selectInput(inputId = 'cluster_vars', label = 'Clustering variables', 
                             multiple =TRUE,
                             selectize = TRUE,
                             selected = c("Long","Lat","elevation","dist_to_coast"), 
                             choices = c("Longitude" = "Long","Latitude" = "Lat","Elevation" = "elevation",
                                         "Distance to coast" = "dist_to_coast")),
                 plotlyOutput("region_boxplot")
                 ),
        tabPanel("Analysis",
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
      circle = TRUE, status = "danger", icon = icon("sliders"), inline = F, width = "450px",
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
      dplyr::select(-region) %>% 
      dplyr::mutate(region = case_when(input$region_option == "UK Economic regions" ~ EU_region,
                                       input$region_option == "clusters" ~ clusters,
                                       input$region_option == "species" ~ `Bird species`,
                                       input$region_option == "year" ~ `Year submitted to PBMS`,
                                       #a == 0 | a == 1 | a == 4 | a == 3 |  c == 4 ~ 3,
                                       TRUE ~ as.factor(1)))
  })
  
  #factpal <- colorFactor("RdYlBu", 1:2, n = 2)
  #factpal <- colorFactor(topo.colors(5), 1:5, n = 5)
  #factpal <- colorFactor("viridis", 1:5, n = 5)
  #factpal <- colorFactor(cbPalette[1:5], 1:5, n = 5)
  
  
  output$mymap <- renderLeaflet({                          
    print(head(PBMSdata2() ))
    print(input$region_option)
    nlevels <- length(unique(PBMSdata2()$region))
    if (input$region_option== 'year'| input$region_option== 'species'){
      factpal <- colorFactor("viridis",levels = sort(unique(PBMSdata2()$region),decreasing=TRUE)) # new factors not quite working
    } else {
      factpal <- colorFactor(cbPalette[1:nlevels], 1:nlevels, n = nlevels)
    }
    
    m = leaflet(data=PBMSdata2() ) %>% addTiles()
    if (input$marker_option == "circles") {
        m = m %>% addCircleMarkers(~Long,~Lat,
                 color = ~factpal(region),
                 popup = ~paste("<b><font color='#0281a4'>UK Predatory Bird Monitoring Scheme</font></b>",
                                "<br/><b>Year:</b>",`Year submitted to PBMS`, 
                                "<br/><b>PBMS ID:</b>",`CARCASS_ID`, 
                                "<br/><b>Species:</b>",`Bird species`,
                                "<br/><b>Common PCB conc.:</b>", PCB, 
                                "<p><font color='#c2c5cc'> &copy;" ,year(Sys.Date()),"UK Centre for Ecology and Hydrology </font></p>") ,
                 label = ~as.character(`Year submitted to PBMS`)) %>% 
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
                                                "<br/><b>PBMS ID:</b>",`CARCASS_ID`, 
                                                "<br/><b>Species:</b>",`Bird species`,
                                                "<br/><b>Common PCB conc.:</b>", PCB, 
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
    m
  })
  
  ####  location clustering ###
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
