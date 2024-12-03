# COSST shiny application

# Just click in "Run App"

# Load packages
library(shiny)
library(tidyverse)
library(raster)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgbif)
library(geodata)
library(CoordinateCleaner)
library(spThin)
library(usdm)
library(dismo)
library(leaflet)
library(viridis)
library(viridisLite)
library(enmSdmX)
library(sf)
library(igraph)

# Define the UI
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .shiny-notification {
        position: fixed;
        top: 50%;
        left: 50%;
        transform: translate(-50%, -50%);
        max-width: 300px;
        width: auto;
      }
    "))
  ),
  titlePanel("Climate-Oriented Seed Sourcing Tool (COSST)"),
  sidebarLayout(
    sidebarPanel(
      width = 3,  # Make the sidebar narrower (3 out of 12 columns)
      textInput("species", "Species name", "Caryocar brasiliense Cambess."),
      numericInput("lat", "Latitude", value = -14.350088458305477),
      numericInput("long", "Longitude", value = -48.400013676905374),
      selectInput("provenance", "Provenance strategy", choices = c("Composite" = 'compo', 
                                                                   "Climate-adjusted" = 'clima.adjus', 
                                                                   "Predictive" = 'predi')),
      conditionalPanel(
        condition = "input.provenance == 'clima.adjus' || input.provenance == 'predi'",
        selectInput("ssp", "Future scenario", choices = c("Optimistic (SSP1)" = 126, "Pessimistic (SSP5)" = 585)),
        selectInput("time", "Time period", choices = c("2021-2040", "2041-2060", "2061-2080", "2081-2100"))
      ),
      actionButton("update", "Go"),
      sliderInput("opacity", "Map opacity", min = 0, max = 1, value = 0.8, step = 0.1),
      
      # CSV File Input
      fileInput("file1", "Upload CSV with seed sourcing sites coordinates", accept = ".csv"),
      
      p("This tool enables users to visualize priority areas for seed collection based on composite, climate-adjusted, and predictive strategies. Colored areas represent the species' geographical range, with warmer colors indicating higher priority for the selected provenance strategy. Please upload a CSV file containing coordinates for the sourcing sites (long, lat) to extract the COSST index at those locations.")
      ),
    mainPanel(
      leafletOutput("cosstMap", width = "100%", height = "800px"),  # Increased map height
      tableOutput("cosst_table")  # Updated table output ID
    )
  ),
)

# Define the server logic
server <- function(input, output) {
  
  observeEvent(input$update, {
    species <- input$species
    site <- data.frame(lat = input$lat, long = input$long)
    provenance <- input$provenance
    ssp <- input$ssp
    time <- input$time
    
    # Initialize Progress Bar
    withProgress(message = "Processing COSST Application...", value = 0, {
      
      # STEP 2: DOWNLOAD AND PROCESS THE BASE MAP, OCCURRENCE AND CLIMATIC DATA
      incProgress(0.15, detail = "Downloading and processing base map and occurrence data...")
      
      {
        
        # Base map
        
        # Make it spatial
        coordinates(site) = ~ long + lat
        
        # Specify the SDM training extent
        exten.train = as(ne_countries(), "Spatial") # World
        
        # Specify the SDM projection extent
        exten.proje = exten.train[exten.train$sovereignt == raster::intersect(site, exten.train)@data[["sovereignt"]], ] # Country of the site
        
        # Occurrence data
        
        # We provided a code to retrieve GBIF data directly from R, though note that for the paper's examples we have combined data from GBIF and SpeciesLink
        
        # Download GBIF data
        gbif = occ_data(scientificName = species,
                        limit = 10000)
        
        # Filter occurrence coordinates (lat, long)
        occur = data.frame(species = species,
                           lat = gbif[["data"]][["decimalLatitude"]],
                           long = gbif[["data"]][["decimalLongitude"]])
        
        # Remove NAs
        occur = drop_na(occur)
        
        # Remove coordinates off planetary limits
        occur = occur[occur$lat > -90 &
                        occur$lat < 90 &
                        occur$long > -180 &
                        occur$long < 180,]
        
        # Clean occurrences with CoodinateCleaner package: https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13152
        occur = clean_coordinates(x = occur,
                                  lon = "long",
                                  lat = "lat",
                                  species = "species",
                                  tests =  c("zero", "capitals", "centroids", "duplicates", "equal",
                                             "gbif", "institutions", "seas"),
                                  verbose = F)
        
        occur = occur %>% # Filter valid occurrences
          dplyr::filter(.summary == T) %>%
          dplyr::select(species, long, lat)
        
        # Spatial thinning
        occur = thin(loc.data = occur, # Dataset
                     verbose = F,
                     long.col = "long", # Longitude column
                     lat.col = "lat", # Latitude column
                     spec.col = "species", # Species column
                     thin.par = 5, # Removing occurrences closer than 5 km
                     reps = 1,
                     locs.thinned.list.return = T,
                     write.log.file = F,
                     write.files = F)
        
        # Convert to data frame
        occur = data.frame(occur)
        
        # Climatic data
        # We provided a code to retrieve WorldClim data directly from R, though note that for the paper's examples we have used data from CHELSA climatologies instead
        
        # Download WorldClim data
        bioclim = worldclim_global(var = "bio", res = 10, path = tempdir()) # Resolution of 10 minutes of a degree
        
        # Rename it
        names(bioclim) = c("bio01", "bio02", "bio03", "bio04", "bio05", "bio06", "bio07",
                           "bio08", "bio09", "bio10", "bio11", "bio12", "bio13", "bio14",
                           "bio15", "bio16", "bio17", "bio18", "bio19")
        
        # Remove BIO18 and BIO19 due to discontinuity issues, see https://doi.org/10.1111/aec.13234
        bioclim = bioclim[[1:17]]
        
        # Crop to the training extent
        biocl.train = raster::crop(x = bioclim,
                                   y = exten.train)
        
        # Compute the variance inflation factor (VIF)
        vif = vifstep(biocl.train)
        
        # Exclude variables with high VIF
        biocl.train = exclude(biocl.train, vif)
        
        # Crop to the projection extent
        biocl.proje = raster::crop(x = bioclim,
                                   y = exten.proje)
        
        # Exclude variables with high VIF
        biocl.proje = exclude(biocl.proje, vif)
        
        # Convert to Raster Stack
        bioclim = as(bioclim, "Raster")
        biocl.train = as(biocl.train, "Raster")
        biocl.proje = as(biocl.proje, "Raster")
        
        # Mask to the projection contour
        biocl.proje = raster::mask(x = biocl.proje,
                                   mask = exten.proje)
        # Clean environment
        rm(gbif)
      }
      
      # STEP 3: BACKGROUND DATA AND OCCURRENCE DATA SPLIT
      incProgress(0.15, detail = "Preparing background and occurrence data...")
      
      {
        # Delimit species accessible area
        
        # Transform to SpatialPointsDataFrame
        coordinates(occur) = ~ Longitude + Latitude
        
        # Set the CRS of the environmental layers
        proj4string(occur) = proj4string(bioclim)
        
        # Drawn circles around occurrences
        acces.area = circles(p = occur,
                             d = 250000, # Circle radius in metres
                             lonlat = T)
        
        # Convert to a polygon
        acces.area = polygons(acces.area)
        
        # Note that for the paper's example we used a sampling bias layer
        
        # # Mask the spatial sampling bias raster with the species accessible area
        # acce.area = raster::mask(backgrou,
        #                          acce.area)
        
        # Mask the climatic raster with the species accessible area
        acces.area = raster::mask(bioclim[[1]],
                                  acces.area)
        
        # Drawn background points
        background = sampleRast(x = acces.area,
                                n = 10000,
                                replace = T,
                                prob = F)
        
        # Convert to data frame
        background = data.frame(background)
        
        # Transform to SpatialPointsDataFrame
        coordinates(background) = ~ x + y
        
        # Set the CRS of the environmental layers
        proj4string(background) = proj4string(bioclim)
        
        # Split training and validation datasets
        samp = sample(x = nrow(occur@coords), # Random sampling
                      size = floor(nrow(occur %>% as.data.frame()) * .2)) # Sampling 20% of the data
        
        # Training dataset
        occur.train = occur[-samp, ]
        
        # Validation dataset
        occur.valid = occur[samp, ]
        
        # Clean environment
        rm(samp, acces.area)
      }
      
      # STEP 4: RUN MAXENT
      incProgress(0.15, detail = "Preparing background and occurrence data...")
      
      {
        # Train a MaxEnt model for the historical period
        maxent = maxent(x = biocl.train, # Predictor layer
                        p = occur.train, # Presences
                        a = background) # Background layer
        
        # Extract permutation importance
        permu.impor = maxent@results %>%
          data.frame() %>%
          rownames_to_column("bioclim") %>%
          filter(str_detect(bioclim, ".permutation.importance")) %>%
          mutate(bioclim = str_replace(bioclim, ".permutation.importance", "")) %>%
          dplyr::rename("permutation.importance" = ".") %>%
          mutate(permutation.importance = permutation.importance/100)
        
        # Project the model
        projection = raster::predict(object = biocl.proje, # Reduced resolution
                                     model = maxent) # Project trained model
        
        # Reduce to two decimals
        projection = round(x = projection,
                           digits = 2)
        
        # Remove decimals
        projection = projection * 100
        
        # Threshold the model
        biocl.valid = raster::extract(x = biocl.train, # Climate data at the validation occurrences
                                      y = occur.valid)
        
        # Extract suitability
        suita = dismo::predict(object = maxent,
                               x = biocl.valid)
        
        # Reduce to two decimals
        suita = round(x = suita,
                      digits = 2)
        
        # Remove decimals
        suita = suita * 100
        
        # Determine threshold
        thres = quantile(suita, .1, na.rm = T)
        
        # Model spatial restriction
        
        # Drawn circles around occurrences
        acces.polyg = circles(p = occur,
                              d = 250000, # Circle radius in metres
                              lonlat = T)
        
        # Convert to a polygon
        acces.polyg = polygons(acces.polyg)
        
        # Binarize
        proje.binar = projection >= thres
        
        # Transforme unsuitable areas into NAs
        proje.binar[proje.binar[] == 0] = NA
        
        # Detect patches
        proje.patch = raster::clump(proje.binar)
        
        # Extract valid patches (i.e., patches touching or falling under the species accessible area)
        valid.patch = raster::extract(proje.patch, acces.polyg) %>%
          unique() %>%
          unlist() %>%
          as.numeric()
        
        # Transform invalid patches into NAs
        proje.patch[!proje.patch[] %in% valid.patch] = NA
        
        # Create the restrict projection object
        proje.restr = !is.na(proje.patch)
        
        # Transform 0 into NAs
        proje.restr[proje.restr[] == 0] = NA
      }
      
      # STEP 5: DOWNLOAD FUTURE CLIMATIC DATA
      incProgress(0.15, detail = "Downloading future climate data...")
      
      {
        # Top four ESM were selected based on the ISIMIP3b bias adjustment fact sheet (page 8): https://www.isimip.org/gettingstarted/isimip3b-bias-adjustment/
        
        # UKESM1-0-LL, SSP3, 2081-2100
        
        # Download CMIP6 data
        biocl.UKESM1 = cmip6_world("UKESM1-0-LL", ssp, time,
                                   var = "bioc", res = 10, path = tempdir()) # Resolution of 10 minutes of a degree
        
        # Convert to Raster Stack
        biocl.UKESM1 = as(biocl.UKESM1, "Raster")
        
        # Exclude variables with high VIF
        biocl.UKESM1 = exclude(biocl.UKESM1, vif)
        
        # Crop and mask
        biocl.UKESM1 = raster::mask(x = raster::crop(x = biocl.UKESM1,
                                                     y = exten.proje),
                                    mask = exten.proje)
        
        # MPI-ESM1-2-HR, SSP3, 2081-2100
        
        # Download CMIP6 data
        biocl.MPI = cmip6_world("MPI-ESM1-2-HR", ssp, time,
                                var = "bioc", res = 10, path = tempdir()) # Resolution of 10 minutes of a degree
        
        # Convert to Raster Stack
        biocl.MPI = as(biocl.MPI, "Raster")
        
        # Exclude variables with high VIF
        biocl.MPI = exclude(biocl.MPI, vif)
        
        # Crop and mask
        biocl.MPI = raster::mask(x = raster::crop(x = biocl.MPI,
                                                  y = exten.proje),
                                 mask = exten.proje)
        
        # IPSL-CM6A-LR, SSP3, 2081-2100
        
        # Download CMIP6 data
        biocl.IPSL = cmip6_world("IPSL-CM6A-LR", ssp, time,
                                 var = "bioc", res = 10, path = tempdir()) # Resolution of 10 minutes of a degree
        
        # Convert to Raster Stack
        biocl.IPSL = as(biocl.IPSL, "Raster")
        
        # Exclude variables with high VIF
        biocl.IPSL = exclude(biocl.IPSL, vif)
        
        # Crop and mask
        biocl.IPSL = raster::mask(x = raster::crop(x = biocl.IPSL,
                                                   y = exten.proje),
                                  mask = exten.proje)
        
        # MRI-ESM2-0, SSP3, 2081-2100
        
        # Download CMIP6 data
        biocl.MRI = cmip6_world("MRI-ESM2-0", ssp, time,
                                var = "bioc", res = 10, path = tempdir()) # Resolution of 10 minutes of a degree
        
        
        # Convert to Raster Stack
        biocl.MRI = as(biocl.MRI, "Raster")
        
        # Exclude variables with high VIF
        biocl.MRI = exclude(biocl.MRI, vif)
        
        # Crop and mask
        biocl.MRI = raster::mask(x = raster::crop(x = biocl.MRI,
                                                  y = exten.proje),
                                 mask = exten.proje)
        
        # Number of independent bioclimatic variables
        biocl.numbe = as.numeric(length(vif@results[["Variables"]]))
        
        # Empty raster stack - mean
        biocl.futur.mean = raster::stack()
        
        # Empty raster stack - standard deviation
        biocl.futur.sd = raster::stack()
        
        for (i in c(1:biocl.numbe)) {
          
          # Calculate mean
          biocl.futur.mean.loop = raster::stack(biocl.UKESM1[[i]],
                                                biocl.IPSL[[i]],
                                                biocl.MPI[[i]],
                                                biocl.MRI[[i]]) %>%
            raster::calc(fun = mean)
          
          # Calculate sd
          biocl.futur.sd.loop = raster::stack(biocl.UKESM1[[i]],
                                              biocl.IPSL[[i]],
                                              biocl.MPI[[i]],
                                              biocl.MRI[[i]], na.rm = T) %>%
            raster::calc(fun = sd)
          
          # Add the layer to to stack - mean
          biocl.futur.mean = raster::stack(biocl.futur.mean, biocl.futur.mean.loop)
          
          # Add the layer to to stack - sd
          biocl.futur.sd = raster::stack(biocl.futur.sd, biocl.futur.sd.loop)
          
        }
        
        # Rename it - mean
        names(biocl.futur.mean) = names(biocl.proje)
        
        # Rename it - sd
        names(biocl.futur.sd) = names(biocl.proje)
        
        # Clean environment
        rm(biocl.UKESM1, biocl.IPSL, biocl.MPI, biocl.MRI,
           biocl.numbe, biocl.futur.mean.loop, biocl.futur.sd.loop)
      }
      
      # STEP 6: RUN COSST
      
      {
        # Extract weights, i.e., MaxEnt permutation importance
        weights = permu.impor$permutation.importance
        
        # Normalize it
        weights = weights / sum(weights)
        
        # Add a CRS
        crs(site) = proj4string(biocl.proje)
        
        # Extract future climate at the site
        biocl.site = raster::extract(x = biocl.futur.mean,
                                     y = site) %>%
          as.numeric()
        
        # Calculate the climatic match
        clima.match = abs(biocl.proje - biocl.site) # Absolute values
        
        # Normalize it
        clima.match = (clima.match - clima.match@data@min) / (clima.match@data@max - clima.match@data@min)
        
        # Inverse it
        clima.match = (clima.match - 1) * -1
        
        # Weight it
        clima.match = clima.match * weights
        
        # Sum it
        clima.match = sum(clima.match)
        
        # Round it
        clima.match = round(clima.match, 5)
        
        # Crop and mask it
        clima.match = raster::mask(x = raster::crop(x = clima.match,
                                                    y = proje.restr),
                                   mask = proje.restr)
        
        # Normalize it
        clima.match = (clima.match - clima.match@data@min) / (clima.match@data@max - clima.match@data@min)
        
        # Climatic uncertainty
        
        # Round it
        clima.match.sd = round(biocl.futur.sd, 5)
        
        # Normalize it
        clima.match.sd = (clima.match.sd - clima.match.sd@data@min) / (clima.match.sd@data@max - clima.match.sd@data@min)
        
        # Round it
        clima.match.sd = round(clima.match.sd, 5)
        
        # Weight it
        clima.match.sd = clima.match.sd * weights
        
        # Sum it
        clima.match.sd = sum(clima.match.sd)
        
        # Round it
        clima.match.sd = round(clima.match.sd, 5)
        
        # Crop and mask it
        clima.match.sd = raster::mask(x = raster::crop(x = clima.match.sd,
                                                       y = proje.restr),
                                      mask = proje.restr)
        
        # Normalize it
        clima.match.sd = (clima.match.sd - clima.match.sd@data@min) / (clima.match.sd@data@max - clima.match.sd@data@min)
        
        # Inverse it
        clima.match.sd = 1 - clima.match.sd
        
        # Extract climate uncertainty at the site
        clima.uncer = raster::extract(x = clima.match.sd,
                                      y = site)
        
        # Euclidean distance
        
        # Calculate distance
        geogr.proxi = raster::distanceFromPoints(object = proje.restr,
                                                 xy = site)
        
        # Crop and mask it
        geogr.proxi = raster::crop(x = raster::mask(x = geogr.proxi,
                                                    mask = proje.restr),
                                   y = proje.restr)
        
        # Normalize it
        geogr.proxi = (geogr.proxi - geogr.proxi@data@min) / (geogr.proxi@data@max - geogr.proxi@data@min)
        
        # Inverse it
        geogr.proxi = 1 - geogr.proxi
        
        # Climate-adjusted
        clima.adjus = (clima.match * clima.uncer) + geogr.proxi
        
        # Normalize it
        clima.adjus = (clima.adjus - clima.adjus@data@min) / (clima.adjus@data@max - clima.adjus@data@min)
        
        # Stack it provenance projections
        prove.proje = raster::stack(geogr.proxi,
                                    clima.adjus,
                                    clima.match)
        
        # Rename it
        names(prove.proje) = c('compo', 'clima.adjus', 'predi')
      }
      
      # STEP 7: EXTRACT COSST INDEX AT SEED SOURCING SITES
      
      {
        
        # Define cosst.supp as a reactive expression
        cosst.supp <- reactive({
          req(input$file1)
          
          # Read input CSV
          seed.supp = read.csv(input$file1$datapath)
          
          # Create a new object
          seed.supp.coor = seed.supp
          
          # Make it spatial
          coordinates(seed.supp.coor) = ~ long + lat 
          
          # Extract COSST index
          cosst.supp = raster::extract(x = prove.proje[[provenance]],
                                       y = seed.supp.coor) %>%
            data.frame() %>%
            rename("cosst" = ".")
          
          # Transform into data frame and merge
          cosst.supp = seed.supp %>%
            data.frame() %>%
            cbind(cosst.supp) %>%
            mutate(cosst.perce = cosst / sum(cosst, na.rm = T) * 100) %>%
            dplyr::rename("cosst.index" = cosst,
                          "seed.mix.contribution.%" = cosst.perce)
          
          # Clean environement
          rm(seed.supp.coor)
          
          return(cosst.supp)
        })
      }
      
      # Convert site to data.frame for leaflet
      site_df <- as.data.frame(site)
      
      output$cosstMap <- renderLeaflet({
        leaflet() %>%
          addTiles() %>%
          addRasterImage(prove.proje[[provenance]], 
                         colors = colorNumeric(viridis_pal(option = "turbo")(256), 
                                               domain = c(0, 1), 
                                               na.color = "transparent"), 
                         opacity = input$opacity) %>%
          addCircleMarkers(data = site_df, 
                           lng = ~long, 
                           lat = ~lat, 
                           color = "black", 
                           radius = 5, 
                           fill = TRUE, 
                           fillColor = "black", 
                           fillOpacity = 1, 
                           popup = ~paste("Site Location: ", lat, ", ", long))
      })
      
      # Observe clicks on the map
      observeEvent(input$cosstMap_click, {
        click <- input$cosstMap_click
        
        # Convert the click to spatial points
        click_coords <- data.frame(long = click$lng, lat = click$lat)
        coordinates(click_coords) <- ~long + lat
        proj4string(click_coords) <- proj4string(prove.proje[[provenance]])
        
        # Extract the raster value at the clicked location
        raster_value <- raster::extract(prove.proje[[provenance]], click_coords)
        
        # Display the raster value in a modal dialog box
        showModal(modalDialog(
          title = "COSST index at clicked location",
          paste("Latitude:", round(click$lat, 5)),
          paste("Longitude:", round(click$lng, 5)),
          paste("COSST index:", ifelse(is.na(raster_value), "No Data", round(raster_value, 5))),
          easyClose = TRUE,
          footer = NULL
        ))
      })
      
      # Render table from cosst.supp reactive expression
      output$cosst_table <- renderTable({
        req(cosst.supp())
        cosst.supp()  # Display cosst.supp data frame
      })
      
      # Final update to progress
      incProgress(1, detail = "COSST Application Completed.")
    })
  })
}

shinyApp(ui = ui, server = server)
