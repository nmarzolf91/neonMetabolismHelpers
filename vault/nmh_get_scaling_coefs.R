#' Download hydraulic scaling equations determined for NEON wadable sites from my personal github repo
#' 
#' @author Nick Marzolf, \email{nick.marzolf@@jonesctr.org}
#' @author Wes Slaughter, \email{wslaughter@@berkeley.edu}

nmh_get_scaling_coefs <- function(){
  
  # this function requries the Field Discharge and Reaeration Field and Lab data products
  # If the nmh_get_neon_data() function is used for those data products, they will be saved here
  dir <- 'data/raw/neon/'
  
  # find the data needed
  dsc_dir <- glue::glue(dir, 'Discharge field collection/')
  rea_dir <- glue::glue(dir, 'Reaeration field and lab collection/')
  
  # this is only for the 24 wadable sites
  sites <- list.files(rea_dir)
  
  # create an empty list for the for loop
  all <- list()
  
  # for loop around sites
  for(i in 1:length(sites)) {
    
    # define sites
    site <- sites[i]
    
    # widths <- feather::read_feather(glue::glue(rea_dir,'{site}/rea_widthFieldData.feather')) 
    
    # mean_widths <- widths %>% 
    #   dplyr::mutate(eventID = paste(siteID, lubridate::date(collectDate), sep = '_')) %>% 
    #   dplyr::arrange(collectDate, widthMeasurementNumber) %>% 
    #   dplyr::group_by(eventID) %>% 
    #   dplyr::summarise(n_meas = n(),
    #                    mean_width = mean(wettedWidth, na.rm = TRUE),
    #                    median_width = median(wettedWidth, na.rm = TRUE),
    #                    min_width = min(wettedWidth, na.rm = TRUE),
    #                    max_width = max(wettedWidth, na.rm = TRUE),
    #                    sd_width = sd(wettedWidth, na.rm = TRUE),
    #                    se_width = sd_width/sqrt(n_meas))
    
    # read in field discharge data
    event_dsc_fieldData <- try(feather::read_feather(
      glue::glue(dsc_dir, site, '/dsc_fieldData.feather')) %>% 
        dplyr::mutate(eventID = paste(siteID, lubridate::date(collectDate), sep = '_'))
    )
    
    # error catching
    if(inherits(event_dsc_fieldData, 'try-error')) {
      cat(glue::glue('No dsc field data from eventID: {site}')) 
    }
    
    # read in cross-section data
    event_dsc_individualFieldData <- try(feather::read_feather(
      glue::glue(dsc_dir, site, '/dsc_individualFieldData.feather')) %>% 
        dplyr::mutate(eventID = paste(siteID, lubridate::date(collectDate), sep = '_'))
    )
    # and error handling
    if(inherits(event_dsc_individualFieldData, 'try-error')) {
      cat(glue::glue('No dsc_individualFieldData from eventID: {site}'))
    }
    
    # event_dsc_fieldData_calc <- try(conv.calc.Q(stageData = event_dsc_fieldData,
    #                                             dischargeData = event_dsc_individualFieldData) %>% 
    #                                   dplyr::mutate(calcQ_m3s = calcQ/1000,
    #                                                 eventID = paste(siteID, lubridate::date(collectDate), sep = '_')) %>% 
    #                                   dplyr::select(eventID, calcQ_m3s)
    # )
    # if(inherits(event_dsc_fieldData_calc, 'try-error')){
    #   cat(glue::glue('{site} has no field discharge data'))
    #   next
    # }
    
    # run another loop around each discharge measuring event
    out <- list()
    for(j in 1:length(unique(event_dsc_individualFieldData$recorduid))){
      
      # define the event
      event <- event_dsc_fieldData[j,]
      
      # pull out the sections from the cross-section dataset
      sections <- event_dsc_individualFieldData %>% 
        dplyr::filter(recorduid %in% event$recorduid) %>% 
        dplyr::arrange(stationNumber)
      
      # pull in discharge
      discharge_m3s <- event$finalDischarge/1000
      
      # save the summary from the cross-section with
      ## width, mean depths, mean velocity, and discharge
      out[[j]] <- sections %>% 
        dplyr::group_by(eventID) %>% 
        dplyr::summarise(n_stations = length(stationNumber),
                         width_m = max(tapeDistance),
                         mean_depth_m = mean(waterDepth, na.rm = TRUE),
                         mean_velocity_ms = mean(averageVelocity, na.rm = TRUE)) %>% 
        dplyr::mutate(discharge_m3s = discharge_m3s,
                      site = site, 
                      date = lubridate::date(substr(eventID, 6, 30)))
    } # end event for loop
    
    # bind the event outputs
    all[[i]] <- do.call(rbind, out)
    
    
  } # end for loop
  
  # bind the site data together
  all <- do.call(rbind, all)
  
  # calculate hydraulic scaling for NEON sites
  
  # matrix with the coefficients + errors to populate within the loop
  mat <- matrix(ncol = 18,
                nrow = length(unique(all$site)))
  colnames(mat) <- c('site', 'a', 'a_se', 'b', 'b_se', 'c', 'c_se', 'd', 'd_se', 'e', 'e_se',
                     'f', 'f_se','r2_width', 'r2_vel', 'r2_depth', 'prod_coefs', 'sum_coefs')
  
  # for loop each site
  for(k in 1:length(unique(all$site))) {
    
    # define the site and subset the dataframe
    site <- unique(all$site)[k]
    
    site_df <- all %>% 
      dplyr::filter(site %in% !!site,
                    mean_velocity_ms > 0)
    
    # power law scaling for depth, width, and velocity
    lm_depth <- lm(data = site_df %>% 
                     filter(is.finite(log(mean_depth_m)),
                            is.finite(log(discharge_m3s))
                     ),
                   log(mean_depth_m) ~ log(discharge_m3s))
    
    c <- as.numeric(summary(lm_depth)$coefficients[1])
    c_se <- as.numeric(summary(lm_depth)$coefficients[3])
    d <- as.numeric(summary(lm_depth)$coefficients[2])
    d_se <- as.numeric(summary(lm_depth)$coefficients[4])
    r2_depth <- summary(lm_depth)$adj.r.squared
    
    lm_width <- lm(data = site_df %>% 
                     dplyr::filter(is.finite(log(width_m)),
                                   is.finite(log(discharge_m3s))
                     ),
                   log(width_m) ~ log(discharge_m3s))
    a <- as.numeric(summary(lm_width)$coefficients[1])
    a_se <- as.numeric(summary(lm_width)$coefficients[3])
    b <- as.numeric(summary(lm_width)$coefficients[2])
    b_se <- as.numeric(summary(lm_width)$coefficients[4])
    r2_width <- summary(lm_width)$adj.r.squared
    
    lm_vel <- lm(data = site_df %>% 
                   dplyr::filter(is.finite(log(mean_velocity_ms)),
                                 is.finite(log(discharge_m3s))
                   ),
                 log(mean_velocity_ms) ~ log(discharge_m3s))
    
    e <- as.numeric(summary(lm_vel)$coefficients[1])
    e_se <- as.numeric(summary(lm_vel)$coefficients[3])
    f <- as.numeric(summary(lm_vel)$coefficients[2])
    f_se <- as.numeric(summary(lm_vel)$coefficients[4])
    r2_vel <- summary(lm_vel)$adj.r.squared
    
    # Power law math
    sum_coefs <- b+d+f
    prod_coefs <- a*c*e
    
    # populate the output matrix
    mat[k,1] <- site
    mat[k,2] <- exp(a)
    mat[k,3] <- a_se
    mat[k,4] <- b
    mat[k,5] <- b_se
    mat[k,6] <- exp(c)
    mat[k,7] <- c_se
    mat[k,8] <- d
    mat[k,9] <- d_se
    mat[k,10] <- exp(e)
    mat[k,11] <- e_se
    mat[k,12] <- f
    mat[k,13] <- f_se
    mat[k,14] <- r2_width
    mat[k,15] <- r2_vel
    mat[k,16] <- r2_depth
    mat[k,17] <- prod_coefs
    mat[k,18] <- sum_coefs
    
  } # end for loop
  
  # prepare the output table
  mat <- data.frame(mat) %>% 
    mutate_at(c(2:ncol(mat)), as.numeric) %>% 
    format(., digits = 4)
  
  readr::write_csv(mat,
                   'data/3_NEON_hydraulic_coefs.csv')
  
  return(mat)
}

