## code to prepare `ebola_gueckedou` dataset goes here

data_tidy <- read_csv(
      "Guinea_modified.csv",
      col_types = readr::cols(
            .default            = readr::col_double(),
            Location            = readr::col_character(),
            `Ebola data source` = readr::col_character(),
            `Indicator type`    = readr::col_character(),
            `Case definition`   = readr::col_character()
            )
      ) %>%
      rename(
            location=Location,
            source_data = `Ebola data source`,
            type_data = `Indicator type`,
            definition_data = `Case definition`
            ) %>%
      pivot_longer(
            `30 December 2013 to 05 January 2014 (2014-W01)`:`18 to 24 May 2015 (2015-W21)`,
            names_to = "time",
            values_to = "cases"
            )

data_gueckedou <- data_tidy %>%
      filter(location == "GUECKEDOU") %>%
      group_by(time) %>%
      summarize(cases = sum(cases))

I_k <- data_gueckedou$cases

date_start <- dmy("30/12/2013")
date_end   <- dmy("18/05/2015")
dates      <- seq(date_start, date_end   , by = 7) # use for plotting
ts         <- seq(0         , length(I_k), by = 7) # use in MCMC algorithm



ebola_gueckedou <- list(
      I_k=I_k,
      dates=dates,
      ts=ts
      )

usethis::use_data(ebola_gueckedou, overwrite = TRUE)
