# Ontogeny_GF
Code used to study the ontogeny of greater flamingos originating from distinct Mediterranean breeding areas 
Reference Information
=====================

Provenance for this README
--------------------------

* File name: README.md
* Authors: Lina Lopez-Ricaurte
* Other contributors: Wouter M.G. Vansteelant, Arnaud Antoine, Olivier Duriez, Frédéric Jiguet, Sergio Nissardi, Davide Scridel, Lorenzo Serra, Stéphan Tillo, Arnaud Béchet, Jocelyn Champagnon
* Date created: 2025-05-15

- - -

Accompanying Paper and Data
---------------------------

* Paper Title: "Gradual ontogenetic shifts in the mobility and space use of long-lived migratory flamingos"
* Paper identifier: doi: https://doi.org/10.1101/2025.06.16.659883
* Dataset Title: Data for the article "Gradual ontogenetic shifts in the mobility and space use of long-lived migratory flamingos"
* Persistent Identifier: 
* Dataset Contributors:
* Creators: Wouter M.G. Vansteelant, Arnaud Antoine, Olivier Duriez, Frédéric Jiguet, Sergio Nissardi, Davide Scridel, Lorenzo Serra, Stéphan Tillo, Arnaud Béchet, Jocelyn Champagnon
* Date of Issue: 2025-05-15

- - -

Methodological Information
==========================

* Methods of data collection: see main manuscript for details

- - -

Instructions for replicating analyses
=====================================

Code/Software
-----------------
Code in this GitHub repository includes all steps to further process open access tracking data and metada and to replicate analyses and visualisatios.

Setup
-----
 * store all codes in a dedicated working directory
 * store tracking and metadata made available through Dryad in a subfolder 'Data' of your working directory 
 * open R Studio and create R project in the same working directory
 * run scripts one by one in the following order:

*1) 00_load_Packges.R:* loads all required packages and additional functions needed to run the analyses 

*2) 01_Fig_1_tracking_data.R:* visualises the tracking data included in this study

*3) 02_Fig2&models3d.R:* segments data into local versus travel events, calculates staging metrics, and models age using different functional relationships. GLMMs to test the effect of age and season on staging metrics and creates figure 2 

*4) 03_Sensitivity_analyses_staging_events.R:* use different temporal thresholds to quantify staging events and creates figure S2

*5) 04_Reproductive_status.R:* defines a 1 km buffer around each breeding colony and calculates the duration flamingos spend in the colony during the reproductive months

*6) 05_Traditional_sites.R:* Segments data into traditional sites and other locations, calculates durations, and models age using different functional relationships. Utilizes GLMMs to test the effect of age and site type on flamingo behavior across seasons. Generates Figure 3


Pre-processed tracking data and metadata on Dryad
==================================================

Details for: GF-fULL_resampled_v20250212.csv
---------------------------------------
* Description: a comma-delimited file containing the tracking data of 83 greater flamingos (Phoenicopterus roseus) tagged as fledglings in Southern France and Italy, resampled to an hourly resolution. 

* Format(s): .csv

* Size(s): 189.46 MB

* Variables:
  * dev: an 8-digit identifier for each device used in this study. Each value also represents a unique individual
  * date: timestamp in UTC
  * colony: origin of each bird
  * date_time: timestamp in UTC
  * LON: longtiude (decimal degrees, WGS84)
  * LAT: latitude (decimal degrees, WGS84)
  * spd_gps: recorded by GPS
  * colony.lat: latitude of the colony of origin (decimal degrees, WGS84)
  * colony.lat: longitude of the colony of origin (decimal degrees, WGS84)
  * year: year
  * sex: Female, Male
  * dist: calculated distance in meters between the current position to the previous one
  * dur: duaration in seconds
  * dir: direction
  * spd: speed m/s
  * spd.kmh: speed in kilometers 
  * julian: julian day
  * month: month
 

Details for 101 fledglings tagged during the breeding seasons of 2015 – 2024 and a list of the traditional sites are located in the supplementary material Table S1 and S2. 

- - -
END OF README
