# Idaho_Predator_CoOccurrence
R code and anonymized data associated with the publication S. B. Bassing, D. E. Ausband, M. Mumma, S. Thompson, M. A. Hurley, and M. R. Falcy. Mammalian predator co-occurrence affected by prey and habitat more than competitor presence at multiple time scales. Ecological Monographs.

This repository archives the code used for analyses described in the manuscript. Camera trap locations and images are sensitive and cannot be provided publicly. Qualified researchers can contact Matt Boone (idfgdatarequests@idfg.idaho.gov), Data Management Lead, Idaho Department of Fish and Game, Idaho Fish and Wildlife Information System, 600 S Walnut, Boise, ID 83712 for full data requests. Query details: Camera trap data collected as part of the North Idaho Predator-Prey Project, collected from July 1 – Sept. 15, 2020 and 2021 in Game Management Units (GMU) 1, 6, and 10A, including camera ID and location coordinates, camera setup style, camera deployment and retrieval dates, and detections of black bears, bobcats, coyotes, mountain lions, wolves, elk, moose, rabbit-hare, and white-tailed deer, including date, time, and trigger method of individual images.

Data description:
Covariate data: Camera_Deployment_Data, Covariates_Smr20 & Covariates_Smr21
NewLocationID = unique identifier for each camera based on deployment location (GMU), camera setup (P vs U), and camera number
GMU = Idaho Department of Fish & Game Game Management Unit where sampling occurred (GMU 1, GMU 6, and GMU 10A)
Setup = camera setup style (ungulate = random location, 1.5 m high and perpendicular to ground; predator = on dirt-bottomed linear feature, >2.5 m high and angled towards ground)
Season = Summer of data collection (Smr20 = summer 2020, Smr21 = summer 2021)
CameraHeight_M = height of camera from ground (m)
CameraFacing = general feature camera is view-shed focuses on (e.g., road, decommissioned road, hiking trail, game trail, random)
perc_forest = percentage of forested habitat within 100m of a camera site
Elevation__10m2 = elevation (m) derived from a Digital Elevation Model (10 m-sq pixel resolution)
TRI = Terrain Ruggedness Index representing terrain variability at camera site, derived from a Digital Elevation Model 
elk_perday = number of independent detection events where elk were detected per day
lagomorphs_perday = number of independent detection events where lagomoprhs (hares and rabbits) were detected per day
moose_perday = number of independent detection events where moose were detected per day
whitetaileddeer_perday = number on independent detection events where white-tailed deer were detected per day

Detection histories: DetectionHist_Smr20 & DetectionHist_Smr21
Contains the binary detection/non-detection data for each species at each camera site during each 7-day sampling occasion in summer, July 1 - September 15, 2020 and 2021
