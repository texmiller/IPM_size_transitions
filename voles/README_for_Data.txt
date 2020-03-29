The data collection was led by Lowell L. Getz. Please refer to http://www.life.illinois.edu/getz/index.html for more information on the data file. The file presented here contains a part of this data as compiled by the authors of the associated manuscript.

The file data.zip contains two files:
1. Males.csv
Each row describes one capture of an animal

Columns:
Date - the date of the trapping session
Session - a unique number referring to the month of trapping
Phase - the assigned density phase in that session
id - the unique ID of the animal
mass - the body mass of the animal
stage - juvenile or adult
repro - 1: juvenile, 2: breeding, 4: non-breeding

2. Population.csv
This file contains one row for each trapping session.

Columns:
Date - the date of the trapping session
Session - a unique number referring to the month of trapping
Phase - the assigned density phase
Captured_individuals - the total number of individuals (m/f) captured in that session
Average_body_mass - the average body mass of all captured individuals (m/f)
Captured_male_juveniles - the number of male newborns that were captured in that session, estimated as total number of individuals between 15 and 25 grams times 3 (due to the restricted mass range) divided by two (to exclude the females newborns).
Captured_breeding_males - the number of reproductively active males
Average_juveniles_mass_(t+1) - average mass of the newborns
Average_breeding_male_mass_(t) - average mass of the reproductively active males
