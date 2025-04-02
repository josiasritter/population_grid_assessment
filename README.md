# population_grid_assessment
This repository contains the code and research data of the following study, evaluating the accuracy of global gridded population datasets in 307 rural areas around the world:

"Global gridded population datasets systematically underrepresent rural population".
Josias Láng-Ritter, Marko Keskinen, Henrikki Tenkanen.
To appear in Nature Communications, 2025.

Steps to follow:
1. multipops_all.py is the main program that estimates the population amounts in the rural areas for five global gridded population datasets.
2. results.geojson is the resulting polygon file that includes the rural areas evaluated, the estimated rural populations, and the reported rural populations.
3. multipops_validation_natcom.py is the program that applies the bias correction to account for underrepresentation of reservoir areas in GeoDAR, carries out the validation procedure against reported rural populations, and generates figures in the Results section of the paper.
