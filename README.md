advanced R code review

title: C&O Occupancy by Little Brown Bats
author: Megan Moran

First part of code is organizing data and putting it in the proper format to run occupancy models in unmarked. The unmarked package needs to include:
1. a detection history (1, 0, or NA for species detection at each site, a matrix) = y
2. site covariates (variables that differ by site - forest cover, elevation, etc., a dataframe or vector) = siteCovs
3. observation covariates (variables that differ by site and observation - often things like temperature, precipitation, etc. in my case I use julian day and Year, a matrix) = obsCovs  

The second part of the code is running occupancy models with unmarked

