## Citation:
Petracca LP, B Gardner, BT Maletzke, and SJ Converse. 2023. Merging integrated population models and individual-based models to project population dynamics of recolonizing species. Biological Conservation.

## Abstract:
Recolonizing species exhibit unique population dynamics, namely dispersal to and colonization of new areas, that have important implications for management. A resulting challenge is how to simultaneously model demographic and movement processes so that recolonizing species can be accurately projected over time and space. We introduce a framework for spatially explicit projection modeling that harnesses the rigorous parameter estimation made possible by an integrated population model (IPM) and the flexible movement modeling made possible by an individual-based model (IBM). Our framework has two components: [1] a Bayesian IPM-driven age- and state-structured population model that governs the population state process and estimation of demographic rates, and [2] an IBM-driven spatial model that allows for the projection of dispersal and habitat colonization. We applied this model framework to estimate current and project future dynamics of gray wolves (Canis lupus) in Washington State, USA. We used data from 74 telemetered wolves and yearly pup and pack counts to parameterize the model, and then projected statewide dynamics over 50 years. Mean population growth was 1.29 (95% Bayesian Credible Interval = 1.26-1.33) during initial recolonization from 2009-2020 and decreased to 1.02 (95% Prediction Interval = 0.98-1.04) in the projection period (2021-2070). Our results suggest that gray wolves have a ~100% probability of colonizing the last of Washington State’s three specified recovery regions by 2030, regardless of alternative assumptions about how dispersing wolves select new territories. Our spatially explicit projection model can be used to project the dynamics of any species for which spatial spread is an important driver of population dynamics.

## Code 
1) [Scripts](./scripts/): The main folder (described [here](./scripts/a_DESCRIPTION.txt)) contains three codes -- one to run the IPM for the data collection period (2009-2020) and two different versions of the projection model for the years 2020-2070. The latter scripts differ based on chosen territory selection method. 

The ["Bonus_Model_Formulation_with_Many_Ages" folder](./scripts/Bonus_Model_Formulation_with_Many_Ages) contains the script for a projection model formulation that is age-structured and removes wolves from the model after reaching 15 years of age. This script uses the categorical RSF territory selection process.

Lastly, this folder has four scripts in the ["Preliminary_Analysis_Steps" folder](./scripts/Preliminary_Analysis_Steps) for processing of GPS collar data for the territory size and RSF analyses. However, these data are not available publicly due to their sensitive nature. Interested parties should contact Donny Martorello at WDFW (Donny.Martorello@dfw.wa.gov).

2) [Functions](./functions/): The main folder (described [here](./functions/a_DESCRIPTION.txt)) contains all functions related to the individual-based movement model. These functions allow for attraction of lone wolves, lethal removals, and dispersal of wolves depending on chosen territory selection method. 

The ["Bonus_Model_Formulation_with_Many_Ages" folder](./functions/Bonus_Model_Formulation_with_Many_Ages) contains the same functions for a projection model formulation that is age-structured and removes wolves from the model after reaching 15 years of age. These functions use the categorical RSF territory selection process.

## Data
Data files used to run the models in this paper are found in the [data](./data) folder. Please see [here](./data/a_DESCRIPTION.txt) for more information.

## Updates
The scripts in this folder are the most recent, as reflected in the Corrigendum to our 2024 paper (https://www.sciencedirect.com/science/article/pii/S0006320724001952). We also fixed an additional index, whereby the impact of this correction (e.g., about one wolf after 50 years) was much smaller than simulation error. 