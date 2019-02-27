# UrbanGW
This repository takes raster and csv input data for a groundwater model of the Valley of Mexico and runs a simulation model to compare or optimize recharge scenarios under multiple objectives.

Dependencies:
```flopy```
```osgeo```
```numpy```

## Data Inputs

![Data](/images/data_processing.png)

The data for the model must be in the following format:

Data Set Required | Format
-------------------- | --------------------
Digital elevation model of model area | raster
Active model cells by model layer | raster
Layer thickness by model layer | raster
Initial head for model area | raster
Percent of grid cell for each land use type (natural, urban, water) | raste
Interpolated monthly precipitation | raster
Pumping data by month or year | csv
Area covered by each municipality | raster
Leaks as percentage of total usage by municipality | csv
List of all potential recharge basin sites | csv
List and characteristics of existing wastewater treatment plants | csv

### Ensure rasters are at the model grid resolution
```
ValleMexico_data_main.py
```

## Simulation Options
Each model run will create an object with the following results:

 ```self.wwtps``` is a list of randomly selected wastewater treatment plants where reuse has been implemented. ```self.basins``` is a list of the row and column where each infiltration basin has been implemented. ```self.mthlyleak``` is an array with the total quantity of leaks in the system per month in m3 cost is the total number of interventions time their weights defined in the model set-up. ```self.cost``` is a summed relative cost of the scenario based on the interventions applied. ```self.wells``` is a dictionary of well objects input into the MODFLOW model which includes positive flow from wastewater treatment plants, leaks, and recharge basins and negative flow from pumping wells. ```self.landuse``` is a dictionary that contains a raster and list of the percentage of land use type (NATURAL, URBAN, or WATER) per model cell for each model phase.

```
main.py
```

### Test Mode


### Scenario Mode

Scenario mode allows the user to compare predetermined scenarios defined above by the quantity of each recharge intervention.

```
plt_scen = True
run_scenarios = True
```

#### Scenario options

#### Plotting


### Optimization Mode


