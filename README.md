# UrbanGW
This repository takes raster and csv input data for a groundwater model of the Valley of Mexico and runs a simulation model to compare or optimize recharge scenarios under multiple objectives. The portfolio of intervention options includes repairs to the distribution network to reduce leaks, increased wastewater treatment and infiltration, and increased infiltration of imported water using recharge basins. These policies can be optimized using the NSGAII MOEA according to four objectives: minimize energy use for pumping, maximize hydraulic head in high subsidence regions, minimize groundwater mounding in urban areas, and minimize the cost of the interventions.

Dependencies: ```flopy``` ```osgeo``` ```numpy``` ```pickle```

To run new simulations you must have MODFLOW installed in ```C:\WRDAPP\MF2005.1_12\bin\mf2005.exe``` which can be downloaded at https://water.usgs.gov/ogw/modflow/mf2005.html#downloads. However, you can still plot saved files when ```run_scenarios``` and ```run_optimization``` are set to ```False```.

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

The user must input the data above and integer choices for interventions to run the simulation model. The choices for interventions are defined as: the number of wastewater treatment plants to rehabilitate for wastewater injection into the aquifer from the existing plants, the number of infiltration basins that will recharge the aquifer using imported water from the potential sites, and the percent of fixed leaks when compared to historical leaks (0 indicates the same level as historical leaks and 100 indicates all leaks are fixed).

Each model run will create an object with the following results:

 ```self.wwtps``` is a list of randomly selected wastewater treatment plants where reuse has been implemented. ```self.basins``` is a list of the row and column where each infiltration basin has been implemented. ```self.mthlyleak``` is an array with the total quantity of leaks in the system per month in m3 cost is the total number of interventions time their weights defined in the model set-up. ```self.cost``` is a summed relative cost of the scenario based on the interventions applied. ```self.wells``` is a dictionary of well objects input into the MODFLOW model which includes positive flow from wastewater treatment plants, leaks, and recharge basins and negative flow from pumping wells. ```self.landuse``` is a dictionary that contains a raster and list of the percentage of land use type (NATURAL, URBAN, or WATER) per model cell for each model phase.

```
main.py
```

### Test Mode

Use this mode to play around with the interventions in a single model. (Plotting capabilities are not yet implemented)

### Scenario Mode

Scenario mode allows the user to compare predetermined scenarios defined by the quantity of each recharge intervention.

```
plt_scen = True
run_scenarios = True
```

#### Plotting
```plt_head_change``` plots a comparison of the head change during the model simulation between each scenario and the first "Historical" scenario. ```plt_scen_objectives``` plots a bar chart comparing each of the three objectives.

### Optimization Mode

The optimize option runs an optimization problem with the model using the recharge decisions and objectives defined above. If the ```run_optimization``` is ```False```, previous results are loaded from the file indicated.
