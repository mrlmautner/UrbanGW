# UrbanGW
This repository takes raster and csv input data for a groundwater model of the Valley of Mexico and runs a simulation model to compare scenarios or optimize under multiple objectives.

Dependencies:
```flopy```
```osgeo```
```numpy```

## Data Inputs

![Data](/images/data_processing.png)

### Ensure rasters are at the model grid resolution
```
ValleMexico_data_main.py
```

## Simulation Options
The model run will create an object with the following results:

 '''self.wwtps''' is a list of randomly selected wastewater treatment plants where reuse has been implemented. '''self.basins''' is a list of the row and column where each infiltration basin has been implemented. '''self.mthlyleak''' is an array with the total quantity of leaks in the system per month in m3 cost is the total number of interventions time their weights defined in the model set-up
    self.cost - a summed relative cost of the scenario based on the interventions applied
    self.wells - dictionary of well objects input into the MODFLOW model which includes positive flow from wastewater treatment plants, leaks, and recharge basins and negative flow from pumping wells
    self.landuse - dictionary that contains a raster and list of the percentage of land use type (NATURAL, URBAN, or WATER) per model cell for each model phase

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


