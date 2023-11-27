
# dynamicsthR: Conservation Action Planning with Dynamic Mathematical Modeling

## Overview

The `dynamicsthR` package is a powerful tool designed for deterministic dynamic mathematical modeling, specifically tailored for guiding conservation action planning. By utilizing integer programming models, this package simulates various threats, aiding decision-makers in determining optimal conservation actions, their locations, and timing.

## Overview

Ensure you have the "remotes" package installed. If not, install it using the following command:

``` r
if (!require(remotes)) install.packages("remotes")
remotes::install_github("josesalgr/dynamicsthR")

library(dynamicsthR)
```

## Usage

The package follows a clear three-step process:

1)  Function: **`create_instance()`**
    -   *`type`*`:` Specifies the simulation type, either "square" for grid simulation or "specific" when using inputs from a .dat file.
    -   `file` (optional): If type is "specific," this parameter specifies the path to the .dat file.

Example:

```{r, echo=FALSE}
instance <- dynamicsthR::create_instance(type = "square")
```

2.  Function: **`model()`**
    -   **`input`**: A list containing inputs generated by the **`create_instance()`** function.
    -   **`levels`**: Levels of threats.
    -   **`periods`**: Planning periods.
    -   **`budget`**: Annual budget for conservation actions.

Example:

```{r, echo=FALSE}
conservation_model <- dynamicsthR::model(
  input = instance,
  levels = 3,
  periods = 5,
  budget = 100000
)
```

3.  Solve the problem: `solve()`
    -   **`solver`**: Specifies the solver to be used, either "gurobi" (requires academic license) or "CBC" (free).

Example:

```{r, echo=FALSE}
solution <- dynamicsthR:solve(model = conservation_model, solver = "gurobi")

```

## **License**

This project is licensed under the [MIT License](https://chat.openai.com/c/LICENSE.md).