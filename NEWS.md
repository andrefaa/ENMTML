# ENMTML 1.0.0

## Block cross-validation improved (24-02-2021)

-   Codes to perform block cross-validation (**part** argument of the **ENMTML** function) were improved to
    speed up the process.

-   Now will be tested 30 grid sizes. The minimum
    and maximum grid sizes to be tested are determined by 2 x environmental layer
    resolution and 100 x environmental layer resolution, respectively.

-   We expect to adapt the **part** argument so
    that users can select the minimum and maximum grid size and the number of grids
    to be tested.

New grid cell size for block cross-validation method .

## Major changes

-   Fix to an issue that caused parallel modelling to fail using macOS\
-   Pseudo-absence allocation debug

## Minor changes

-   Tests updated\
-   Vignettes included\
-   Creation of a website with guides to run the function and a FAQ\
-   Readme update
