## CI tests for MYNN-EDMF
* Cases include:
    * clr, land point (43.40°N, 73.02°W) on Aug 23, 2024 

* Input:

    * Input file for clr case is generated using WRFv4.5.1

* To run the CI tests offline:
    ```
    export NFHOME=/path/to/your/netcdf/dir
    cd tests
    make
    ./driver
    ```
