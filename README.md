README Document for Reproducing Results in
==========================================
"Identification and Estimation of Discrete Choice Models with Unobserved Choice Sets"
=============================================
Victor H. Aguiar
vaguiarl@sfu.ca

Nail Kashaev
nkashaev@uwo.ca

Software
========

Version 1.5.3 of the `Julia` programming language was used to code the analysis files. (Version 1.10.1 was used for plotting the results.) For details about how to install `Julia` on different platforms and how to make `Julia` programs executable from the command line see <https://julialang.org/downloads/platform/>. After installation of `Julia`, run `using Pkg` followed by `Pkg.instantiate()` in the `Julia` terminal after setting the replication folder as the main one.

Simulations and estimations use `KNITRO 12.3.0.`.  

Data
====

Our empirical application is based on data from Nielsen Consumer LLC and marketing databases provided through the NielsenIQ Datasets at the Kilts Center for Marketing Data Center at The University of Chicago Booth School of Business.
We cannot include it to this package. For the information on acquiring the data for academic researchers see <https://www.chicagobooth.edu/research/kilts/datasets/nielseniq-nielsen>. 

`Data/extract_2016_2018_34.csv` only contains the household and market identifiers sufficient for reconstructing the data.

Hardware
========

- The code was run on Mac Studio (M2 Ultra with 128 Gb of RAM), MacOS 14.3.1 (23D60).


Content
=======

-   `Application`  -- this folder contains the analysis files to replicate the results presented in Section 5 of the main text and Online Appendix C.

-   `Simulations`  -- this folder contains the analysis files to replicate the results in Online Appendix B.

-   `Manifest.toml` and `Project.toml`  -- toml files with all necessary `Julia` packages.



Below, we describe the contents of every folder.

`Applications`
============

-   `Data` --  this folder contains one csv file (`extract_2016_2018_34.csv`) with all the data needed for estimation. 
    In the file, we provide all entries except for `household_code`, `tripidy`, and `ziphhn`, are set to 0. This information should be enough to reconstruct the data from the Nielsen HomeScan.
-   `Results_prelim` --  this folder contains the output from running the code in `decomp_app.jl`.
-   `Tables_and_figures` --  this folder contains all tables and figures used in the main paper and the online appendix.
-   `dataandest_app.jl` -- this file constructs tables and figures.
-   `decomp_app.jl` -- this file estimates consideration sets and the distribution of choices conditional on consideration sets.
-   `elasticities_app.jl` -- this file estimates the elasticities reported in Table 2 in the main text and Table 14 in thr online appendix.
-   `functions_dataandest_app.jl` -- this file contains functions used in `dataandest_app.jl`.
-   `functions_decomp_app.jl` -- this file contains functions used in `decomp_app.jl`.
-   `functions_elasticities_app.jl` -- this file contains functions used in `elasticities_app.jl`.


`Simulations`
============

-   `Recovering_FM` -- this folder contains the files needed to replicate Tables 3-10.
    - `Results_sim2` -- this folder contains the results of estimation: F_model_N.csv and M_model_N.csv, where model $\in\{nested, core\}$ and N $\in\{2000,5000,10000,50000\}$. 
    - `Tables` -- this folder contains the resulting tables.
    - `sim2_functions.jl` -- this file contains functions used in `sim2_main.jl`.
    - `sim2_main.jl` -- this file conducts simulations.
    - `sim2_results.jl`-- the file generates all tables.

-   `Recovering_Sets` -- the folder contains the files needed to replicate Tables 1-2.
    - `Results_sim1` -- this folder contains the results of estimation: Ms1_model_N.csv and Ms2_model_N.csv, where model $\in\{nested, core\}$ and N $\in\{2000,5000,10000,50000\}$. 
    - `Tables` -- this folder contains the resulting tables.
    - `sim1_functions.jl` -- this file contains functions used in `sim1_main.jl`.
    - `sim1_main.jl` -- this file conducts simulations.
    - `sim1_results.jl`-- the file generates all tables.
