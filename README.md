**Code for "Just in Time or Just a Guess? Addressing Challenges in Validating Prediction Models Based on Longitudinal Data"**

This repository contains the code used for the analyses presented in the paper *Just in Time or Just a Guess? Addressing Challenges in Validating Prediction Models Based on Longitudinal Data.*

The code provided here replicates the analysis and results discussed in the paper. Additionally, a Shiny app that complements our paper is **available in a separate repository**: <https://github.com/AnnaLangener/Justintime>.

-   **MAIN SCRIPT: Paper_Results.Rmd**: Markdown file providing a step-by-step walk through of the results and figures presented in the paper.

Simulation Study:

-   **Run_Simulation.R**: Script used to generate the results from the simulation study and execute the necessary functions.

-   **Simulation_Functions_notshiny.R**: Contains the core functions used in the simulation study.

    -   **simulation_results_0.csv**: Dataset containing the results from the simulations.

    -   **simulation_results_baseline.csv**: Dataset containing the results from the simulations for the baseline model.

Empirical Example: (Qwantify App dataset (esm_clean.csv from <https://osf.io/sxfrx/files/osfstorage>)

-   **EmpiricalExample.R**: Was used to prepare the results for the empirical binary example.

-   **Simulation_UploadData.R:** Contains the functions to run different cross-validation splits and to shuffle the data for binary outcomes (used in EmpiricalExample.R)

-   **EmpiricalExample_continious.R**: Was used to prepare the results for the empirical regression example.

-   **Simulation_continous.R:** Contains the functions to run different cross-validation splits and to shuffle the data for continious outcomes (used in EmpiricalExample_continious.R)

    -   final_results_df_5_b.csv: Results Binary Example, "complete" (N = 241, T = 50)

    -   final_results_df_4_b.csv: Results Binary Example, N = 100, T = 60

    -   final_results_df_3_b.csv: Results Binary Example, N = 30, T = 100

    -   final_results_df_2_b.csv: Results Binary Example, N = 20, T = 150

    -   final_results_df_5.csv: Results Regression Example, "complete" (N = 241, T = 50)

    -   final_results_df_4.csv: Results Regression Example, N = 100, T = 60

    -   final_results_df_3.csv: Results Regression Example, N = 30, T = 100

    -   final_results_df_2.csv: Results Regression Example, N = 20, T = 150

Figures:

-   **Figures/**: Directory containing all figures used in the paper.
