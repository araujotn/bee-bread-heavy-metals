# Heavy metal contamination in the larval food of a solitary bee species 🌸🐝

## Description
This project investigates heavy metal contamination in *Xylocopa frontalis* bee bread, a solitary bees species that occurs at the Brazilian savanna (Cerrado). The goal is to analyze larval food samples, identify heavy metal concentration, and investigaste the difference between the yellow passion fruit crops and proteted savanna samples.


## Project Structure
- **/data**: Datasets used for analysis.
- **/src**: Scripts and dataset for data analysis.
- **/results**: Generated figures.

## How to Use

1. Download the script file and dataset in the `/src` folder and save them in a project folder named `bee-bread-heavy-metals`:
   ```
   bee-bread-heavy-metals/
   ├── data/
   ├── src/
   │   ├── heavy_metals_analysis.R
   │   └── heavy_metals_data.csv
   └── results/
   ```

2. Open R or RStudio and set the working directory to the project folder:
   ```r
   setwd("path/to/bee-bread-heavy-metals")
   ```

3. Install the required R packages if they are not already installed. You can run the following code in R:
   ```r
   packages <- c("robustbase", "car", "sdmTMB", "glmmTMB", "dplyr", "DHARMa", "ggplot2", "ggpubr", "gridExtra", "vegan")
   installed <- packages %in% rownames(installed.packages())
   if(any(!installed)) install.packages(packages[!installed])
   ```

4. Run the analysis scripts in the `/src` folder.


## Statistical Methods
The analytical approach was designed to first identify overall differences in the heavy metal profiles between site types (yellow passion fruit crop and Brazilian savanna areas) and then to model the behavior of each individual metal.

1.  **Overall Profile Comparison**:
    * A **Permutational Multivariate Analysis of Variance (PERMANOVA)** was performed to test for significant differences in the complete heavy metal profile between the crop and savanna sites.

2.  **Individual Metal Analysis**:
    * **Generalized Linear Models (GLMs)** were fitted to investigate the effect of site type on each metal's concentration. To account for the statistical properties of each metal's dataset, different model families and link functions were used:
        * **Robust GLM (Gamma family, log link)**: Applied to `Al`, `Cu`, `Fe`, `Sn`, and `Zn`. This approach is robust to outliers and suitable for the positive, right-skewed distribution of concentration data.
        * **GLM (Gaussian family)**: Used for `Cr`, as its model residuals were consistent with a normal distribution.
        * **GLM (Binomial family, logit link)**: Used for `Ni`, which was treated as a binary (presence/absence) variable due to its absence in all crop site samples.
        * **Two-Part Hurdle Model (Delta-Gamma)**: Applied to `Cd` and `Pb`, which had zero-inflated distributions. This model separately assesses:
            1.  The probability of the metal's presence (Binomial model, logit link).
            2.  The concentration level, conditional on its presence (Gamma model, log link).


## Contributing
Contributions are welcome! Please open an issue or submit a pull request via e-mail (thayane.n.a@gmail.com).

## Funding
This study was funded by scholarships from Conselho Nacional de Desenvolvimento Científico e Tecnológico (CNPq) and Coordenação de Aperfeiçoamento de Pessoal de Nível Superior (CAPES). Also received financial support from project Programa de Pesquisas Ecológicas de Longa Duração - Triângulo Mineiro e Sudeste de Goiás (PELD/TMSG) (441225/2016-0 and 441142/2020-6); and by Fundação de Amparo à Pesquisa do Estado de Minas Gerais (FAPEMIG) (APQ 04815-17). 

## References
Anderson, SC, Ward, EJ, English, PA and Barnett, LAK. (2022) sdmTMB: an R package for fast, flexible, and user-friendly generalized linear mixed effects models with spatial and spatiotemporal random fields. bioRxiv. DOI: 10.1101/2022.03.24.485545

Araújo, TN, Junqueira, CN, Castro-Melo, ALS, Rocha-Filho, LC, Santos, DQ and Augusto, SC. (2025) A scientific note on the heavy metal contamination in the larval food of Xylocopa frontalis (Apidae) at the Brazilian savanna. Apidologie. DOI: 

Auguie, B. (2017) gridExtra: Miscellaneous Functions for "Grid" Graphics. R package version 2.3. URL: https://CRAN.R-project.org/package=gridExtra

Brooks, ME, Kristensen, K, van Benthem, KJ, Magnusson, A, Berg, CW, Nielsen, A, Skaug, HJ, Maechler, M and Bolker, B.M. (2017) glmmTMB Balances Speed and Flexibility Among Packages for Zero-inflated Generalized Linear Mixed Modeling. The R Journal, 9(2), 378-400. DOI: <u>10.32614/RJ-2017-066</u>

Fox, J and Weisberg, S. (2019) An {R} Companion to Applied Regression, Third Edition. Thousand Oaks CA: Sage. URL: https://socialsciences.mcmaster.ca/jfox/Books/Companion/

Hartig, F. (2022) DHARMa: Residual Diagnostics for Hierarchical (Multi-Level / Mixed) Regression Models. R package version 0.4.6. URL: https://CRAN.R-project.org/package=DHARMa

Kassambara, A. (2020) ggpubr: 'ggplot2' Based Publication Ready Plots. R package version 0.4.0. URL: https://CRAN.R-project.org/package=ggpubr

Maechler, M, Rousseeuw, P, Croux, C, Todorov, V, Ruckstuhl, A, Salibian-Barrera, M, Verbeke, T, Koller, M, Conceicao, ELT and di Palma, MA. (2023) robustbase: Basic Robust Statistics. R package version 0.95-1. URL: http://CRAN.R-project.org/package=robustbase

Oksanen, J, Blanchet, FG, Friendly, M, Kindt, R, Legendre, P, McGlinn, D, Minchin, PR, O'Hara, RB, Simpson, GL, Solymos, P, Stevens, MHH, Szoecs, E and Wagner, H. (2020) vegan: Community Ecology Package. R package version 2.5-7. URL: https://CRAN.R-project.org/package=vegan

Wickham, H. (2016) ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York.

Wickham, H, François, R, Henry, L and Müller, K. (2021) dplyr: A Grammar of Data Manipulation. R package version 1.0.7. URL: https://CRAN.R-project.org/package=dplyr