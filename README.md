# Heavy metals contamination in bees larval food

## Description
This project investigates heavy metal contamination in Xylocopa frontalis bee bread, a solitary bees species that occurs at the Brazilian savanna (Cerrado). The goal is to analyze larval food samples, identify heavy metal concentration, and investigaste the difference between the yellow passion fruit crops and proteted savanna samples.


## Project Structure
- **/data**: Datasets used for analysis.
- **/src**: Scripts and source code for data analysis.
- **/results**: Analysis results and generated figures.

## How to Use
git commit -m "first commit"
1. Clone this repository:
   ```sh
   git clone https://github.com/araujotn/bee-bread-heavy-metals.git
   ```

2. Open R or RStudio and set the working directory to the project folder:
   ```r
   setwd("path/to/Heavy-metal-contamination-in-bees-larval-food")
   ```

3. Install the required R packages (example):
   ```r
   install.packages(c("robustbase", "car", "sdmTMB", "glmmTMB", "dplyr", "DHARMa", "ggplot2", "ggpubr", "gridExtra", "vegan"))
   ```

4. Run the analysis scripts in the `/src` folder:
   ```r
   source("src/heavy-metals-analysis.R")
   ```

## Contributing
Contributions are welcome! Please open an issue or submit a pull request.

## Funding
This study was funded by scholarships from Conselho Nacional de Desenvolvimento Científico e Tecnológico (CNPq) and Coordenação de Aperfeiçoamento de Pessoal de Nível Superior (CAPES). Also received financial support from project Programa de Pesquisas Ecológicas de Longa Duração - Triângulo Mineiro e Sudeste de Goiás (PELD/TMSG) (441225/2016-0 and 441142/2020-6); and by Fundação de Amparo à Pesquisa do Estado de Minas Gerais (FAPEMIG) (APQ 04815-17). 
