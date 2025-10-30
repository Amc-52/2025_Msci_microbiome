# 2025_Msci_microbiome
Repository for scripts and workflows relating to my masters project: Exploring gut microbiome-aGVHD interactions using metabolic modelling in paediatric leukaemia.

## Project Abstract
Acute graft-versus-host disease (aGVHD) is a major and potentially fatal complication of allogeneic hematopoietic stem cell transplantation HSCT, a curative treatment for childhood leukaemia. A wealth of research implicating gut microbiome community composition in the pathogenesis and severity of aGVHD exists, however the complexities of this relationship have not yet been resolved. Constraint-based modelling was performed on publicly available 16S data from a cohort of paediatric leukaemia patients undergoing allogeneic hematopoietic stem cell transplantation (HSCT) to produce personalised metabolic models (PMMs) of patient gut microbiota. Random Forest models were created to identify important bacteria, metabolic reactions and subsystems in the onset of aGVHD. Cross feeding, growth and metabolic reactions were then simulated using agent-based modelling to gain insight into cell-level dynamics of the gut microbiomes of patients with aGVHD. This workflow shows the potential of using modern techniques on older data, and argues for the implementation of a knowledge base where these PMMs could be stored and analysed on a larger scale.

## Methods Overview
Briefly- raw 16S rRNA reads were accessed from NCBI and preprocessed. These reads were imported into RStudio where further processing was performed as well as taxonomic identification to the genus level. Mathmatical models of bacterial metabolism were generated with the COBRA toolbox and AGORA2 reconstructions (Heirendt, L. et al. 2019; Heinken, A. et al. 2023), and back in RStudio, agent based modelling with the package BacArena (Bauer E, Zimmermann J 2025) was used to investigate bacterial community dynamics. Various other analyses such as Random Forest Modelling was also performed. 

## References
Heirendt, L. et al. (2019) ‘Creation and analysis of biochemical constraint-based models using the COBRA Toolbox v.3.0’, Nature Protocols, 14(3), pp. 639–702. Available at: https://doi.org/10.1038/s41596-018-0098-2.

Bauer E, Zimmermann J (2025). _BacArena: Modeling Framework for Cellular Communities in their Environments_. R package version 1.8.2, commit 7e6f4b5b412cedc1e6b5f6965e790813f8c90aca, <https://github.com/euba/bacarena>.

Heinken, A. et al. (2023) ‘Genome-scale metabolic reconstruction of 7,302 human microorganisms for personalized medicine’, Nature Biotechnology, 41(9), pp. 1320–1331. Available at: https://doi.org/10.1038/s41587-022-01628-0.

