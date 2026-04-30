# conda activate TCGA_env

library(TCGAbiolinks)
library(dplyr)

# Ottieni la lista completa dei progetti GDC
projects <- getGDCprojects()

# Filtra solo quelli del programma TARGET
target_projects <- projects %>%
  filter(grepl("TARGET", project_id))

# Visualizza i progetti TARGET
print(target_projects[, c("project_id", "name", "disease_type")])
