# conda activate TCGA_env

library(TCGAbiolinks)
library(dplyr)

# Get the full list of GDC projects
projects <- getGDCprojects()

# Filter only those belonging to the TARGET program
target_projects <- projects %>%
  filter(grepl("TARGET", project_id))

# Display TARGET projects
print(target_projects[, c("project_id", "name", "disease_type")])
