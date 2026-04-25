# ==============================================================================
# Script: generate_mock_data.R
# Purpose: Generates synthetic datasets to ensure code reproducibility 
#          while maintaining data privacy/ethics standards.
# ==============================================================================

set.seed(42) # For reproducibility
n_total <- 200

# 1. Generate Mock Demographic Data
# Matching variables: Id.number, Start.of.flwp, Gender, Age, Birthday, 
# Date.of.Death, Date.of.Migration, Exposure
mock_demographics <- data.frame(
  id_number = 1:n_total,
  Start.of.flwp = sample(20000101:20181231, n_total, replace = TRUE),
  Gender = sample(1:2, n_total, replace = TRUE),
  Age = sample(30:85, n_total, replace = TRUE),
  Birthday = as.Date('1940-01-01') + sample(0:15000, n_total, replace = TRUE),
  Date.of.Death = sample(c("20191010", "NA"), n_total, replace = TRUE, prob = c(0.05, 0.95)),
  Date.of.Migration = sample(c("20150505", "NA"), n_total, replace = TRUE, prob = c(0.02, 0.98)),
  Exposure = sample(c("MCL", "control"), n_total, replace = TRUE, prob = c(0.2, 0.8))
)

# 2. Generate Mock Register Data
# Matching variables: id.number, ICD.codes, Date.of.ICD.codes
# We include some codes starting with C18, C19, C20 to simulate CRC events.
icd_pool <- c("C189", "C199", "C209", "C349", "I109", "E119", "Z511")
n_register_entries <- 150

mock_register <- data.frame(
  id.number = sample(1:n_total, n_register_entries, replace = TRUE),
  ICD.codes = sample(icd_pool, n_register_entries, replace = TRUE),
  Date.of.ICD.codes = as.Date('2010-01-01') + sample(0:3000, n_register_entries, replace = TRUE)
)

# 3. Export to match original file formats (Tab-separated)
write.table(mock_demographics, "Demographic_data.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(mock_register, "Register_data.txt", sep="\t", row.names=FALSE, quote=FALSE)

cat("Successfully created mock datasets: Demographic_data.txt and Register_data.txt\n")
