
#===================================================================
# Matched Cohort Analysis: MCL and Colorectal Cancer as Second Primary Malignancy
#===================================================================
# Author: Dr. Mohammad Haroon Khan
# Date: March 27, 2026
# R Version: 4.5.1
# Purpose: Assess CRC risk in MCL patients vs. population controls
# ===================================================================

# --- Load Required Libraries ---
library(dplyr)
library(lubridate)
library(survival)

#================================================================ 
# STEP 1: PREPROCESSING & DATA CLEANING 
#================================================================
# 1.1 Read the two raw files exactly as provided
demo <- read.table("Demographic_data.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE, 
                   na.strings = c("NA", "N/A", ""), quote = "", fill = TRUE, comment.char = "")
register <- read.table("Register_data.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, 
                       na.strings = c("NA", "N/A", ""), quote = "", fill = TRUE)
View(demo)
View(register)

# 1.2 Quick raw inspection
cat("Raw Demographic rows:", nrow(demo), "\n")
cat("Raw Register rows:", nrow(register), "\n")
cat("Exposure levels before cleaning:", paste(unique(demo$Exposure), collapse = ", "), "\n")

#================================================================
# STEP 2: DATA CLEANING
# ================================================================
# 2.1 Standardise column names (R is case-sensitive)
names(demo)[names(demo) == "id_number"] <- "id_number"
names(demo)[names(demo) == "Start.of.flwp"] <- "start_flwp"
names(demo)[names(demo) == "Date.of.Death"] <- "death_date"
names(demo)[names(demo) == "Date.of.Migration"] <- "mig_date"
names(register)[names(register) == "id.number"] <- "id_number"

# 2.2 Convert all date columns (YYYYMMDD or YYYY-MM-DD format) using base R
demo <- demo %>%
    mutate(
        start_flwp = as.Date(as.character(start_flwp), format = "%Y%m%d"),
        Birthday   = as.Date(Birthday, format = "%Y-%m-%d"),
        death_date = as.Date(death_date, format = "%Y%m%d"),
        mig_date   = as.Date(mig_date, format = "%Y%m%d")
    )
register <- register %>%
    mutate(Date.of.ICD.codes = as.Date(Date.of.ICD.codes, format = "%Y-%m-%d"))

# 2.3 Clean Exposure variable (standardise spelling)
demo <- demo %>%
    mutate(Exposure = trimws(tolower(Exposure))) %>%
    mutate(Exposure = case_when(
        Exposure == "mcl" ~ "MCL",
        Exposure == "control" ~ "control",
        TRUE ~ NA_character_
    ))

# 2.4 Remove any rows with missing critical identifiers
demo <- demo %>% filter(!is.na(id_number) & !is.na(start_flwp) & !is.na(Exposure))
register <- register %>% filter(!is.na(id_number) & !is.na(Date.of.ICD.codes))

# 2.5 Identify first CRC event (C18, C19 or C20) per person AFTER start of follow-up 
# (we keep only the earliest valid date)
crc_events <- register %>%
    mutate(ICD_list = strsplit(ICD.codes, "\\s+")) %>%   
    # better splitting on any whitespace
    rowwise() %>%
    mutate(
        is_CRC = any(sapply(ICD_list, function(x) {
            code <- trimws(x)
            # Match C18, C180–C189, C19, C199, C20, C209 etc.
            grepl("^C18[0-9]?$|^C19[0-9]?$|^C20[0-9]?$", code)
        }))
    ) %>%
    filter(is_CRC) %>%
    group_by(id_number) %>%
    summarise(
        crc_date = min(Date.of.ICD.codes, na.rm = TRUE),
        .groups = "drop"
    )

# 2.6 Merge again with the corrected CRC dates
demo <- left_join(demo, crc_events, by = "id_number")

# 2.7 Recalculate event and follow-up (same as before)
end_fu <- as.Date("2019-12-31")
demo <- demo %>%
    mutate(
        censor_date=pmin(death_date, mig_date, end_fu, na.rm = TRUE),
        event_date=ifelse(!is.na(crc_date) & crc_date <= censor_date, crc_date, NA),
        event_date=as.Date(event_date, origin = "1970-01-01"),
        time_days=as.numeric(censor_date - start_flwp)
    ) %>%
    filter(time_days > 0) %>%
    mutate(
        event = as.integer(!is.na(event_date)),
        time_years = time_days / 365.25
    )

# Check again
cat("\n=== AFTER FIX ===\n")
cat("Persons with CRC event:", sum(demo$event, na.rm = TRUE), "\n")
cat("MCL patients with CRC:", sum(demo$event[demo$Exposure == "MCL"], na.rm = TRUE), "\n")
cat("Controls with CRC:", sum(demo$event[demo$Exposure == "control"], na.rm = TRUE), "\n")

# ================================================================
# 3. FINAL ANALYSIS: Incidence Rates + Cox Models
# ================================================================
# 3.1. Summary statistics by group

summary_stats <- demo %>%
    group_by(Exposure) %>%
    summarise(n_persons = n(),
              person_years = round(sum(time_years, na.rm = TRUE), 1),
              n_events = sum(event, na.rm = TRUE),
              .groups = "drop"
    ) %>%
    mutate(
        ir_per1000 = round(n_events / person_years * 1000, 2),
        # Byar 95% CI for Poisson rate
        ir_lower = round((qchisq(0.025, 2 * n_events) / 2 / person_years) * 1000, 2),
        ir_upper = round((qchisq(0.975, 2 * n_events + 2) / 2 / person_years) * 1000, 2),
        ir_ci = paste0(ir_per1000, " (", ir_lower, "–", ir_upper, ")")
    )

print("=== SUMMARY STATISTICS ===")
print(summary_stats)

# 3.2. Cox models
cox_unadj <- coxph(Surv(time_years, event) ~ Exposure, data = demo)
cox_adj <- coxph(Surv(time_years, event) ~ Exposure + Age + as.factor(Gender), data = demo)

# 3.3. Safe extraction of HR and 95% CI (avoids subscript out of bounds)
# Unadjusted
unadj_summary <- summary(cox_unadj)
hr_unadj <- unadj_summary$conf.int[1, "exp(coef)"]      

# first coefficient after intercept
lower_unadj <- unadj_summary$conf.int[1, "lower .95"]
upper_unadj <- unadj_summary$conf.int[1, "upper .95"]

# Adjusted
adj_summary <- summary(cox_adj)
hr_adj <- adj_summary$conf.int[1, "exp(coef)"] # Exposure coefficient is first
lower_adj <- adj_summary$conf.int[1, "lower .95"]
upper_adj <- adj_summary$conf.int[1, "upper .95"]

# 4. Build the exact table requested in the Instructions.docx
final_table <- data.frame(
    Metric = c("Number of persons", "Number of person years", "Number of events",
               "Incidence rate per 1000 PY (95% CI)", "Unadjusted HR (95% CI)",
               "Adjusted* HR (95% CI)"),
    Controls = c(
        summary_stats$n_persons[summary_stats$Exposure == "control"],
        summary_stats$person_years[summary_stats$Exposure == "control"],
        summary_stats$n_events[summary_stats$Exposure == "control"],
        summary_stats$ir_ci[summary_stats$Exposure == "control"],
        "1 (Ref.)",
        "1 (Ref.)"
    ),
    `MCL patients` = c(
        summary_stats$n_persons[summary_stats$Exposure == "MCL"],
        summary_stats$person_years[summary_stats$Exposure == "MCL"],
        summary_stats$n_events[summary_stats$Exposure == "MCL"],
        summary_stats$ir_ci[summary_stats$Exposure == "MCL"],
        paste0(round(hr_unadj, 2), " (", round(lower_unadj, 2), "–", round(upper_unadj, 2), ")"),
        paste0(round(hr_adj, 2), " (", round(lower_adj, 2), "–", round(upper_adj, 2), ")")
    ),
    check.names = FALSE
)
print("=== FINAL RESULTS TABLE ===")
print(final_table, right = FALSE)

# 5. Save for my document
write.table(final_table, file = "Haroon_SPM.txt", sep = "\t", 
            row.names = FALSE, quote = FALSE)

cat("\nTable successfully saved as 'Haroon_SPM.txt'\n")
cat("Copy this into Word, format as a nice table, and add interpretation.\n")
