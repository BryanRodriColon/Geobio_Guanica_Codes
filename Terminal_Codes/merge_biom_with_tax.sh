#!/bin/bash
# Change to the working directory
cd /kuhpc/work/sturm/b324r112/PROJECTS/Puerto_Rico/Amplicon/Sci_Rep_data_2019/07_export_data/biom_manipulation

# (Optional) Ensure your local biom-format executables are in PATH:
export PATH=$WORK/.local/bin:$PATH



# 1. Convert the BIOM file to a TSV with taxonomy using biom-format.
biom convert -i 050525_feature-table.biom -o feature-table.tsv --to-tsv --table-type="OTU table" --header-key taxonomy



# 2. Run the R processing code.
module load R

Rscript - <<'EOF'

# Set working directory to ensure outputs are saved in the correct folder.
setwd("/kuhpc/work/sturm/b324r112/PROJECTS/Puerto_Rico/Amplicon/Sci_Rep_data_2019/07_export_data/biom_manipulation")
cat("Current working directory: ", getwd(), "\n")

# Load required packages
library(stringr)
library(dplyr)

# Load the feature table TSV; skip the comment header line.
feature_table <- read.delim("feature-table.tsv", header = TRUE, skip = 1, stringsAsFactors = FALSE, comment.char = "")

# Remove a leading '#' from the first column header if present.
colnames(feature_table)[1] <- sub("^#", "", colnames(feature_table)[1])


# Load the taxonomy file
taxonomy <- read.delim("050525_taxonomy.tsv", header = TRUE, stringsAsFactors = FALSE)

### **Remove "Eukaryota",  NA Entries and taxonomy with less than %70 confidence**
taxonomy <- taxonomy %>%
  filter(Confidence >= 0.70 & !is.na(Taxon) & !grepl("d__Eukaryota", Taxon))

# Define a function to split a taxonomy string and pad if needed.
split_and_pad_tax <- function(tax_str, expected = 7) {
  parts <- unlist(strsplit(tax_str, ";"))
  parts <- str_trim(parts)  # remove extra whitespace
  if (length(parts) < expected) {
    parts <- c(parts, rep(NA, expected - length(parts)))
  }
  names(parts) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  return(parts)
}

# Define a helper to replace blanks or entries that are only the prefix.
fix_taxon_value <- function(x, prefix) {
  x_lower <- tolower(x)
  if (is.na(x) || x == "" || x == paste0(prefix, "__") ||
      grepl(paste0("^", prefix, "__$"), x)) {
    return(paste0(prefix, "__Unclassified"))
  } else {
    return(x)
  }
}

# Process each taxonomy string from the taxonomy file.
fixed_tax_list <- lapply(taxonomy$Taxon, function(x) {
  parts <- split_and_pad_tax(x, expected = 7)
  parts["Domain"]  <- fix_taxon_value(parts["Domain"],  "d")
  parts["Phylum"]  <- fix_taxon_value(parts["Phylum"],  "p")
  parts["Class"]   <- fix_taxon_value(parts["Class"],   "c")
  parts["Order"]   <- fix_taxon_value(parts["Order"],   "o")
  parts["Family"]  <- fix_taxon_value(parts["Family"],  "f")
  parts["Genus"]   <- fix_taxon_value(parts["Genus"],   "g")
  parts["Species"] <- fix_taxon_value(parts["Species"], "s")
  return(parts)
})

# Remove Eukaryota and NA domains **before merging**
taxonomy <- taxonomy[!(taxonomy$Taxon == "Eukaryota" | is.na(taxonomy$Taxon)), ]



# Recombine the fixed taxonomic levels into a single string and Add the fixed taxonomy as a new column to the data frame.
taxonomy$new_taxonomy <- sapply(fixed_tax_list, function(x) paste(x, collapse = ";"))

taxonomy <- data.frame(Feature.ID = taxonomy$Feature.ID,
                         taxonomy = taxonomy$new_taxonomy,
                         stringsAsFactors = FALSE)



#Rename the ID columns for consistency.
colnames(taxonomy)[colnames(taxonomy) == "Feature.ID"] <- "X.OTU.ID"



# Filter the feature table to include only ASVs that passed the taxonomy filtering
feature_table <- feature_table[feature_table[[1]] %in% taxonomy$X.OTU.ID, ]

#Merge the feature table and the updated taxonomy by matching the ASV/taxon IDs.
feature_table_w_taxonomy <- merge(feature_table, taxonomy, 
                      by.x = colnames(feature_table)[1], by.y = "X.OTU.ID",
                      all.x = TRUE, sort = FALSE)

# If merge creates duplicate taxonomy columns, keep the new one.
if("taxonomy.y" %in% colnames(feature_table_w_taxonomy)){
  feature_table_w_taxonomy$taxonomy <- feature_table_w_taxonomy$taxonomy.y
  feature_table_w_taxonomy$taxonomy.x <- NULL
  feature_table_w_taxonomy$taxonomy.y <- NULL
}

# Ensure that the first column header is exactly "#OTU ID".
colnames(feature_table_w_taxonomy)[1] <- "#OTU ID"


# Write the comment line as the first row.
cat("# Constructed from biom file\n", file = "050525_feature_table_w_taxonomy.tsv")


# Append the merged table (with header) to the file.
write.table(feature_table_w_taxonomy, file = "050525_feature_table_w_taxonomy.tsv", sep = "\t", 
            row.names = FALSE, quote = FALSE, append = TRUE)


EOF

# 3. Convert the updated TSV back into a BIOM file.
biom convert -i 050525_feature_table_w_taxonomy.tsv -o 050525_feature_table_w_taxonomy.biom --to-json --table-type="OTU table" --process-obs-metadata taxonomy

echo "Pipeline complete: Updated BIOM file created as 050525_feature_table_w_taxonomy.biom"
