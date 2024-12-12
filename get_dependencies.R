##############################################
# File for getting dependencies

# December 10, 2024
##############################################


if (!requireNamespace("here", quietly = TRUE)) {
    install.packages("here", dependencies = TRUE, repos = "http://cran.us.r-project.org")
}

# set working directory to current/source directory
library(here)
setwd(here::here())

r_files <- list.files(path = getwd(), pattern = "\\.R$", recursive = TRUE, full.names = TRUE)

# Function to extract package names from `library()` calls
extract_packages <- function(file) {

    lines <- readLines(file, warn = FALSE)


    matches <- regmatches(lines, gregexpr("(?<=library\\()[-\\.\\w]+(?=\\))", lines, perl = TRUE))


    unique(unlist(matches))
}

all_packages <- sort(unique(unlist(lapply(r_files, extract_packages))))

# Print required packages
cat("Required packages: \n")
print(all_packages)
cat("\n")

# save list of required packages to dependencies.txt
write.table(all_packages, file = "dependencies.txt", sep = "\n", row.names = FALSE, col.names = FALSE)


