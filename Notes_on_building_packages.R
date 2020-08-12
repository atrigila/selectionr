
#Create datasets
sample.alignment.1 <- readLines("/Users/Usuario/Desktop/AAK1_final_mask_align_NT.aln")

sample.alignment.2 <- readLines("/Users/Usuario/Desktop/ACSL4_final_unmask_align_NT.aln")

use_data(sample.alignment.2)


# Generate package documentation
devtools::document("/Users/Usuario/selectionr/")

# Examine the contents of the man directory
dir("/Users/Usuario/selectionr/man")

# View the documentation for the fast.classify function
help("fast.classify")

# Check your package
check("/Users/Usuario/selectionr/")

# Add dependencies to the DESCRIPTION file


# Build the package
build("datasummary")

