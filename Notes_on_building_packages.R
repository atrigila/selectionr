
#Create datasets
sample.alignment.1 <- readLines("/Users/Usuario/Desktop/AAK1_final_mask_align_NT.aln")

sample.alignment.2 <- readLines("/Users/Usuario/Desktop/ACSL4_final_unmask_align_NT.aln")

use_data(sample.alignment.2)

#Generate vignette
library(devtools)
use_vignette("selectionr_vignette", title = "selectionr")


# Generate package documentation
devtools::document("/Users/Usuario/selectionr/")

# Examine the contents of the man directory
dir("/Users/Usuario/selectionr/man")

# View the documentation for the fast.classify function
help("ensembl.orthologue.download")

# Check your package
check("/Users/Usuario/selectionr/")

# Add dependencies to the DESCRIPTION file



# Build the package
build("/Users/Usuario/selectionr/")


# Continuous integration
#You can run use_github() and use_travis() to set up your package for use with GitHub and Travis CI.
# Travis CI can be set up to run checks every time you update your code.

#Unit tests

# Set up the test framework
## You can set up a test framework in a package using the function use_testthat().

use_testthat("/Users/Usuario/selectionr/")

# Look at the contents of the new folder which has been created
dir("/Users/Usuario/selectionr/tests")

# Create a summary of the iris dataset using your data_summary() function
iris_summary <- data_summary(iris)

# Count how many rows are returned
summary_rows <- nrow(iris_summary)

# Use expect_equal to test that calling data_summary() on iris returns 4 rows
expect_equal(summary_rows, 4)

# Install package
setwd("..")


install("cats")





# Use context() and test_that() to group the tests below together
context("Test data_summary()")

# To run tests in packages you need to group these individual tests together. You do this using a function test_that(). You can use this to group together expectations that test specific functionality.

# You can use context() to collect these groups together. You usually have one context per file. An advantage of doing this is that it makes it easier to work out where failing tests are located.

test_that("data_summary() handles errors correctly", {

  # Create a vector
  my_vector <- 1:10

  # Use expect_error()
  expect_error(data_summary(my_vector))

  # Use expect_warning()
  expect_warning(data_summary(airquality, na.rm = FALSE))

})


# Run the tests on the datasummary package
# With your tests scripts saved in the package structure you can always easily re-run your tests using the test() function in devtools. This function looks for all tests located in the tests/testhat or inst/tests directory with filenames beginning with test- and ending in .R, and executes each of them. As with the other devtools functions, you supply the path to the package as the first argument to the test() function.
test("datasummary")




