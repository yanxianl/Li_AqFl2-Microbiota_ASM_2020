library(Maaslin2)

input_data <- system.file(
  'extdata','example1_features.txt', package="Maaslin2")
input_metadata <-system.file(
  'extdata','example1_metadata.txt', package="Maaslin2")
fit_data <- Maaslin2(
  input_data, input_metadata,'maaslin2_example_output')
