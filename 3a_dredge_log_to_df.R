# Author: Alex Morin

## This script provides a function to turn a text log of dredge outputs into a
## dataframe that tracks the warning status for each fit model.
## Usage:: dredge_log_to_df("log.txt")

## TODO: Makes strong assumptions about format of the output. May fail if the
## text is 'contaminated' with non-dredge warnings calls. Could write rule to
## check first and last ids before string matching. 
## TODO: Consider redesign of warnings_to_df(), which could use more segmented
## functions with assertions for each possible line.
## TODO: Functionality to isolate variables for each model call


library(stringr)


get_model_ids <- function(log_lines) {
  # Assumes dredge outputs distinct model calls starting with a digit identifier. 
  # Returns a string vector of these isolated digits.
  string <- str_extract(log_lines, "^[:digit:]+")
  string[!is.na(string)] 
}


init_warning_df <- function(ids) {
  # Initialize an n x 4 dataframe where n is the number of models
  
  data.frame(
    Model_id = ids,
    Has_warning = rep(FALSE, length(ids)),
    Probability_warning = rep(FALSE, length(ids)),
    Converge_warning = rep(FALSE, length(ids)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}


detect_warning <- function(line) {
  # Helper to detect if any warning was found and returns TRUE or FALSE
  str_detect(line, "^Warning:")
}


is_prob_warning <- function(line) {
  # Helper to detect if a probability fit warning was raised
  str_detect(line, ".*fitted probabilities numerically 0 or 1 occurred.*")
}


is_converge_warning <- function(line) {
  # Helper to detect if a convergence error was raised
  str_detect(line, ".*algorithm did not converge.*")
  
}


warnings_to_df <- function(log_df, log_lines) {
  # Iterates over each line of log lines and flags warnings for each model call
  # if a warning was generated.
  
  id <- ""
  
  for (line in log_lines) {
    
    if (str_detect(line, "^ ")) {   # skip if carry over from model call
      next
    }
    
    if (str_detect(line, "^[:digit:]+")) {
      id <- str_extract(line, "^[:digit:]+")
      
    } else if (detect_warning(line)) {
      log_df[log_df$Model_id == id, "Has_warning"] <- TRUE
      
      if (is_converge_warning(line)) {
        log_df[log_df$Model_id == id, "Converge_warning"] <- TRUE
      }
      
      if (is_prob_warning(line)) {
        log_df[log_df$Model_id == id, "Probability_warning"] <- TRUE
      }
    }
  }
  
  return (log_df)
}


dredge_log_to_df <- function(log_file) {
  # Main function to load a dredge log file and return an n x 4 dataframe,
  # where n is the amount of discrete model calls
  
  log_lines <- readLines(log_file)
  ids <- get_model_ids(log_lines)
  warnings_to_df(log_df = init_warning_df(ids), log_lines)
}
