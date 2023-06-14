

update_record <- function(record_path, text) {
  write(text, 
        file = record_path, 
        append = TRUE)
}