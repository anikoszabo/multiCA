## functions to run nuweb on all files in the package

nuweb_files <- function(pkg){
  path_nuweb <- file.path(pkg$path, "nuweb")
  files <- dir(path_nuweb, pattern = "\\.[Ww]$", full.names = TRUE)
  names(files) <- basename(files)
  withr::with_collate("C", sort(files))
}

run_nuweb <- function(file, path){
  call <- paste("nuweb -p", path, file, collapse=" ")
  system(call)
}

nuweb <- function(pkg = ".", path=NULL){
  pkg <- as.package(pkg)
  if (is.null(path)) {
    path <- file.path(pkg$path, "nuweb")
  }
  files <- nuweb_files(pkg)
  if (length(files) == 0) 
    return()
  cat(paste0("Running ", length(files), " nuweb files in ", pkg$package, "\n"))
  lapply(files, run_nuweb, path=path)
}