#' Find the path to the compiled models
#'
#' @return a character string
#' @export
#'
#' @examples find_compiled_model()
find_compiled_model <- function(){
  dyn_lib <- system.file(paste("libs/seabirdeb", .Platform$dynlib.ext, sep = ""), package="seabirdeb")
  if (dyn_lib == ""){
    dyn_lib <- system.file(paste("libs/", .Platform$r_arch, "/seabirdeb", .Platform$dynlib.ext, sep = ""), package="seabirdeb")
  }
  return(dyn_lib)
}
