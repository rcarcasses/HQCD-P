# this function helps to load external solutions and data
#' @export
#' @param entity It can be a solution for A,Phi or a normalizable mode, etc.
#' @param baseDir The base directory from which the files will be look for
#' @param afterLoad A function that can be used to further format the imported data
loadExt <- function(entity, baseDir = '', afterLoad = function(x) x) {
  # check if there is a rds file in the working directory with the name of the
  # entity, if there is one, load it and return the object inside
  handlers <- list(
    list(extension = '.rds', handler = readRDS),
    list(extension = '.csv', handler = read.csv)
    )
  # try all handlers, remove the NULL ones
  obj <- Filter(Negate(is.null), lapply(handlers, function(h) loadExtType(entity, baseDir, h$extension, h$handler)))
  # get the first one and apply afterLoad, then return the result
  if(length(obj) > 0) {
    flog.trace('Using %s data from file %s', entity, obj[[1]]$fileName)
    afterLoad(obj[[1]]$data)
  }
  else
    NULL
}

# check if exist a file with name 'entity.extension', if it does,
# then call handler to load its data, else returns NULL
loadExtType <- function(entity, baseDir, extension, handler, ...) {
  fileName <- paste0(baseDir, entity, extension)
  if(file.exists(fileName)) {
    flog.trace('Loading %s data from file %s', entity, fileName)
    list(data = handler(fileName, ...), fileName = fileName)
  }
  else
    NULL
}
