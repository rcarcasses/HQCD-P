#' @export
cacheQ <- TRUE

#' @export
startRedis <- function(host = 'localhost', port = 6379) {
  cacheQ <-TRUE
  cat('connecting to redis...')
  rredis::redisConnect(host = host, port = port, nodelay = FALSE)
  cat('DONE\n')
}

#' @export
stopRedis <- function() {
  system('docker stop redis')
}

#' @export
cache <- function(f, ...) {
  funName <- as.character(substitute(f))
  if(!cacheQ) {
    flog.debug(paste('[RedisCache] function ', funName, ' won\'t be cached, caching is disabled'))
    return(f)
  }

  extraKey <- ''
  if(length(list(...)) > 0) {
    extraKey <- paste0(list(...), collapse = ',')
    if(extraKey != '')
      flog.debug(paste('[RedisCache] Using extra key *', extraKey, '* for labeling cached values for', funName))
  }

  function(...) {
    # beware of potential issues with partial argument matching
    # see: https://stackoverflow.com/questions/15264994/prevent-partial-argument-matching
    arguments <- list(...)
    key <- paste0(format(Filter(is.numeric, arguments), digits = 16), collapse = ',')
    key <- paste0(key, '-', funName, '-', extraKey, collapse = '')

    # cat('key is ', key, '\n')
    val <- rredis::redisGet(key)  # check if it has been already computed
    if(!is.null(val))
      return(val)

    #otherwise compute the value, store and return it
    val <- do.call(f, arguments)
    if(!is.null(val)){
      set <- paste0('fset-', funName)
      rredis::redisSAdd(set, key)    # all the functions keys are stored in a set
      rredis::redisSet(key, val)
    }

    val
  }
}
#' Clears all the cached values
#' @export
clearFunCache <- function(funName) {
  if(is.function(funName))
    funName <- as.character(substitute(funName))

  s <- paste0('fset-', funName)
  keys <- rredis::redisSMembers(s)
  for (k in keys)
    rredis::redisDelete(k)

  rredis::redisDelete(s)  # remove all the keys
}

#' Clears all the cached values
#' @export
clearCache <- function() {
  rredis::redisFlushAll();
}
