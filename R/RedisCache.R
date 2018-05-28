#' @export
cacheEnv <- new.env()
assign('cacheQ', TRUE, envir = cacheEnv)
assign('use', 'inmemory', envir = cacheEnv)
assign('cacheValues', list(), envir = cacheEnv)
assign('maxInMemoryEntries', 1e3, envir = cacheEnv)

#' @export
startRedis <- function(host = 'localhost', port = 6379) {
  rredis::redisConnect(host = host, port = port, nodelay = FALSE)
}

#' @export
setCache <- function(cacheQ = TRUE, use = 'redis') {
  assign('cacheQ', cacheQ, envir = cacheEnv)
  assign('use', use, envir = cacheEnv)
}

#' @export
cache <- function(f, ...) {
  funName <- as.character(substitute(f))
  if(!get('cacheQ', envir = cacheEnv)) {
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
    numArgs <- Filter(is.numeric, arguments)
    key <- do.call(paste, c(list(), mapply(function(v, n) paste0(n, '=', v), numArgs, names(numArgs))))
    key <- paste0(key, '-', funName, '-', extraKey, collapse = '')

    if(get('use', envir = cacheEnv) == 'redis') {
      startRedis()
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
      rredis::redisClose()
      val
    } else {
      # use in memory cache
      if(is.null(cacheValues[[key]])) {
        cacheValues[[key]] <- do.call(f, arguments)
        # if cached values are too big remove the older entries
        if(length(cacheValues) > get('maxInMemoryEntries', envir = cacheEnv))
          cacheValues[[1]] <- NULL
      }

      assign('cacheValues', cacheValues, envir = cacheEnv)
      cacheValues[[key]]
    }
  }
}

#' Clears all the cached values
#' @export
clearFunCache <- function(funName) {
  startRedis()
  if(is.function(funName))
    funName <- as.character(substitute(funName))

  s <- paste0('fset-', funName)
  keys <- rredis::redisSMembers(s)
  for (k in keys)
    rredis::redisDelete(k)

  rredis::redisDelete(s)  # remove all the keys
  rredis::redisClose()
}

#' Clears all the cached values
#' @export
clearCache <- function() {
  startRedis()
  rredis::redisFlushAll()
  rredis::redisClose()
}
