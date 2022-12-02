# Taken from: https://www.r-bloggers.com/finding-my-dropbox-in-r/

get_dropbox <- function () {
  
  # install.packages("rjson")
  if (Sys.info()['sysname'] == 'Darwin') {
    info <- RJSONIO::fromJSON(
      file.path(path.expand("~"),'.dropbox','info.json'))
  }
  if (Sys.info()['sysname'] == 'Windows') {
    info <- RJSONIO::fromJSON(
      if (file.exists(file.path(Sys.getenv('APPDATA'), 'Dropbox','info.json'))) {
        file.path(Sys.getenv('APPDATA'), 'Dropbox', 'info.json')
      } else {
        file.path(Sys.getenv('LOCALAPPDATA'),'Dropbox','info.json')
      }
    )
  }
  
  dir <- info$personal$path
  return(dir)
}
