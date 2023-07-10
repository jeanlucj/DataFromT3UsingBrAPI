library(httr)

brapi = list(
  server = NULL,
  version = "v2",
  get = function(call, ...) {
    brapi$request("GET", call, ...)
  },
  post = function(call, ...) {
    brapi$request("POST", call, ...)
  },
  put = function(call, ...) {
    brapi$request("PUT", call, ...)
  },
  request = function(method, call, ..., query=list(), body=list(), page=0, pageSize=10, token=NULL) {
    
    # Check to make sure the server is defined
    # Add https protocol if not already provided
    # remove trailing slash
    if ( is.null(brapi$server) ) stop("The BrAPI server is not defined!\nAn example on how to define the BrAPI server:\nbrapi$server = 'wheat.triticeaetoolbox.org'")
    if ( !startsWith(brapi$server, "http") ) brapi$server = paste0("https://", brapi$server)
    brapi$server = sub("/+$", "", brapi$server)
    
    # Check to make sure the call is defined
    # remove leading slash
    if ( !hasArg(call) ) stop("The BrAPI call is required!")
    call = sub("^/+", "", call)
    
    # Build the full URL
    url = paste(brapi$server, "brapi", brapi$version, call, sep="/")
    
    # Add page parameters to the query list, for GET requests
    if ( method == "GET" ) {
      query$page = page
      query$pageSize = pageSize
    }
    
    # Add token as Authorization header, if provided
    config = list()
    if ( !is.null(token) ) {
      config = add_headers(Authorization = paste("Bearer", token, sep=" "))
    }
    
    # Make the Request
    resp = VERB(
      method, url, config,
      query = query,
      body = body,
      encode = "json"
    )
    
    # Print Response Info
    cat(
      sprintf("Response [%s]", resp$url),
      sprintf("  %s", http_status(resp)$message),
      sprintf("  Content Type: %s", http_type(resp)),
      sep = "\n"
    )
    warn_for_status(resp)
    content = content(resp)
    
    # Check for error message in the metadata
    if ( "metadata" %in% names(content) && "status" %in% names(content$metadata) ) {
      for ( status in content$metadata$status ) {
        if ( status$messageType == "ERROR" ) {
          warning(status$message)
        }
      }
    }
    
    # Parse the content
    return(list(
      response = resp,
      status = http_status(resp),
      content = content
    ))
  }
)