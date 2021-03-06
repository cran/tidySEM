#' @title Extract dictionary from tidy_sem
#' @description Provides access to the \code{dictionary} element of a
#' \code{tidy_sem} object. This can be used to return or assign to the
#' \code{dictionary} element.
#' @param x Object of class tidy_sem.
#' @return data.frame
#' @rdname dictionary
#' @examples
#' dict <- tidy_sem(iris, split = "\\.")
#' dictionary(dict)
#' @export
dictionary <- function(x) UseMethod("dictionary")

#' @method dictionary tidy_sem
#' @export
dictionary.tidy_sem <- function(x) x[["dictionary"]]

#' @method dictionary data_dict
#' @export
dictionary.data_dict <- dictionary.tidy_sem

#' @rdname dictionary
#' @param value A valid value for \code{dictionary(x)}.
#' @export
`dictionary<-` <- function(x, value){
  UseMethod("dictionary<-")
}

#' @method dictionary<- tidy_sem
#' @export
`dictionary<-.tidy_sem` <- function(x, value)
{
  x$dictionary <- value
  x
}

#' @title Extract syntax from tidy_sem
#' @description Provides access to the \code{syntax} element of a
#' \code{tidy_sem} object. This can be used to return or assign to the
#' \code{syntax} element.
#' @param x Object of class tidy_sem.
#' @return data.frame
#' @examples
#' dict <- tidy_sem(iris, split = "\\.")
#' dict <- add_paths(dict, Sepal.Width ~~ Sepal.Length)
#' syntax(dict)
#' @rdname syntax
#' @export
syntax <- function(x) UseMethod("syntax")

#' @method syntax tidy_sem
#' @export
syntax.tidy_sem <- function(x) x[["syntax"]]

#' @rdname syntax
#' @param value A valid value for \code{syntax(x)}.
#' @export
`syntax<-` <- function(x, value){
  UseMethod("syntax<-")
}

#' @method syntax<- tidy_sem
#' @export
`syntax<-.tidy_sem` <- function(x, value)
{
  x$syntax <- value
  x
}

#' @title Extract data from tidy_sem
#' @description Provides access to the \code{data} element of a
#' \code{tidy_sem} object. This can be used to return or assign to the
#' \code{data} element.
#' @param x Object of class tidy_sem.
#' @return data.frame
#' @rdname get_data
#' @examples
#' dict <- tidy_sem(iris, split = "\\.")
#' get_data(dict)
#' @export
get_data <- function(x) UseMethod("get_data")

#' @method get_data tidy_sem
#' @export
get_data.tidy_sem <- function(x) x[["data"]]

#' @method get_data data_dict
#' @export
get_data.data_dict <- get_data.tidy_sem

#' @rdname get_data
#' @param value A valid value for \code{get_data(x)}.
#' @export
`get_data<-` <- function(x, value){
  UseMethod("get_data<-")
}

#' @method get_data<- tidy_sem
#' @export
`get_data<-.tidy_sem` <- function(x, value)
{
  x$data <- value
  x
}
