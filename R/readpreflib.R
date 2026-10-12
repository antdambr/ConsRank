#' Read a PrefLib data file
#'
#' Read a file of preferences in the format of PrefLib (www.preflib.org) and
#' return the rankings in the format used by the functions of the package: a
#' matrix with one row per distinct ranking and one column per object, and the
#' frequency of each ranking.
#'
#' @param file The name of a PrefLib data file, a URL (for example the raw link
#'   of a file of the PrefLib-Data repository on GitHub) or a connection
#' @param unranked What to do with the objects that a judge does not rank, in the
#'   files of incomplete orders (soi, toi) and of categorical preferences (cat):
#'   \code{"last"} (default) ties them in the last position, after the ranked
#'   objects, as in the files that PrefLib obtains by adding the unranked
#'   alternatives at the bottom; \code{"na"} leaves them missing (NA)
#' @param expand Logical. If FALSE (default), the rankings are returned once each
#'   with their frequencies in wk; if TRUE, each ranking is repeated as many times
#'   as its frequency, one row per judge, and wk is NULL
#'
#' @return A list containing the following components:
#' \tabular{lll}{
#'   X \tab \tab the rankings: a matrix with one row per distinct ranking (one per
#'   judge with expand = TRUE) and one column per object, named with the names of
#'   the alternatives; the entries are the positions of the objects (1 = first
#'   position, objects tied in the same bucket have the same position, the
#'   positions of each ranking are 1, 2, 3, ... without gaps)\cr
#'   wk \tab \tab the frequency of each ranking (NULL with expand = TRUE)\cr
#'   datatype \tab \tab the data type of the file: "soc", "soi", "toc", "toi" or
#'   "cat"\cr
#'   categories \tab \tab the names of the categories, in order (cat files only,
#'   NULL otherwise)\cr
#'   header \tab \tab the metadata of the file: a named list with the values of
#'   the header lines
#' }
#'
#' @details
#' PrefLib files begin with header lines starting with "#" (file name, title, data
#' type, number of alternatives, names of the alternatives, number of voters, ...),
#' followed by one line for each distinct preference, of the form
#' \code{count: preference}. The data types are: strict orders complete (soc),
#' strict orders incomplete (soi), orders with ties complete (toc), orders with
#' ties incomplete (toi) and categorical preferences (cat). In the orders, the
#' alternatives are separated by commas from the first to the last, and tied
#' alternatives are grouped in braces: \code{13: 1, {4, 3}, 2} means that 13
#' voters put alternative 1 first, 4 and 3 tied in second position and 2 last.
#' In the categorical preferences each line lists the categories in order, from
#' the most preferred one, and an empty category is \code{{}}: the categories
#' become the buckets of a ranking with ties, and the empty categories are
#' skipped. Weighted matching files (wmd) do not contain rankings and are not
#' read.
#'
#' The function checks that the sum of the frequencies is equal to the number of
#' voters declared in the header. With \code{unranked = "last"}, two different
#' incomplete orders can give the same ranking (for example the order of all the
#' objects but one, and the full order): such rankings are merged and their
#' frequencies are summed. The rankings are in the order in which they first
#' appear in the file.
#'
#' With \code{unranked = "na"} the result can be used with \code{\link{tau_x}}
#' and \code{\link{combinpmatr}}, which treat a pair with a missing object as in
#' Emond and Mason (2002); \code{\link{kemenyd}} and the median ranking functions
#' need complete rankings, as given by \code{unranked = "last"}.
#'
#' @examples
#' # a small file in the PrefLib format: 4 objects, orders with ties
#' tf <- tempfile(fileext = ".toi")
#' writeLines(c(
#'   "# FILE NAME: example.toi",
#'   "# TITLE: A small example",
#'   "# DATA TYPE: toi",
#'   "# NUMBER ALTERNATIVES: 4",
#'   "# ALTERNATIVE NAME 1: a",
#'   "# ALTERNATIVE NAME 2: b",
#'   "# ALTERNATIVE NAME 3: c",
#'   "# ALTERNATIVE NAME 4: d",
#'   "# NUMBER VOTERS: 10",
#'   "# NUMBER UNIQUE ORDERS: 3",
#'   "5: 1, {4, 3}, 2",
#'   "3: 2, 1",
#'   "2: 3, 1, 2, 4"), tf)
#' P <- readpreflib(tf)
#' P$X
#' P$wk
#' CR <- consrank(P$X, wk = P$wk, algorithm = "BB")
#' CR$Consensus
#' # the objects not ranked by a judge left missing
#' readpreflib(tf, unranked = "na")$X
#' \dontrun{
#' # a file of the PrefLib-Data repository: the sushi data (top-10 of 100 sushi)
#' P <- readpreflib(paste0("https://raw.githubusercontent.com/PrefLib/",
#'   "PrefLib-Data/main/datasets/00014%20-%20sushi/00014-00000002.soi"))
#' dim(P$X)
#' }
#'
#' @references
#' Mattei, N., and Walsh, T. (2013). PrefLib: A library for preferences
#' http://www.preflib.org. In Algorithmic Decision Theory (ADT 2013), Lecture
#' Notes in Computer Science, vol. 8176, 259-270. Springer.
#'
#' Emond, E. J., and Mason, D. W. (2002). A new rank correlation coefficient with
#' application to the consensus ranking problem. Journal of Multi-Criteria
#' Decision Analysis, 11(1), 17-28.
#'
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#'
#' @seealso \code{\link{consrank}}, \code{\link{mcbo}} and
#'   \code{\link{consranktopk}} for the median ranking of the data
#' @seealso \code{\link{tabulaterows}} for the frequencies of the rankings in a
#'   data matrix
#'
#' @keywords PrefLib
#' @keywords Preference data
#'
#' @export

readpreflib <- function(file, unranked = c("last", "na"), expand = FALSE) {

  # the choice for the unranked objects, "last" by default
  unranked <- match.arg(unranked)
  # expand must be TRUE or FALSE
  if (!is.logical(expand) || length(expand) != 1 || is.na(expand)) {
    stop("expand must be TRUE or FALSE")
  }

  # read all the lines of the file (a file name, a URL or a connection)
  lines <- readLines(file, warn = FALSE, encoding = "UTF-8")
  # remove the spaces at the end of the lines and the empty lines
  lines <- sub("[[:space:]]+$", "", lines)
  lines <- lines[nzchar(lines)]

  # ---- header: the lines starting with "#" ----
  ishead <- startsWith(lines, "#")
  # a PrefLib file in the current format always begins with the header
  if (!any(ishead)) {
    stop("the file has no header lines starting with '#': it is not a PrefLib ",
         "file in the current format (see https://preflib.github.io/PrefLib-Jekyll/",
         "format). Files in the old PrefLib format are not supported")
  }
  hdr <- lines[ishead]
  # the name of each metadata: the text between "#" and the first ":"
  keys <- toupper(trimws(sub("^#([^:]*):.*$", "\\1", hdr)))
  # the value of each metadata: the text after the first ":"
  vals <- trimws(sub("^#[^:]*:?", "", hdr))
  # the metadata as a named list
  header <- as.list(vals)
  names(header) <- keys

  # ---- data type ----
  datatype <- tolower(header[["DATA TYPE"]])
  # without the metadata, the extension of the file name gives the data type
  if (length(datatype) == 0 || !nzchar(datatype)) {
    datatype <- if (is.character(file)) tolower(sub("^.*\\.", "", basename(file))) else ""
  }
  # matching data are graphs, not rankings
  if (datatype == "wmd") {
    stop("the file contains weighted matching data (wmd), not rankings: readpreflib ",
         "reads only the data types soc, soi, toc, toi and cat")
  }
  # any other data type is unknown
  if (!datatype %in% c("soc", "soi", "toc", "toi", "cat")) {
    stop("unknown data type '", datatype, "': readpreflib reads only the data ",
         "types soc, soi, toc, toi and cat (line '# DATA TYPE:' of the header)")
  }

  # ---- alternatives ----
  # the number of alternatives declared in the header
  M <- suppressWarnings(as.integer(header[["NUMBER ALTERNATIVES"]]))
  if (length(M) == 0 || is.na(M) || M < 1) {
    stop("the header does not give the number of alternatives (line ",
         "'# NUMBER ALTERNATIVES:'): it is required to read the file")
  }
  # the lines with the names of the alternatives, "ALTERNATIVE NAME i"
  isname <- grepl("^ALTERNATIVE NAME [0-9]+$", keys)
  # by default the alternatives are named with their numbers
  items <- as.character(seq_len(M))
  # the number of each named alternative
  idx <- as.integer(sub("^ALTERNATIVE NAME ", "", keys[isname]))
  # the names given in the header replace the numbers
  ok <- idx >= 1 & idx <= M
  items[idx[ok]] <- vals[isname][ok]
  # PrefLib does not allow two alternatives with the same name; if it happens,
  # the names are made unique so that the columns can be told apart
  if (anyDuplicated(items)) {
    warning("some alternatives have the same name in the header: the names have ",
            "been made unique")
    items <- make.unique(items)
  }

  # the names of the categories, in order (cat files only)
  categories <- NULL
  if (datatype == "cat") {
    # the lines "CATEGORY NAME j"
    iscat <- grepl("^CATEGORY NAME [0-9]+$", keys)
    # the number of each category
    cidx <- as.integer(sub("^CATEGORY NAME ", "", keys[iscat]))
    # the names in the order of the categories
    categories <- vals[iscat][order(cidx)]
  }

  # ---- preferences: the lines "count: preference" ----
  dat <- lines[!ishead]
  if (length(dat) == 0) stop("the file contains no preferences after the header")
  # every line must contain ":" between the count and the preference
  bad <- which(!grepl(":", dat, fixed = TRUE))
  if (length(bad) > 0) {
    stop("line ", which(!ishead)[bad[1]], " of the file is not of the form ",
         "'count: preference'")
  }
  # the position of the first ":" in each line
  colon <- regexpr(":", dat, fixed = TRUE)
  # the count: the text before ":"
  counts <- suppressWarnings(as.numeric(substr(dat, 1, colon - 1)))
  if (anyNA(counts) || any(counts <= 0)) {
    stop("line ", which(!ishead)[which(is.na(counts) | counts <= 0)[1]],
         " of the file does not begin with a positive count")
  }
  # the preference: the text after ":"
  pref <- substring(dat, colon + 1)
  # each preference split into its parts: a bucket in braces or a single
  # alternative, in order from the first to the last
  tok <- regmatches(pref, gregexpr("\\{[^}]*\\}|[^,{}[:space:]]+", pref))

  # number of lines of preferences
  N <- length(dat)
  # for each line, the alternatives and the position of their bucket
  parsed <- lapply(tok, function(tk) {
    # the alternatives of each part, without braces and spaces
    ids <- strsplit(gsub("[{}[:space:]]", "", tk), ",", fixed = TRUE)
    # the parts that contain at least one alternative (empty categories excluded)
    len <- lengths(ids)
    # the position of each non empty part: 1, 2, 3, ... without gaps
    pos <- cumsum(len > 0)
    # each alternative with the position of its part
    list(id = unlist(ids), pos = rep(pos, len))
  })
  # the line of each alternative, the alternative and its position
  nper <- vapply(parsed, function(p) length(p$id), integer(1))
  row <- rep(seq_len(N), nper)
  id <- suppressWarnings(as.integer(unlist(lapply(parsed, `[[`, "id"))))
  pos <- unlist(lapply(parsed, `[[`, "pos"))
  # the alternatives must be numbers between 1 and the number of alternatives
  bad <- is.na(id) | id < 1 | id > M
  if (any(bad)) {
    stop("line ", which(!ishead)[row[which(bad)[1]]], " of the file contains an ",
         "alternative that is not a number between 1 and ", M)
  }
  # an alternative cannot appear twice in the same preference
  dup <- duplicated(cbind(row, id))
  if (any(dup)) {
    stop("line ", which(!ishead)[row[which(dup)[1]]], " of the file contains the ",
         "same alternative twice")
  }

  # ---- the matrix of the rankings ----
  # one row per line of the file, one column per alternative, NA = not ranked
  X <- matrix(NA_real_, nrow = N, ncol = M, dimnames = list(NULL, items))
  # the position of each ranked alternative
  X[cbind(row, id)] <- pos

  # the files of complete orders must rank all the alternatives in every line
  if (datatype %in% c("soc", "toc") && anyNA(X)) {
    warning("the file declares complete orders (", datatype, ") but some ",
            "preferences do not rank all the alternatives: the missing ones are ",
            "treated as unranked (argument unranked)")
  }
  # unranked objects tied in the last position, after the ranked ones
  if (unranked == "last" && anyNA(X)) {
    # the last position of each ranking (0 if nothing is ranked)
    last <- suppressWarnings(apply(X, 1, max, na.rm = TRUE))
    last[!is.finite(last)] <- 0
    # the missing entries take the next position
    miss <- which(is.na(X), arr.ind = TRUE)
    X[miss] <- last[miss[, 1]] + 1
  }

  # ---- checks with the header ----
  # the number of voters declared in the header
  nv <- suppressWarnings(as.numeric(header[["NUMBER VOTERS"]]))
  if (length(nv) == 1 && !is.na(nv) && !isTRUE(all.equal(sum(counts), nv))) {
    warning("the sum of the counts (", sum(counts), ") is different from the ",
            "number of voters declared in the header (", nv, ")")
  }
  # the number of distinct preferences declared in the header
  nu <- suppressWarnings(as.numeric(c(header[["NUMBER UNIQUE ORDERS"]],
                                      header[["NUMBER UNIQUE PREFERENCES"]])[1]))
  if (length(nu) == 1 && !is.na(nu) && nu != N) {
    warning("the number of preferences in the file (", N, ") is different from ",
            "the number of unique preferences declared in the header (", nu, ")")
  }

  # ---- output ----
  if (expand) {
    # one row per judge: the frequencies must be whole numbers
    if (any(counts != round(counts))) {
      stop("the counts are not whole numbers: use expand = FALSE and the weights wk")
    }
    # each ranking repeated as many times as its frequency
    X <- X[rep(seq_len(N), counts), , drop = FALSE]
    wk <- NULL
  } else {
    # a key for each ranking, to merge the rankings that are equal
    key <- apply(X, 1, paste, collapse = ",")
    # the group of each ranking, numbered in order of first appearance
    grp <- match(key, unique(key))
    # the frequencies of the distinct rankings
    wk <- as.vector(rowsum(counts, grp, reorder = TRUE))
    # one row for each distinct ranking
    X <- X[!duplicated(grp), , drop = FALSE]
  }
  # no row names
  rownames(X) <- NULL

  list(X = X, wk = wk, datatype = datatype, categories = categories,
       header = header)
}
