#' Plot rankings on a permutation polytope of 3 or 4 objects
#'
#' @param X Sample of rankings (matrix), either one row for each judge or the distinct rankings
#'   returned by tabulaterows() with their frequencies in Wk
#' @param L Labels of the objects (character vector)
#' @param Wk Frequency weights for each ranking (numeric vector). If NULL, each row of X counts 1
#' @param nobj Number of objects: 3 or 4 (integer)
#' @param engine Rendering engine for 4-object plots: "auto", "rgl", or "plotly" (character)
#' @param size Scale factor of the symbols (default 1): values larger than 1 enlarge them
#'
#' @return For nobj=3: invisible NULL (creates base graphics plot)
#'         For nobj=4: invisible NULL for rgl, plotly object for plotly engine
#'
#' @details
#' For 3 objects: uses base R graphics (always available)
#' For 4 objects:
#'   - engine="auto": tries rgl first, falls back to plotly if unavailable
#'   - engine="rgl": requires 'rgl' package (may need XQuartz on macOS)
#'   - engine="plotly": requires 'plotly' package (works on all platforms)
#'
#' Each observed ranking is drawn as a circle (3 objects) or a sphere (4 objects) whose area is
#' proportional to its frequency, i.e. the sum of the weights Wk of the rows of X equal to it (rows
#' with ties are matched whatever the coding of the ties, e.g. (1, 3, 3) and (1, 2, 2)); the most
#' frequent ranking has the largest symbol, whose size is controlled by \code{size}. If X is NULL,
#' all the points of the polytope are shown. With the plotly engine the frequency is shown when the
#' mouse is over a point.
#'
#' @examples
#' data(BU)
#' polyplotnew(BU[, 1:3], Wk = BU[, 4])
#' \donttest{
#' data(German)
#' tab <- tabulaterows(German)
#' if (requireNamespace("rgl", quietly = TRUE)) {
#'   polyplotnew(tab$X, Wk = tab$Wk, nobj = 4, engine = "rgl")
#' }
#' }
#'
#' @references 
#' Thompson, G. L. (1993). Generalized permutation polytopes and exploratory 
#' graphical methods for ranked data. The Annals of Statistics, 1401-1430.
#'
#' Heiser, W. J., and D'Ambrosio, A. (2013). Clustering and prediction of 
#' rankings within a Kemeny distance framework. In Algorithms from and for 
#' Nature and Life (pp. 19-31). Springer International Publishing.
#'
#' @author Antonio D'Ambrosio, Sonia Amodio
#'
#' @export

polyplotnew <- function(X = NULL, L = NULL, Wk = NULL, nobj = 3, engine = "auto", size = 1) {

  # Input validation
  if (!nobj %in% c(3, 4)) {
    stop("nobj must be either 3 or 4")
  }

  if (!engine %in% c("auto", "rgl", "plotly")) {
    stop("engine must be 'auto', 'rgl', or 'plotly'")
  }

  if (!is.numeric(size) || length(size) != 1 || size <= 0) {
    stop("size must be a positive number")
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # 3 OBJECTS: Use base graphics
  # ═══════════════════════════════════════════════════════════════════════════
  
  if (nobj == 3) {
    return(plot_polytope_3d(X, L, Wk, size))
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # 4 OBJECTS: Choose rendering engine
  # ═══════════════════════════════════════════════════════════════════════════
  
  if (engine == "auto") {
    engine <- choose_engine()
  }
  
  if (engine == "rgl") {
    return(plot_polytope_4d_rgl(X, L, Wk, size))
  } else {
    return(plot_polytope_4d_plotly(X, L, Wk, size))
  }
}


# ═══════════════════════════════════════════════════════════════════════════
# HELPER: Choose best available engine
# ═══════════════════════════════════════════════════════════════════════════

choose_engine <- function() {
  # Try rgl first
  if (requireNamespace("rgl", quietly = TRUE)) {
    # Test if rgl actually works (may fail on Mac without XQuartz)
    tryCatch({
      rgl::open3d()
      rgl::close3d()
      message("Using rgl for 3D rendering")
      return("rgl")
    }, error = function(e) {
      message("rgl not available (may need XQuartz on macOS)")
    })
  }
  
  # Fallback to plotly
  if (requireNamespace("plotly", quietly = TRUE)) {
    message("Using plotly for 3D rendering")
    return("plotly")
  }
  
  # No engine available
  stop(
    "No 3D rendering engine available.\n",
    "Please install either:\n",
    "  - rgl: install.packages('rgl')\n",
    "  - plotly: install.packages('plotly')\n",
    "Note: rgl on macOS requires XQuartz: brew install --cask xquartz"
  )
}


# ═══════════════════════════════════════════════════════════════════════════
# 3D PLOT (3 objects) - Base graphics
# ═══════════════════════════════════════════════════════════════════════════

plot_polytope_3d <- function(X, L, Wk, size = 1) {
  
  # Define polytope structure
  ranks <- rbind(
    c(1, 2, 3),  # ABC
    c(1, 2, 2),  # A(BC)
    c(1, 3, 2),  # ACB
    c(1, 2, 1),  # (AC)B
    c(2, 3, 1),  # CAB
    c(2, 2, 1),  # C(AB)
    c(3, 2, 1),  # CBA
    c(2, 1, 1),  # (BC)A
    c(3, 1, 2),  # BCA
    c(2, 1, 2),  # B(AC)
    c(2, 1, 3),  # BAC
    c(1, 1, 2),  # (AB)C
    c(1, 1, 1)   # (ABC)
  )
  
  # Polytope coordinates
  coord <- rbind(
    c(0.000000, 0.408248),     # 1: ABC
    c(0.176777, 0.306186),     # 2: A(BC)
    c(0.353553, 0.204124),     # 3: ACB
    c(0.353553, 0.000000),     # 4: (AC)B
    c(0.353553, -0.204124),    # 5: CAB
    c(0.176777, -0.306186),    # 6: C(AB)
    c(0.000000, -0.408248),    # 7: CBA
    c(-0.176777, -0.306186),   # 8: (BC)A
    c(-0.353553, -0.204124),   # 9: BCA
    c(-0.353553, 0.000000),    # 10: B(AC)
    c(-0.353553, 0.204124),    # 11: BAC
    c(-0.176777, 0.306186),    # 12: (AB)C
    c(0.000000, 0.000000)      # 13: (ABC)
  )
  
  # Generate labels
  rr <- labelsn(ranks, 3, if (is.null(L)) 1:3 else L, labs = if (is.null(L)) 2 else 1)
  
  # Observed rankings and their frequencies
  pw <- polytope_weights(X, Wk, ranks)

  # Create base plot
  plot(coord, ylim = c(-0.5, 0.5), xlim = c(-0.5, 0.5),
       axes = FALSE, ann = FALSE, asp = 1)

  # Draw polytope edges
  draw_polytope_edges_3d(coord)

  # Adjust text coordinates for readability
  tcoord <- adjust_text_coordinates_3d(coord)

  # Plot points
  plot_points_3d(coord, pw, rr, tcoord, size)
  
  invisible(NULL)
}


# ═══════════════════════════════════════════════════════════════════════════
# HELPER: Draw 3D polytope edges
# ═══════════════════════════════════════════════════════════════════════════

draw_polytope_edges_3d <- function(coord) {
  # Main perimeter
  lines(coord[1:11, ])
  lines(c(coord[11, 1], coord[1, 1]), c(coord[11, 2], coord[1, 2]))
  
  # Internal dashed lines
  lines(c(coord[2, 1], coord[8, 1]), c(coord[2, 2], coord[8, 2]), lty = 2)
  lines(c(coord[4, 1], coord[10, 1]), c(coord[4, 2], coord[10, 2]), lty = 2)
  lines(c(coord[6, 1], coord[12, 1]), c(coord[6, 2], coord[12, 2]), lty = 2)
}


# ═══════════════════════════════════════════════════════════════════════════
# HELPER: Adjust text coordinates for 3D plot
# ═══════════════════════════════════════════════════════════════════════════

adjust_text_coordinates_3d <- function(coord) {
  tcoord <- coord
  
  # Top and center
  top_center <- c(1, 2, 3, 11, 12, 13)
  tcoord[top_center, 2] <- tcoord[top_center, 2] + 0.1
  
  # Bottom
  bottom <- c(5, 6, 7, 8, 9)
  tcoord[bottom, 2] <- tcoord[bottom, 2] - 0.1
  
  # Sides
  tcoord[4, 1] <- tcoord[4, 1] + 0.1   # Right
  tcoord[10, 1] <- tcoord[10, 1] - 0.1 # Left
  
  return(tcoord)
}


# ═══════════════════════════════════════════════════════════════════════════
# HELPER: Plot points for 3D
# ═══════════════════════════════════════════════════════════════════════════

plot_points_3d <- function(coord, pw, rr, tcoord, size = 1) {

  indplot <- pw$ind

  if (is.null(pw$w)) {
    # all the points of the polytope
    points(coord[indplot, 1], coord[indplot, 2],
           pch = 16, cex = 0.8 * size, col = "blue")
  } else {
    # circles with area proportional to the frequencies
    sizes <- pmax(0.5, 6 * size * relsize(pw$w))
    points(coord[indplot, 1], coord[indplot, 2],
           pch = 21, cex = sizes, col = "blue",
           bg = grDevices::adjustcolor("blue", alpha.f = 0.6))
  }

  # Add labels
  text(tcoord[indplot, , drop = FALSE], rr[indplot, ])
}


# ═══════════════════════════════════════════════════════════════════════════
# 4D PLOT - rgl version
# ═══════════════════════════════════════════════════════════════════════════

plot_polytope_4d_rgl <- function(X, L, Wk, size = 1) {

  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("rgl package is required but not available")
  }

  # Get polytope structure
  structure <- create_polytope_4d_structure()

  # Generate labels
  rr <- labelsn(structure$ranks, 4,
                if (is.null(L)) 1:4 else L,
                labs = if (is.null(L)) 2 else 1)

  # Observed rankings and their frequencies
  pw <- polytope_weights(X, Wk, structure$ranks)

  # Initialize 3D plot
  rgl::plot3d(structure$coord, type = 'p', xlab = '', ylab = '', zlab = '',
              aspect = TRUE, box = FALSE, axes = FALSE, col = 1)

  # Draw all edges
  draw_polytope_edges_4d_rgl(structure)

  # Plot points/spheres
  radius <- plot_spheres_rgl(structure$coord, pw, size)

  # Add text labels, moved outwards (from the center of the polytope) beyond the spheres
  xyz <- structure$coord[pw$ind, , drop = FALSE]
  dir <- sweep(xyz, 2, structure$coord[nrow(structure$coord), ])
  len <- sqrt(rowSums(dir^2))
  len[len == 0] <- 1
  rgl::text3d(xyz + dir / len * (radius + 0.15), texts = rr[pw$ind, ],
              col = 1, cex = 0.7)

  invisible(NULL)
}


# ═══════════════════════════════════════════════════════════════════════════
# 4D PLOT - plotly version
# ═══════════════════════════════════════════════════════════════════════════

plot_polytope_4d_plotly <- function(X, L, Wk, size = 1) {
  
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("plotly package is required but not available")
  }
  
  # Get polytope structure
  structure <- create_polytope_4d_structure()
  
  # Generate labels
  rr <- labelsn(structure$ranks, 4, 
                if (is.null(L)) 1:4 else L, 
                labs = if (is.null(L)) 2 else 1)
  
  # Observed rankings and their frequencies
  pw <- polytope_weights(X, Wk, structure$ranks)

  # Create plotly object
  p <- plotly::plot_ly()

  # Add edges
  p <- add_polytope_edges_plotly(p, structure)

  # Add points
  p <- add_points_plotly(p, structure$coord, pw, rr, size)
  
  # Configure layout
  p <- p |>plotly::layout(
    scene = list(
      xaxis = list(showticklabels = FALSE, showgrid = FALSE, 
                   zeroline = FALSE, title = ""),
      yaxis = list(showticklabels = FALSE, showgrid = FALSE, 
                   zeroline = FALSE, title = ""),
      zaxis = list(showticklabels = FALSE, showgrid = FALSE, 
                   zeroline = FALSE, title = ""),
      aspectmode = "cube"
    ),
    showlegend = FALSE
  )
  
  return(p)
}


# ═══════════════════════════════════════════════════════════════════════════
# HELPER: Create 4D polytope structure
# ═══════════════════════════════════════════════════════════════════════════

create_polytope_4d_structure <- function() {
  
  # Hexagon A first
  E1 <- rbind(
    c(0.5, 0.5, 1.4142135),
    c(1.0, 1.0, 0.70710677),
    c(1.5, 0.5, 0.0),
    c(1.5, -0.5, 0.0),
    c(1.0, -1.0, 0.70710677),
    c(0.5, -0.5, 1.4142135)
  )
  MA <- colMeans(E1)
  E1T <- rbind(
    c(1.25, 0.75, 0.3536),
    c(0.75, -0.75, 1.0607),
    c(0.5, 0, 1.4142),
    c(1.5, 0, 0),
    c(1.25, -0.75, 0.3536),
    c(0.75, 0.75, 1.0607)
  )
  ranksA <- rbind(
    c(1,2,3,4), c(1,3,2,4), c(1,4,2,3), c(1,4,3,2), c(1,3,4,2), c(1,2,4,3),
    c(1,2,2,2), c(1,3,2,3), c(1,2,3,2), c(1,2,3,3), c(1,3,2,2), c(1,3,3,2), c(1,2,2,3)
  )
  
  # Hexagon B first
  E2 <- rbind(
    c(-1.5, -0.5, 0.0),
    c(-1.5, 0.5, 0.0),
    c(-1.0, 1.0, 0.70710677),
    c(-0.5, 0.5, 1.4142135),
    c(-0.5, -0.5, 1.4142135),
    c(-1.0, -1.0, 0.70710677)
  )
  MB <- colMeans(E2)
  E2T <- rbind(
    c(-0.5, 0, 1.4142),
    c(-1.5, 0, 0),
    c(-1.25, 0.75, 0.3536),
    c(-0.75, -0.75, 1.0607),
    c(-1.25, -0.75, 0.3536),
    c(-0.75, 0.75, 1.0607)
  )
  ranksB <- rbind(
    c(4,1,3,2), c(4,1,2,3), c(3,1,2,4), c(2,1,3,4), c(2,1,4,3), c(3,1,4,2),
    c(2,1,2,2), c(2,1,3,3), c(3,1,2,2), c(3,1,2,3), c(2,1,3,2), c(3,1,3,2), c(2,1,2,3)
  )
  
  # Hexagon C first
  E3 <- rbind(
    c(-1.0, 1.0, -0.70710677),
    c(-0.5, 0.5, -1.4142135),
    c(0.5, 0.5, -1.4142135),
    c(1.0, 1.0, -0.70710677),
    c(0.5, 1.5, 0.0),
    c(-0.5, 1.5, 0.0)
  )
  MC <- colMeans(E3)
  E3T <- rbind(
    c(0, 0.5, -1.4142),
    c(0, 1.5, 0),
    c(0.75, 1.25, -0.3536),
    c(-0.75, 0.75, -1.0607),
    c(0.75, 0.75, -1.0607),
    c(-0.75, 1.25, -0.3536)
  )
  ranksC <- rbind(
    c(4,2,1,3), c(4,3,1,2), c(3,4,1,2), c(2,4,1,3), c(2,3,1,4), c(3,2,1,4),
    c(2,2,1,2), c(3,3,1,2), c(2,2,1,3), c(2,3,1,3), c(3,2,1,2), c(2,3,1,2), c(3,2,1,3)
  )
  
  # Hexagon D first
  E4 <- rbind(
    c(-1.0, -1.0, -0.70710677),
    c(-0.5, -0.5, -1.4142135),
    c(0.5, -0.5, -1.4142135),
    c(1.0, -1.0, -0.70710677),
    c(0.5, -1.5, 0.0),
    c(-0.5, -1.5, 0.0)
  )
  MD <- colMeans(E4)
  E4T <- rbind(
    c(0, -0.5, -1.4142),
    c(0, -1.5, 0),
    c(-0.75, -1.25, -0.3536),
    c(0.75, -0.75, -1.0607),
    c(-0.75, -0.75, -1.0607),
    c(0.75, -1.25, -0.3536)
  )
  ranksD <- rbind(
    c(4,2,3,1), c(4,3,2,1), c(3,4,2,1), c(2,4,3,1), c(2,3,4,1), c(3,2,4,1),
    c(2,2,2,1), c(3,3,2,1), c(2,2,3,1), c(3,2,3,1), c(2,3,2,1), c(3,2,2,1), c(2,3,3,1)
  )
  
  # Squares
  ESQ <- rbind(
    c(0, 0.5, 1.4142), c(0, -0.5, 1.4142),
    c(-0.5, 0, -1.4142), c(0.5, 0, -1.4142),
    c(1.25, 0.75, -0.3536), c(0.75, 1.25, 0.3536),
    c(-1.25, 0.75, -0.3536), c(-0.75, 1.25, 0.3536),
    c(1.25, -0.75, -0.3536), c(0.75, -1.25, 0.3536),
    c(-1.25, -0.75, -0.3536), c(-0.75, -1.25, 0.3536)
  )
  ranksSQ <- rbind(
    c(1,1,2,3), c(1,1,3,2), c(3,2,1,1), c(2,3,1,1),
    c(1,3,1,2), c(1,2,1,3), c(3,1,1,2), c(2,1,1,3),
    c(1,3,2,1), c(1,2,3,1), c(3,1,2,1), c(2,1,3,1)
  )
  
  # Middle points
  M_AB_CD <- colMeans(rbind(c(0,0.5,1.4142), c(-0.5,0,1.4142), 
                            c(0,-0.5,1.4142), c(0.5,0,1.4142)))
  M_AC_BD <- colMeans(rbind(c(1.25,0.75,0.3536), c(1.25,0.75,-0.3536), 
                            c(0.75,1.25,-0.3536), c(0.75,1.25,0.3536)))
  M_BC_AD <- colMeans(rbind(c(-0.75,1.25,-0.3536), c(-1.25,0.75,-0.3536), 
                            c(-1.25,0.75,0.3536), c(-0.75,1.25,0.3536)))
  
  EMID <- rbind(
    M_AB_CD, M_AC_BD, M_BC_AD,
    M_AB_CD * -1, M_AC_BD * -1, M_BC_AD * -1,
    MA * -1, MB * -1, MC * -1, MD * -1
  )
  
  rankres <- rbind(
    c(1,1,2,2), c(1,2,1,2), c(2,1,1,2),
    c(2,2,1,1), c(2,1,2,1), c(1,2,2,1),
    c(2,1,1,1), c(1,2,1,1), c(1,1,2,1), c(1,1,1,2),
    c(1,1,1,1)
  )
  
  # Combine all
  EE <- rbind(E1, MA, E1T, E2, MB, E2T, E3, MC, E3T, E4, MD, E4T, ESQ, EMID)
  coord <- rbind(EE, colMeans(EE))
  ranks <- rbind(ranksA, ranksB, ranksC, ranksD, ranksSQ, rankres)
  
  # Store edge definitions
  edges <- list(
    hexagons = list(
      E1 = E1, E2 = E2, E3 = E3, E4 = E4,
      E1T = E1T, E2T = E2T, E3T = E3T, E4T = E4T
    ),
    squares = ESQ
  )
  
  return(list(coord = coord, ranks = ranks, edges = edges))
}


# ═══════════════════════════════════════════════════════════════════════════
# HELPER: Draw 4D polytope edges (rgl)
# ═══════════════════════════════════════════════════════════════════════════

draw_polytope_edges_4d_rgl <- function(structure) {
  
  E1 <- structure$edges$hexagons$E1
  E2 <- structure$edges$hexagons$E2
  E3 <- structure$edges$hexagons$E3
  E4 <- structure$edges$hexagons$E4
  E1T <- structure$edges$hexagons$E1T
  E2T <- structure$edges$hexagons$E2T
  E3T <- structure$edges$hexagons$E3T
  E4T <- structure$edges$hexagons$E4T
  ESQ <- structure$edges$squares
  
  # Helper to draw hexagon
  draw_hex <- function(E, ET, col_main = 1, col_tie = 'gray') {
    for (i in 1:5) rgl::segments3d(E[c(i, i+1), ], lwd = 1, col = col_main)
    rgl::segments3d(E[c(1, 6), ], lwd = 1, col = col_main)
    rgl::segments3d(ET[c(1, 2), ], lwd = 0.5, col = col_tie)
    rgl::segments3d(ET[c(3, 4), ], lwd = 0.5, col = col_tie)
    rgl::segments3d(ET[c(5, 6), ], lwd = 0.5, col = col_tie)
  }
  
  # Draw all hexagons
  draw_hex(E1, E1T)
  draw_hex(E2, E2T)
  draw_hex(E3, E3T)
  draw_hex(E4, E4T)
  
  # Draw connecting squares
  square_connections <- list(
    c(E1[1,], E2[4,]), c(E1[6,], E2[5,]),
    c(E3[2,], E4[2,]), c(E3[3,], E4[3,]),
    c(E1[5,], E4[5,]), c(E1[4,], E4[4,]),
    c(E2[1,], E4[1,]), c(E2[6,], E4[6,]),
    c(E1[2,], E3[5,]), c(E1[3,], E3[4,]),
    c(E2[2,], E3[1,]), c(E2[3,], E3[6,])
  )
  
  for (seg in square_connections) {
    rgl::segments3d(matrix(seg, 2, 3, byrow=TRUE), lwd = 1, col = 1)
  }
  
  # Draw internal gray segments
  internal_segs <- list(
    c(ESQ[1,], ESQ[2,]), c(ESQ[3,], ESQ[4,]),
    c(ESQ[5,], ESQ[6,]), c(ESQ[7,], ESQ[8,]),
    c(ESQ[9,], ESQ[10,]), c(ESQ[11,], ESQ[12,]),
    c(E1T[1,], E3T[3,]), c(E2T[3,], E3T[6,]),
    c(E1T[3,], E2T[1,]), c(E1T[5,], E4T[6,]),
    c(E2T[5,], E4T[3,]), c(E3T[1,], E4T[1,]),
    c(ESQ[1,], E3T[2,]), c(ESQ[6,], E2T[6,]),
    c(ESQ[8,], E1T[6,]), c(ESQ[2,], E4T[2,]),
    c(ESQ[7,], E4T[5,]), c(ESQ[5,], E4T[4,]),
    c(ESQ[12,], E1T[2,]), c(ESQ[10,], E2T[4,]),
    c(ESQ[4,], E1T[4,]), c(ESQ[9,], E3T[5,]),
    c(ESQ[3,], E2T[2,]), c(ESQ[11,], E3T[4,])
  )
  
  for (seg in internal_segs) {
    rgl::segments3d(matrix(seg, 2, 3, byrow=TRUE), lwd = 0.5, col = 'gray')
  }
}


# ═══════════════════════════════════════════════════════════════════════════
# HELPER: Add edges to plotly - VERSIONE COMPLETA
# ═══════════════════════════════════════════════════════════════════════════

add_polytope_edges_plotly <- function(p, structure) {
  
  # Extract edges
  E1 <- structure$edges$hexagons$E1
  E2 <- structure$edges$hexagons$E2
  E3 <- structure$edges$hexagons$E3
  E4 <- structure$edges$hexagons$E4
  E1T <- structure$edges$hexagons$E1T
  E2T <- structure$edges$hexagons$E2T
  E3T <- structure$edges$hexagons$E3T
  E4T <- structure$edges$hexagons$E4T
  ESQ <- structure$edges$squares
  
  # Helper to add segment
  add_seg <- function(p, pt1, pt2, color = "black", width = 2) {
    p |>plotly::add_trace(
      x = c(pt1[1], pt2[1]),
      y = c(pt1[2], pt2[2]),
      z = c(pt1[3], pt2[3]),
      type = "scatter3d",
      mode = "lines",
      line = list(color = color, width = width),
      showlegend = FALSE,
      hoverinfo = "none"
    )
  }
  
  # Helper to draw hexagon
  draw_hex_plotly <- function(p, E, ET) {
    for (i in 1:5) p <- add_seg(p, E[i,], E[i+1,], "black", 2)
    p <- add_seg(p, E[1,], E[6,], "black", 2)
    p <- add_seg(p, ET[1,], ET[2,], "gray", 1)
    p <- add_seg(p, ET[3,], ET[4,], "gray", 1)
    p <- add_seg(p, ET[5,], ET[6,], "gray", 1)
    return(p)
  }
  
  # ═══════════════════════════════════════════════════════════════════════
  # DRAW ALL HEXAGONS (4 hexagons with main + internal edges)
  # ═══════════════════════════════════════════════════════════════════════
  
  p <- draw_hex_plotly(p, E1, E1T)
  p <- draw_hex_plotly(p, E2, E2T)
  p <- draw_hex_plotly(p, E3, E3T)
  p <- draw_hex_plotly(p, E4, E4T)
  
  # ═══════════════════════════════════════════════════════════════════════
  # CONNECTING SQUARES (black lines connecting the hexagons)
  # ═══════════════════════════════════════════════════════════════════════
  
  square_connections <- list(
    list(E1[1,], E2[4,]), list(E1[6,], E2[5,]),  # E1-E2 connections
    list(E3[2,], E4[2,]), list(E3[3,], E4[3,]),  # E3-E4 connections
    list(E1[5,], E4[5,]), list(E1[4,], E4[4,]),  # E1-E4 connections
    list(E2[1,], E4[1,]), list(E2[6,], E4[6,]),  # E2-E4 connections
    list(E1[2,], E3[5,]), list(E1[3,], E3[4,]),  # E1-E3 connections
    list(E2[2,], E3[1,]), list(E2[3,], E3[6,])   # E2-E3 connections
  )
  
  for (pair in square_connections) {
    p <- add_seg(p, pair[[1]], pair[[2]], "black", 2)
  }
  
  # ═══════════════════════════════════════════════════════════════════════
  # OTHER HEXAGONS (gray internal structure - ESQ squares)
  # ═══════════════════════════════════════════════════════════════════════
  
  esq_pairs <- list(
    list(ESQ[1,], ESQ[2,]),   # Square 1
    list(ESQ[3,], ESQ[4,]),   # Square 2
    list(ESQ[5,], ESQ[6,]),   # Square 3
    list(ESQ[7,], ESQ[8,]),   # Square 4
    list(ESQ[9,], ESQ[10,]),  # Square 5
    list(ESQ[11,], ESQ[12,])  # Square 6
  )
  
  for (pair in esq_pairs) {
    p <- add_seg(p, pair[[1]], pair[[2]], "gray", 1)
  }
  
  # ═══════════════════════════════════════════════════════════════════════
  # CONNECTIONS BETWEEN E*T (hexagon ties) and ESQ (squares)
  # These are ALL the remaining gray internal edges
  # ═══════════════════════════════════════════════════════════════════════
  
  internal_connections <- list(
    # E1T-E3T connections
    list(E1T[1,], E3T[3,]),
    # E2T-E3T connections
    list(E2T[3,], E3T[6,]),
    # E1T-E2T connections
    list(E1T[3,], E2T[1,]),
    # E1T-E4T connections
    list(E1T[5,], E4T[6,]),
    # E2T-E4T connections
    list(E2T[5,], E4T[3,]),
    # E3T-E4T connections
    list(E3T[1,], E4T[1,]),
    
    # ESQ-E3T connections
    list(ESQ[1,], E3T[2,]),
    # ESQ-E2T connections
    list(ESQ[6,], E2T[6,]),
    # ESQ-E1T connections
    list(ESQ[8,], E1T[6,]),
    # ESQ-E4T connections
    list(ESQ[2,], E4T[2,]),
    list(ESQ[7,], E4T[5,]),
    list(ESQ[5,], E4T[4,]),
    # ESQ-E1T connections
    list(ESQ[12,], E1T[2,]),
    # ESQ-E2T connections
    list(ESQ[10,], E2T[4,]),
    # ESQ-E1T connections
    list(ESQ[4,], E1T[4,]),
    # ESQ-E3T connections
    list(ESQ[9,], E3T[5,]),
    # ESQ-E2T connections
    list(ESQ[3,], E2T[2,]),
    # ESQ-E3T connections
    list(ESQ[11,], E3T[4,])
  )
  
  for (pair in internal_connections) {
    p <- add_seg(p, pair[[1]], pair[[2]], "gray", 1)
  }
  
  return(p)
}


# ═══════════════════════════════════════════════════════════════════════════
# HELPER: Plot spheres (rgl)
# ═══════════════════════════════════════════════════════════════════════════

plot_spheres_rgl <- function(coord, pw, size = 1) {

  if (is.null(pw$w)) {
    # all the points of the polytope
    radii <- rep(0.04 * size, length(pw$ind))
  } else {
    # spheres whose cross-section is proportional to the frequencies
    # (the edges of the polytope have length 1)
    radii <- pmax(0.03, 0.3 * size * relsize(pw$w))
  }
  rgl::spheres3d(coord[pw$ind, , drop = FALSE], col = "royalblue", radius = radii, alpha = 0.6)

  invisible(radii)
}


# ═══════════════════════════════════════════════════════════════════════════
# HELPER: Add points (plotly)
# ═══════════════════════════════════════════════════════════════════════════

add_points_plotly <- function(p, coord, pw, rr, size = 1) {

  indplot <- pw$ind

  if (is.null(pw$w)) {
    # all the points of the polytope
    sizes <- 5 * size
    hover <- rr[indplot, ]
  } else {
    # markers with area proportional to the frequencies
    sizes <- pmax(4, 40 * size * relsize(pw$w))
    hover <- paste0(rr[indplot, ], ": ", pw$w)
  }

  p <- p |>plotly::add_trace(
    x = coord[indplot, 1],
    y = coord[indplot, 2],
    z = coord[indplot, 3],
    type = "scatter3d",
    mode = "markers+text",
    marker = list(size = sizes, color = "royalblue", opacity = 0.7, sizemode = "diameter"),
    text = rr[indplot, ],
    hovertext = hover,
    textposition = "top center",
    textfont = list(size = 10, color = "black"),
    showlegend = FALSE,
    hoverinfo = "text"
  )

  return(p)
}


# ═══════════════════════════════════════════════════════════════════════════
# HELPER: Observed points of the polytope and their frequencies
# ═══════════════════════════════════════════════════════════════════════════

# The rows of X are put in dense form (1, 2, ... with tied objects sharing the
# same value, so that e.g. (1, 3, 3) becomes (1, 2, 2), the coding of ranks) and
# matched to the rows of ranks; the weights of equal rankings are summed (each
# row counts 1 if Wk is NULL). Returns the indices of the observed points of the
# polytope (ind) and their frequencies (w); if X is NULL, all the points and w = NULL.

polytope_weights <- function(X, Wk, ranks) {

  if (is.null(X)) {
    return(list(ind = seq_len(nrow(ranks)), w = NULL))
  }

  X <- as.matrix(X)
  if (ncol(X) != ncol(ranks)) {
    stop("X must have ", ncol(ranks), " columns")
  }
  if (is.null(Wk)) {
    Wk <- rep(1, nrow(X))
  }
  Wk <- as.numeric(Wk)
  if (length(Wk) != nrow(X)) {
    stop("Wk must have one value for each row of X")
  }

  dense <- t(apply(X, 1, function(r) match(r, sort(unique(r)))))
  key <- match(apply(dense, 1, paste, collapse = " "),
               apply(ranks, 1, paste, collapse = " "))
  if (anyNA(key)) {
    warning(sum(is.na(key)), " row(s) of X not plotted (missing values)")
  }

  ok <- !is.na(key)
  w <- numeric(nrow(ranks))
  w[sort(unique(key[ok]))] <- rowsum(Wk[ok], key[ok])[, 1]
  ind <- which(w > 0)

  list(ind = ind, w = w[ind])
}


# relative size of the symbols: the area is proportional to the frequency and
# the most frequent ranking has relative size 1

relsize <- function(w) {
  sqrt(w / max(w))
}


# ═══════════════════════════════════════════════════════════════════════════
# labelsn() helper (unchanged from original)
# ═══════════════════════════════════════════════════════════════════════════

labelsn <- function(x, m, label = 1:m, labs) {
  
  if (!is(x, "matrix")) {
    XX <- matrix(x, ncol = length(x))
  } else {
    XX <- x
  }
  
  nj <- nrow(XX)
  X <- XX
  
  let <- if (labs == 1) label else LETTERS[label]
  
  out <- character(nj)
  
  for (i in 1:nj) {
    ord <- rank(X[i, ])
    orders <- tapply(let, ord, sort)
    
    names1 <- character(length(orders))
    for (j in 1:length(orders)) {
      if (length(orders[[j]]) > 1) {
        names1[j] <- paste0('(', paste(orders[[j]], collapse = ' '), ')')
      } else {
        names1[j] <- as.character(orders[[j]])
      }
    }
    out[i] <- paste(names1, collapse = ' ')
  }
  
  matrix(out, nrow = nj)
}
