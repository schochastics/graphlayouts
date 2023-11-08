ensure_igraph <- function(g) {
    if (!igraph::is_igraph(g)) {
        stop("g must be an igraph object", call. = FALSE)
    }
}

ensure_connected <- function(g) {
    if (!igraph::is_connected(g, mode = "weak")) {
        stop("only connected graphs are supported.", call. = FALSE)
    }
}

get_bbox <- function(xy) {
    lbottom <- c(min(xy[, 1]), min(xy[, 2]))
    rtop <- c(max(xy[, 1]), max(xy[, 2]))
    c(lbottom, rtop)
}

mv_to_null <- function(xy) {
    bbox <- get_bbox(xy)
    xy[, 1] <- xy[, 1] - bbox[1]
    xy[, 2] <- xy[, 2] - bbox[2]
    xy
}

scale_to_100 <- function(x) {
    a <- min(x)
    b <- max(x)
    100 / (b - a) * x - 100 / (b - a) * a
}

interpolate_cent <- function(cent, x) {
    a <- min(cent)
    b <- max(cent)
    alpha <- 100 / (b - a)
    beta <- -100 / (b - a) * a
    (x - beta) / alpha
}

map_to_angle_range <- function(xy, arange) {
    angles <- atan2(xy[, 2], xy[, 1]) / pi * 180
    angles[angles < 0] <- abs(angles[angles < 0]) + 180
    radii <- sqrt(rowSums(xy^2))
    angles <- normalise(angles, to = arange)
    angles <- angles * pi / 180
    cbind(radii * cos(angles), radii * sin(angles))
}

normalise <- function(x, from = range(x), to = c(0, 1)) {
    x <- (x - from[1]) / (from[2] - from[1])
    if (!identical(to, c(0, 1))) {
        x <- x * (to[2] - to[1]) + to[1]
    }
    x
}

get_seed <- function() {
    if (exists(".Random.seed", .GlobalEnv)) {
        return(.GlobalEnv$.Random.seed)
    } else {
        return(NULL)
    }
}

restore_seed <- function(oldseed) {
    if (!is.null(oldseed)) {
        .GlobalEnv$.Random.seed <- oldseed
    } else {
        rm(".Random.seed", envir = .GlobalEnv)
    }
}

iso_project <- function(xyz, a = 35.264, b = 45) {
    a <- a * pi / 180
    b <- b * pi / 180
    T1 <- matrix(c(1, 0, 0, 0, cos(a), sin(a), 0, -sin(a), cos(a)), 3, 3, byrow = TRUE)
    T2 <- matrix(c(cos(b), 0, -sin(b), 0, 1, 0, sin(b), 0, cos(b)), 3, 3, byrow = TRUE)
    trans <- T1 %*% T2 %*% t(xyz)
    coords2D <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 0), 3, 3, byrow = TRUE) %*% trans
    coords2D <- t(coords2D)
    return(coords2D[, 1:2])
}

degree_lvl <- function(g) {
    if (!"lvl" %in% igraph::vertex_attr_names(g)) {
        stop("level information should be stored in a vertex attribute called 'lvl'")
    }
    A <- igraph::as_adj(g, "both", sparse = FALSE)
    lvl_mat <- outer(igraph::V(g)$lvl, igraph::V(g)$lvl, function(x, y) x == y)
    cent <- rbind(rowSums(A * lvl_mat), rowSums(A * !lvl_mat))
    rownames(cent) <- c("intra", "inter")
    cent
}

optim_isolates <- function(g, xyz) {
    deg <- degree_lvl(g)
    idx <- which(deg[1, ] == 0)
    if (length(idx) > 0) {
        neigh <- igraph::neighborhood(g, order = 1, nodes = idx, mode = "all", mindist = 1)
        xyz[idx, c(1, 3)] <- do.call("rbind", lapply(neigh, function(id) cbind(mean(xyz[id, 1]), mean(xyz[id, 3]))))
    }
    xyz
}

optim_rotation <- function(g, xyz) {
    D <- igraph::distances(g)
    W <- 1 / D^2
    smin <- stress3D(xyz, W, D)
    # smin <- stress(xyz[,c(1,3)],W,D)
    amax <- 0
    idx <- which(igraph::V(g)$lvl == 2)
    for (alpha in seq(0, 360, 5)) {
        xyz_new <- xyz
        xyz_new[idx, c(1, 3)] <- layout_rotate(xyz_new[idx, c(1, 3)], alpha)
        stemp <- stress3D(xyz_new, W, D)
        # stemp <- stress(xyz_new[,c(1,3)],W,D)
        if (stemp < smin) {
            amax <- alpha
        }
    }
    xyz[idx, c(1, 3)] <- layout_rotate(xyz[idx, c(1, 3)], amax)
    xyz
}

optim_level <- function(g, lvl, xy) {
    A <- igraph::as_adj(g, "both")
    Ainter <- A[igraph::V(g)$lvl != lvl, igraph::V(g)$lvl == lvl]
    adjList <- apply(Ainter, 1, function(x) which(x == 1))

    xy2 <- do.call("rbind", lapply(adjList, function(id) cbind(mean(xy[id, 1]), mean(xy[id, 2]))))
    idx <- is.na(xy2[, 1])

    mx <- mean(xy[, 1], na.rm = TRUE)
    my <- mean(xy[, 2], na.rm = TRUE)
    r <- max(sqrt((xy[, 1] - mx)^2 + (xy[, 2] - my)^2))

    if (length(idx) > 0) {
        xy2[idx, 1] <- stats::runif(n = sum(idx), min = mx - r, max = mx + r)
        xy2[idx, 2] <- sample(c(-1, 1), sum(idx), replace = T) * sqrt(r^2 - (xy2[idx, 1] - mx)^2) + my
    }
    xy2
}

find_lastD <- function(DList, i, j, k) {
    for (l in seq((k - 1), 1)) {
        if (l == 0) {
            next()
        }
        if (!is.infinite(DList[[l]][i, j])) {
            return(list(value = DList[[l]][i, j], index = l))
        }
    }
    return(list(value = Inf, index = NULL))
}

find_nextD <- function(DList, i, j, k) {
    for (l in seq((k + 1), length(DList))) {
        if (l > length(DList)) {
            break()
        }
        if (!is.infinite(DList[[l]][i, j])) {
            return(list(value = DList[[l]][i, j], index = l))
        }
    }
    return(list(value = Inf, index = NULL))
}

adjust_value <- function(lastD, nextD, k, tlast, tnext, n) {
    if (!is.infinite(lastD$value) && !is.infinite(nextD$value)) {
        beta <- (k - tlast) / (tnext - tlast)
        return((1 - beta) * lastD$value + beta * nextD$value + 1)
    } else if (is.infinite(lastD$value) && !is.infinite(nextD$value)) {
        return(nextD$value + 1)
    } else if (!is.infinite(lastD$value) && is.infinite(nextD$value)) {
        return(lastD$value + 1)
    } else {
        return(sqrt(n))
    }
}

adjust_dist <- function(DList) {
    n <- nrow(DList[[1]])
    for (k in seq_along(DList)) {
        for (i in 1:n) {
            for (j in 1:n) {
                if (is.infinite(DList[[k]][i, j])) {
                    lastD <- find_lastD(DList, i, j, k)
                    nextD <- find_nextD(DList, i, j, k)
                    DList[[k]][i, j] <- adjust_value(lastD, nextD, k, lastD$index, nextD$index, n)
                }
            }
        }
    }
    DList
}
