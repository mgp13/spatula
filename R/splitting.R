#' @export
split_grid <- function(grid_point) {
  bbox <- grid_point$bbox
  xmid <- mean(bbox[1:2])
  ymid <- mean(bbox[3:4])
  ## NOTE: the 1e-10 term avoids duplicates
  bboxes_new <- list(
    c(bbox[1], xmid - 1e-10, bbox[3], ymid - 1e-10),
    c(xmid, bbox[2], bbox[3], ymid - 1e-10),
    c(xmid, bbox[2], ymid, bbox[4]),
    c(bbox[1], xmid - 1e-10, ymid, bbox[4])
  )
  
  ## Make four new grid points
  res <- map(bboxes_new, function(bbox_test) {
    res <- list(
      tx = dplyr::filter(
        grid_point$tx,
        between(x, bbox_test[1], bbox_test[2]) &
          between(y, bbox_test[3], bbox_test[4])
      ),
      # tx = grid_point$tx[between(x, bbox_test[1], bbox_test[2]) & between(y, bbox_test[3], bbox_test[4])],
      bbox = bbox_test,
      bbox_geom = st_rectangle(bbox_test[1], bbox_test[2], bbox_test[3], bbox_test[4])
    )
    res$n <- nrow(res$tx)
    return(res)
  })
}

#' @export
fixedWidthGrid <- function(tx,
                           windowWidth = 100) {
  # splits a tx dataframe into a grid of fixed width
  windowWidth <- as.integer(windowWidth)
  # generate x coordinates for slide windows
  slideWindows <- seq(min(tx$x),
                      max(tx$x) + windowWidth,
                      windowWidth)
  slideWindows[1] <- min(tx$x)
  # initialize data structure to hold fixed-width grid - each element corresponds to one region/window and is named as "min_x_coordinate"_"max_x_coordinate". each element of this list has 4 sub-components as described below.
  grid <- vector(mode = list, length = length(slideWindows))
  for (i in 1:(length(slideWindows) - 1)) {
    # initialize structure with 4 components:
    # region_idx = region/window identifier
    # n = number of transcripts in this region/window
    # tx = data frame of transcripts in the region/window
    # bbox = bbox of the window
    grid[i] <- vector(mode = list, length = 3)
    names(grid[i]) <- c('region_idx', 'n', 'tx', 'bbox', 'density')
    grid[i]['region_idx'] <-
      paste0(slideWindows[i], "_", slideWindows[i + 1])
    grid[i]['tx'] <-
      tx[dplyr::between(x, slideWindows[i], slideWindows[i + 1])]
    grid[i]['n'] <- nrow(grid[i]['tx'])
    grid[i]['bbox'] <- st_rectangle(min(grid[i]['tx']$x),
                                    max(grid[i]['tx']$x),
                                    min(grid[i]['tx']$y),
                                    max(grid[i]['tx']$y))
    grid[i]['density'] <- grid[i]['n'] / st_area(grid[i]['bbox'])
  }
  rm(tx)
  gc()
  return(grid)
}
#' @export
split_tx <-
  function(tx,
           max_tx,
           max_voxels,
           fixedWidth = FALSE,
           windowWidth = 100) {
    if (!fixedWidth) {
      ## Initialize grid with all transcripts
      grid <- list(list(tx = tx,
                        bbox = c(
                          min(tx$x), max(tx$x), min(tx$y), max(tx$y)
                        )))
      grid[[1]]$n <- nrow(grid[[1]]$tx)
      grid[[1]]$bbox_geom <-
        st_rectangle(grid[[1]]$bbox[1], grid[[1]]$bbox[2], grid[[1]]$bbox[3], grid[[1]]$bbox[4])
      
      ## Keep splitting grid points until each has at most max_tx transcripts
      .i <- 0
      while (TRUE) {
        .i <- .i + 1
        grids_split <- which(map_int(grid, 'n') > max_tx)
        if (length(grid) >= max_voxels) {
          break
        } else if (length(grids_split) > 0) {
          i <- grids_split[1]
          grid <- append(grid, split_grid(grid[[i]]))
          grid[[i]] <- NULL
          grids_split <- which(map_int(grid, 'n') > max_tx)
        } else {
          break
        }
      }
      ## For QC purposes, compute the transcript density of each region
      for (i in seq_len(length(grid))) {
        grid[[i]]$density <- grid[[i]]$n / st_area(grid[[i]]$'bbox_geom')
      }
      return(grid)
    } else {
      return(fixedWidthGrid(tx, windowWidth = 100))
    }
  }
