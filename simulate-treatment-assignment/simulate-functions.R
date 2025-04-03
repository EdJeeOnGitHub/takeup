buffer.clusters <- function(cluster.points.data, .width, proj4 = kenya.proj4) {
  cluster.points.data %>%
    # Convert to sf; if already sf, this effectively passes it through
    st_as_sf() %>%
    # Remove duplicates based on geometry
    dplyr::filter(!duplicated(st_geometry(.))) %>%
    # Transform in order to calculate distances/areas in meters
    st_transform(crs = proj4) %>%
    # Buffer each feature by .width (in crs units, e.g. meters)
    st_buffer(dist = .width, join_style = "ROUND")
}


generate.cluster.pop.area.tester <- function(min.pop.prop = 2, school.buffer.radius = 500) { 
  school.pop.area <- pi * (school.buffer.radius^2)
  
  function(clusters, ...) {
    if (is.null(clusters) || empty(clusters)) return(FALSE)
    
    gArea(clusters, byid = TRUE) >= min.pop.prop * school.pop.area
  }
}

generate.cluster.num.schools.tester <- function(primary.schools, min.num.schools = 2) {
  inner.primary.schools <- buffer.clusters(primary.schools, 1000) 
  
  function(clusters, ...) {
    contain.matrix <- gContains(clusters, inner.primary.schools, byid = TRUE) 
    
    colSums(contain.matrix) %>% 
      magrittr::is_weakly_greater_than(min.num.schools)
  }
}

generate.cluster.num.schools.tester2 <- function(primary.schools, min.num.schools = 1, min.area.frac = 0.5, school.buffer.radius = 500) {
  inner.primary.schools <- buffer.clusters(primary.schools, school.buffer.radius) 
  school.pop.area <- pi * (school.buffer.radius^2)
  
  function(cluster, cluster.schools.ids, ...) {
    length(cluster.schools.ids) > 0 && (inner.primary.schools %>% 
                                                   magrittr::extract(.$cluster.id %in% cluster.schools.ids & str_to_upper(.$county) == str_to_upper(cluster$county), ) %>% {
                                                     if (!empty(.)) gIntersection(., cluster, byid = TRUE)
                                                   } %>% {
                                                     !is.null(.) && (gArea(., byid = TRUE) %>% 
                                                                       is_weakly_greater_than(min.area.frac * school.pop.area) %>% 
                                                                       sum %>% 
                                                                       is_weakly_greater_than(min.num.schools))
        })
    }
}

generate.cluster.num.schools.tester2_sf <- function(primary.schools, 
                                                    min.num.schools = 1, 
                                                    min.area.frac = 0.5, 
                                                    school.buffer.radius = 500) {
  # 1) Buffer the primary.schools by 'school.buffer.radius'.
  #    Ensure 'primary.schools' is an sf object. If it's sp, use st_as_sf() first.
  inner.primary.schools_sf <- primary.schools %>%
    st_buffer(dist = school.buffer.radius)

  # 2) Pre-calculate the area of a circle with radius = school.buffer.radius
  school.pop.area <- pi * (school.buffer.radius^2)

  # 3) Return a function that, for each cluster, checks how many schools meet the area criterion
  function(cluster, cluster.schools.ids, ...) {
    # a) If no schools in this cluster, return FALSE early.
    if (length(cluster.schools.ids) == 0) {
      return(FALSE)
    }
    

    # b) Filter the buffered schools to match cluster IDs and county
    subset_schools_sf <- inner.primary.schools_sf %>%
      filter(
        cluster.id %in% cluster.schools.ids,
        str_to_upper(county) == str_to_upper(cluster$County)  # or str_to_upper(...)
      )
    
    # If none match, return FALSE
    if (nrow(subset_schools_sf) == 0) {
      return(FALSE)
    }

    # c) Intersect the matching schools with the cluster geometry
    #    (assuming 'cluster' is an sf with exactly one row)
    intersected_sf <- st_intersection(subset_schools_sf, st_transform(cluster, st_crs(subset_schools_sf)))

    # If intersection is empty, return FALSE
    if (nrow(intersected_sf) == 0) {
      return(FALSE)
    }

    # d) Calculate each intersected area, check if it's >= min.area.frac * school.pop.area
    area_vec <- st_area(intersected_sf)
    min_area <- min.area.frac * school.pop.area
    units(min_area) <- "m^2"
    meets_area_threshold <- area_vec >= min_area

    # e) Count how many intersected polygons meet the threshold and compare to min.num.schools
    meets_count <- sum(meets_area_threshold)
    meets_count >= min.num.schools
  }
}


generate.cluster.schools.pop.tester <- function(primary.schools, min.pop.prop = 2, school.buffer.radius = 500) {
  inner.primary.schools <- buffer.clusters(primary.schools, school.buffer.radius) %>% 
    gUnaryUnion
  
  school.pop.area <- pi * (school.buffer.radius^2)
  
  function(clusters, ...) {
    gIntersection(clusters, inner.primary.schools, byid = TRUE) %>% {
      if (is.null(.)) {
        return(FALSE)
      }
      
      gArea(., byid = TRUE) %>% 
        magrittr::is_weakly_greater_than(min.pop.prop * school.pop.area)
    }
  }
}

ggplot.clusters <- function(selected.clusters, 
                            .rct.schools.data = rct.schools.data, 
                            pilot.locations = NULL, 
                            proj4 = kenya.proj4, 
                            maptype = "roadmap", 
                            source = "google",
                            include.cluster.ids = TRUE, 
                            suppress.selected.clusters = FALSE,
                            caption = NULL,
                            ...) {
  all.clusters <- gUnaryUnion(selected.clusters) 
 
  if (!is.null(.rct.schools.data)) {
    targetable.schools <- .rct.schools.data %>% 
      spTransform(proj4) %>% 
      gWithin(all.clusters, byid = TRUE) %>% 
      drop
  } else {
    targetable.schools <- NULL
  }
  
  all.clusters %<>% spTransform(wgs.84)
  
  map.obj <- all.clusters %>% 
    tidy %>% 
    make_bbox(long, lat, data = .) %>% 
      # get_map(maptype = maptype, source = source, ...) %>% 
      # get_openstreetmap(...) %>% 
      get_stamenmap(maptype = maptype, ...) %>% 
      ggmap()

  if (!suppress.selected.clusters) {
    cluster_polygons <- selected.clusters |>
      spTransform(wgs.84) |>  
      sf::st_as_sf()
    
    map.obj <- map.obj + geom_sf(aes(color = county), alpha = 0.5, linetype = "dashed", inherit.aes = FALSE, data = cluster_polygons)
  }
  
  cluster_center <- spTransform(rct.schools.data, wgs.84)
  
  if (!suppress.selected.clusters) {
    cluster_center %<>% magrittr::extract(., .$cluster.id %in% selected.clusters$cluster.id,)
    
    if (include.cluster.ids) {
      map.obj <- map.obj + geom_text(aes(lon, lat, label = cluster.id), data = as.data.frame(cluster_center))
    } else {
      map.obj <- map.obj + geom_sf(inherit.aes = FALSE, shape = 3, data = sf::st_as_sf(cluster_center)) 
    }
  } else {
    map.obj <- map.obj + geom_sf(inherit.aes = FALSE, shape = 3, color = alpha("darkred", 0.5), data = sf::st_as_sf(cluster_center))
  }
  
  map.obj <- map.obj + geom_point(aes(lon, lat), shape = 3, size = 1, data = pilot.locations) +
    labs(x = "", y = "", caption = caption) +
    scale_color_discrete("County") +
    theme(legend.position = "bottom")
  
  return(map.obj)
}

plot.clusters <- function(selected.clusters, selected.buffers, cluster.ids = FALSE, school.radius = 500) {
  all.clusters <- gUnaryUnion(selected.clusters)
  
  rct.schools.data %>% 
    spTransform(CRS(kenya.proj4)) %>% {
      # buffer.clusters(., school.radius) %>% {
        # magrittr::extract(., !.$cluster.id %in% selected.clusters$cluster.id, ) %>% {
          magrittr::extract(., !drop(gWithin(., all.clusters, byid = TRUE)), ) %>% 
            buffer.clusters(school.radius) %>% 
            plot(border="grey")
          
          magrittr::extract(., drop(gWithin(., all.clusters, byid = TRUE)), ) %>%
            buffer.clusters(school.radius) %>% 
            plot(add = TRUE, border="blue")
      #   }
      # }
      
      magrittr::extract(., .$cluster.id %in% selected.clusters$cluster.id, ) %>% { 
        if (!cluster.ids) plot(., add = TRUE, col = "red")
        if (cluster.ids) text(coordinates(.), labels = .$cluster.id, col = alpha("blue", 0.8), lwd = 2)
      } 
    }

  spTransform(subcounties.adm.data, CRS(kenya.proj4)) %>% 
    gUnaryUnion %>% 
    gIntersection(selected.clusters, byid = TRUE) %>% 
    plot(add = TRUE, col=alpha("red", 0.2))
  
  if (!missing(selected.buffers)) {
    selected.buffers %>% 
      plot(add = TRUE, lty = "dotted", border = alpha("black", 0.5))
  }
}

tu.plot.clusters <- function(yy, ...) {
  yy %>%
    magrittr::extract(.$selected, ) %>%
    plot.clusters(...)
  
  spTransform(rct.schools.data, CRS(kenya.proj4)) %>% {
    schools.data <- .
    purrr::pmap(list(c("bracelet.airtime", "control.ink"), 
                     # c(4000, 2500),
                     attr(yy, "ca.outer.radius")[c("bracelet.airtime", "control.ink")],
                     c("darkgreen", "red")),
                function(grp, radius, color) { 
                  try(schools.data[schools.data$cluster.id %in% yy$cluster.id[yy$cluster.group == grp], ] %>%
                        buffer.clusters(radius) %>%
                        plot(add = TRUE, lty = "dotted", border = alpha(color, 1)), silent = TRUE) 
                })
  }
  
  invisible(NULL)
}

anon_cluster_size_tester <- function(cluster, ...) cluster %>% gArea(byid = TRUE) %>% is_weakly_greater_than(pi * (5000^2) * 0.5)

anon_cluster_size_tester_sf <- function(cluster_sf, ...) {
  # st_area() returns area for each feature in 'cluster_sf'
  # Compare each area to half the area of a circle with radius=5000
  st_area(cluster_sf) >= pi * (5000^2) * 0.5
}

get.rct.clusters <- function(sp.pot.data, rct.area, schools.data, 
                             num.clusters = NULL, # c(control = 30, social.incentive = 60, airtime = 30)
                             ca.outer.radius = 5000, # c(control = 2000, social.incentive = 4000, airtime = 5000) 
                             ca.inner.radius = ca.outer.radius, # c(control = 2000, social.incentive = 2000, airtime = 4000) 
                             cluster.group.order = c("random", "desc.outer.radius"),
                             #school.radius = 500,
                             cluster.size.tester = anon_cluster_size_tester,
                             proj4 = kenya.proj4,
                             verify.shapes = FALSE,
                             plot.rct.clusters = FALSE, seed = 393) {
  set.seed(seed)
  # Unique IDs
  if (is.null(sp.pot.data$cluster.id)) {
    sp.pot.data$cluster.id <- factor(seq_len(nrow(sp.pot.data)))
  }
  
  rand.cluster.group.order <- switch(match.arg(cluster.group.order), random = TRUE, desc.outer.radius = FALSE)
  
  stopifnot(length(dplyr::setdiff(names(ca.outer.radius), names(num.clusters))) == 0)
  stopifnot(all(ca.outer.radius >= ca.inner.radius))
  
  ca.outer.radius %<>% sort(decreasing = TRUE)
  num.clusters  <- num.clusters[names(ca.outer.radius)]
  
  create.buffers <- . %>% 
    purrr::map(~ buffer.clusters(sp.pot.data, ., proj4)) %>% 
    purrr::map(function(buffers) { rownames(buffers) <- buffers$cluster.id; buffers })
  
  outer.rct.buffers <- create.buffers(ca.outer.radius)

  inner.rct.buffers <- create.buffers(ca.inner.radius)[[1]]  
  
  # polygon.ids <- . %>% `@`(polygons) %>% purrr::map(~ .@ID) %>% unlist

  polygon.ids <- . %>% {\(x) row.names(x)}()

  
  # Inner function for possible splitting by (sub)counties 
  inner.get.rct.clusters <- function(inner.buffers, outer.buffers) {
    
    if (plot.rct.clusters && !is.null(rct.area)) {
      rct.area %>%
        st_as_sf() %>%
        st_transform(kenya.proj4) %>%
        plot()
    }
    
    prep.inner.buffers <- function(.inner.buffers) {
      # avail.area <- buffer.clusters(schools.data, school.radius) %>% 
      #   gUnaryUnion %>% 
      #   gIntersection(rct.area) 
      
      if (!is.null(rct.area)) {
        avail.area <- rct.area
        
        .inner.buffers <- .inner.buffers[c(gIntersects(avail.area, .inner.buffers, byid = TRUE)), ]
        .inner.buffers@polygons <- gIntersection(.inner.buffers, avail.area, byid = TRUE)@polygons 
        
        stopifnot(gIsValid(.inner.buffers))
      }
      
      .inner.buffers$selected <- FALSE
      .inner.buffers$cluster.group <- NA
      
      return(.inner.buffers)
    }
   
    inner.buffers %<>% prep.inner.buffers 
    
    rct.clusters.history <- NULL
    dropped.clusters.history <- NULL
    cluster.index <- num.clusters
    
    inner.buffers.intersect.mat = t(st_intersects(outer.buffers[[1]], inner.buffers, sparse = FALSE))
    inner.buffers.intersect.mat <- inner.buffers.intersect.mat[outer.buffers[[1]]$cluster.id %in% inner.buffers$cluster.id, ]
    colnames(inner.buffers.intersect.mat) <- inner.buffers$cluster.id
    rownames(inner.buffers.intersect.mat) <- inner.buffers$cluster.id
    diag(inner.buffers.intersect.mat) <- FALSE
    if (!is.null(schools.data)) {

      contained.schools.mat <- schools.data %>% 
        st_as_sf() %>%
        st_transform(kenya.proj4) %>%
        st_intersects(
          inner.buffers,
          sparse = FALSE
        )  %>%
        magrittr::extract(schools.data$cluster.id %in% inner.buffers$cluster.id, )

      colnames(contained.schools.mat) <- inner.buffers$cluster.id
      # colnames(contained.schools.mat) <- schools.data$cluster.id
      rownames(contained.schools.mat) <- inner.buffers$cluster.id
    } else {
      contained.schools.mat <- NULL
    }
    
    while (is.null(num.clusters) || (sum(cluster.index) > 0)) {
      print("Starting new cluster")
      print(str_glue("Seed: {seed}, Cluster index: {sum(cluster.index)}"))
     
      if (is.null(cluster.index)) {
        current.cluster.group <- NULL
      } else if (rand.cluster.group.order) {
        current.cluster.group <- sample(length(cluster.index), 1, prob = cluster.index * ca.outer.radius)
      } else {
        current.cluster.group <- min(which(cluster.index > 0)) 
      }

      new_clust_options = inner.buffers %>%
        magrittr::extract(!.$selected, ) %>% 
        pull(cluster.id)

      if (length(new_clust_options) > 1) {
        new.rct.cluster = sample(
          new_clust_options,
          1
        )
      }  
      if (length(new_clust_options) == 1) {
        new.rct.cluster = new_clust_options[1]
      } 

      if (length(new_clust_options) == 0) {
        new.rct.cluster = NULL
      } 
      
      if (is.null(new.rct.cluster)) {
        break
      }
    
      new.cluster <- inner.buffers %>% 
        magrittr::extract(.$cluster.id == new.rct.cluster, )
      cl_test = cluster.size.tester(new.cluster, colnames(contained.schools.mat)[contained.schools.mat[new.rct.cluster,]])  

      if (!cl_test) {
        inner.buffers %<>% magrittr::extract(.$cluster.id != new.rct.cluster, )
        
        dropped.clusters.history %<>% c(., new.rct.cluster)  
        
        next  
      }
      
      original.new.cluster <- outer.buffers[[if (!is.null(current.cluster.group)) current.cluster.group else 1]] %>% 
        magrittr::extract(.$cluster.id == new.rct.cluster, )
      
      stopifnot(!is.null(original.new.cluster))
    
      intersected.inner.buffers <- colnames(inner.buffers.intersect.mat)[inner.buffers.intersect.mat[new.rct.cluster, ]] %>% 
        dplyr::intersect(inner.buffers$cluster.id) #%>% 
        # dplyr::setdiff(new.rct.cluster)
      if (any(inner.buffers$cluster.id %in% intersected.inner.buffers)) {

        old.inner.polygons <- st_geometry(
          inner.buffers[inner.buffers$cluster.id %in% intersected.inner.buffers, ])
        saved.intersect.ids <- rlang::duplicate(intersected.inner.buffers)
        

        differenced.polygons <- inner.buffers %>%
          st_as_sf() %>%
          filter(cluster.id %in% intersected.inner.buffers) %>%
          st_difference(st_union(st_as_sf(original.new.cluster)))
        
        if (!is.null(differenced.polygons)) {
          differenced.polygon.ids <- polygon.ids(differenced.polygons)
          # Match on cluster.id to ensure correct ordering:
          idx <- which(inner.buffers$cluster.id %in% differenced.polygon.ids)
          match_ids <- match(inner.buffers$cluster.id[idx], differenced.polygons$cluster.id)

          st_geometry(inner.buffers)[idx] <- st_geometry(differenced.polygons)[match_ids]

          buffers.to.remove <- dplyr::setdiff(intersected.inner.buffers, differenced.polygon.ids)
        } else {
          buffers.to.remove <- intersected.inner.buffers
        }
      }


        # Define the helper function
        check_population_loss <- function(buffers) {
          # Iterate over each row (each selected buffer)
          for (i in seq_len(nrow(buffers))) {
            # Get the current cluster id
            cluster_id <- buffers$cluster.id[i]
            # Determine which columns (i.e. schools) are intersecting for this cluster
            cols_to_test <- colnames(contained.schools.mat)[inner.buffers.intersect.mat[cluster_id, ]]
            # If the cluster size tester returns FALSE, then population is lost.
            if (!cluster.size.tester(buffers[i, ], cols_to_test)) {
              return(TRUE)
            }
          }
          return(FALSE)
        }

        selected_buffers = inner.buffers %>%
          filter(
            selected == TRUE & cluster.id %in% intersected.inner.buffers
          )

        # Extract the relevant subset of inner.buffers
        # selected_buffers <- magrittr::extract(inner.buffers,
        #                                       inner.buffers$selected & (inner.buffers$cluster.id %in% intersected.inner.buffers))

        if (nrow(selected_buffers) == 0) {
          any.rct.population.lost <- FALSE
        } else {
          # Call the function
          any.rct.population.lost <- check_population_loss(selected_buffers)
        }




      # Check if the new cluster would disqualify any of the already selected clusters 
      if (any(rct.clusters.history %in% buffers.to.remove) || any.rct.population.lost) { 
        idx <- inner.buffers$cluster.id %in% intersected.inner.buffers
        st_geometry(inner.buffers)[idx] <- st_geometry(old.inner.polygons)

        
        # Only exclude if there isn't a potentially smaller cluster in other smaller out radius cluster groups
        if (is.null(cluster.index) || max(which(cluster.index > 0)) <= current.cluster.group) { 
          inner.buffers %<>% magrittr::extract(.$cluster.id != new.rct.cluster, )
          
          dropped.clusters.history %<>% c(., new.rct.cluster)  
        }
        
        next  
      } 
      
      inner.buffers$selected[inner.buffers$cluster.id == new.rct.cluster] <- TRUE
      
      inner.buffers <- inner.buffers[!inner.buffers$cluster.id %in% buffers.to.remove, ]
      dropped.clusters.history %<>% c(., buffers.to.remove)
      if (!is.null(current.cluster.group)) {
        inner.buffers$cluster.group[inner.buffers$cluster.id == new.rct.cluster] <- names(cluster.index)[current.cluster.group]
        cluster.index[current.cluster.group] %<>% subtract(1)
      }
      
      if (plot.rct.clusters) {
        inner.buffers %>% magrittr::extract(.$selected, ) %>% plot(add = TRUE, col = scales::alpha("grey", 0.2))
        plot(new.cluster, add = TRUE, border = "red", col = scales::alpha("red", 0.2))
      }
     
      rct.clusters.history %<>% c(., new.rct.cluster)

    }
    
    stopifnot(!anyDuplicated(rct.clusters.history))
    
    attr(inner.buffers, "ca.inner.radius") <- ca.inner.radius
    attr(inner.buffers, "ca.outer.radius") <- ca.outer.radius
    attr(inner.buffers, "cluster.group.order") <- cluster.group.order
    # attr(inner.buffers, "school.radius") <- school.radius
 
    return(inner.buffers) 
  }
  
  inner.get.rct.clusters(inner.rct.buffers, outer.rct.buffers)
}

get.cluster.villages.data <- function(rct.clusters, school.buffer.radius = 1000, min.area.frac = 0.6, .rct.schools.data = rct.schools.data) {


  school.area <- pi * (school.buffer.radius^2)
  
  all.rct.clusters <- st_union(rct.clusters)
  
  schools.within <- .rct.schools.data %>% 
    st_as_sf() %>%
    st_transform(kenya.proj4) %>%
    st_within(
      rct.clusters,
      sparse = FALSE
    ) %>%
    set_rownames(.rct.schools.data$cluster.id) %>% 
    set_colnames(rct.clusters$cluster.id) %>% 
    magrittr::extract(rowSums(.) > 0, )
  
  pot.pts <- .rct.schools.data %>% 
    magrittr::extract(.$cluster.id %in% rct.clusters$cluster.id, ) %>% 
    st_as_sf() %>%
    st_transform(kenya.proj4)

  school.regions <- .rct.schools.data %>% 
    magrittr::extract(.$cluster.id %in% rownames(schools.within), ) %>% 
    buffer.clusters(school.buffer.radius) %>% 
    st_intersection(., all.rct.clusters) 

  min_area = min.area.frac * school.area
  units(min_area) <- "m^2"

  ret.data <- school.regions %>% 
    st_area() %>% 
    magrittr::is_weakly_greater_than(min_area) %>% 
    magrittr::extract(., .) %>%
    names %>% 
    magrittr::extract(schools.within, ., ) %>% 
    adply(1, . %>% { data.frame(pot.cluster.id = names(.)[.]) }, .id = "targeted.cluster.id") %>% 
    left_join(rct.clusters[, c("cluster.id", "cluster.group", "County")], by = c("pot.cluster.id" = "cluster.id"))

rct.clusters %>%
  colnames()

    rct.clusters %>%
      select(cluster.id, cluster.group, county)
    
    
    left_join(.rct.schools.data[, c("cluster.id", "county")], by = c("targeted.cluster.id" = "cluster.id")) %>% 
    filter(county.x == county.y) %>% 
    select(-starts_with("county"))
  
  dist.mat <- spTransform(.rct.schools.data, kenya.proj4) %>% 
    magrittr::extract(.$cluster.id %in% ret.data$targeted.cluster.id,) %>% { 
      set_rownames(set_colnames(gDistance(., pot.pts, byid = T), .$cluster.id), pot.pts$cluster.id)
    }
  
  ret.data %>% 
    group_by(pot.cluster.id) %>% 
    mutate(dist = dist.mat[first(pot.cluster.id), as.character(targeted.cluster.id)],
           village.dist.cat = (ifelse(dist <= (2500/2), "close", "far")),
           cluster.dist.cat = ifelse(all(village.dist.cat == "close"), "close", ifelse(all(village.dist.cat == "far"), "far", "mixed"))) %>% 
    ungroup
}
