#' Used to format 3D data
#'
#' @description Function used to obtain specific object, e.g. celltype or gene, and bind values with
#' transformed coordinates.
#' @param df Data.frame with object values, e.g. xCell scores, gene expression, topics, ...
#' @param trans Transformed x, y coordinates centered around (0, 0).
#' @param object Specific object to analyze, e.g. celltype, gene, or topic.
#' @return Object specific values and coordinates.
format_3D_data <- function(df, trans, object){
  object.values = df[rownames(df) == object, ]
  object.values = object.values[colnames(df) %in% rownames(trans)]
  df = df[, colnames(df) %in% rownames(trans)]
  # Filter out objects (genes/cells/topics) with 0 counts in df
  df = df[rowSums(df) != 0,]
  trans = trans[rownames(trans) %in% colnames(df), ]
  # bind spots with object values
  object.values = cbind(trans, object.values)
  return(object.values)
}

#' Scale to custom range
#'
#' @description Scale numeric vector to fit a certain interval.
#' @param x Input numeric vector.
#' @param a Lower limit.
#' @param b Upper limit.
#' @return Scaled numeric vector.
#' @examples
#' library(STanalysis3D)
#'
#' # Scale vector to fit an interval between 0 and 1
#' x <- seq(1, 100)
#' scaled.x <- scale2range(x, a = 0, b = 1)
#' @export
scale2range <- function(x, a, b){
  res <- (b-a)*(x - min(x))/(max(x) - min(x)) + a
  return(res)
}

#' Plot interactive 3D scatter
#'
#' @description Function used to create a 3D scatter of ST data using consecutive tissue sections.
#' @param object Gene name, celltype, topic etc. specifying the object you want to visualize. This will depend on the
#' content of exp.list which holds the values you want to visualize. "object" must be located in the row names of the
#' exp.list data.frames.
#' @param exp.list List of data.frames containing values of interest, e.g. gene expression values, topic proportions,
#' clusters, xCell scores etc.
#' @param trans.list List with transformed x, y coordinates for each replicate section.
#' @param section.input.list List of section.input objects, i.e. data.frames with dot pattern x, y coordinates and
#' corresponding grid cells.
#' @param max Integer value specifying maximum value of colorscale.
#' @param colorscale Set colorscale for 3D scatter (e.g. Jet, Portland, Viridis and grey, see plotly documentation)
#' @param background.col Sets background color using rgb channels. Default 'rgb(0, 0, 0)' = "black". Use function col2rgb("color") to
#' obtain rgb channel values for custom colors ("white" = 'rgb(255, 255, 255)).
#' @param marker.size Numeric value specifying marker size (point size).
#' @param range Set range of colorscale. range = "zero_to_one" will scale values to [0, 1] interval.
#' @param opacity Numeric between 0-1 specifying marker opacity.
#' @param show.top Numeric between 0-1 specifying which subset of data to show. show.top = 0.5 will show all markers with
#' values larger than 50*max(value).
#' @param scalefactor Integer which adds empty distance outside tissue section along z axis. Larger value will make sections
#' appear closer to each other.
#' @param showtickx Logical specifying whther or not to show x axis values.
#' @param showticky Logical specifying whther or not to show y axis values.
#' @param showtickz Logical specifying whther or not to show z axis values.
#' @param add_sections Repeat tissue sections 3 times to increase 3D effect. NOTE: 60% of the data
#' points in each section will be removed to improve performance
#' @param size.scale.factor use value parameter (e.g. xCell score) to set marker size for each dot. size.scale.factor
#' is used as a multiplier, i.e. {markersize = value*scalefactor} for each datapoint. For example, xCell scores
#' ranging between 0 -> 0.2 can be scaled up to a maximum marker size of 1 by using a scalefactor of 5. Scaling
#' the dots this way will make datapoints with low values less visible and make the high value volumes more distinct
#' in the scatter plot. This scaling can be used in combination with show.top to clean "less important" space from
#' data points which are concealing the high value clusters.
#' @param cam Set camera angle of plot object.
#' @param skip.every Used to subset final data output.
#' @param return.df Logical specifying whether or not the data should be returned. This will override the plot functionality.
#' @importFrom plotly plot_ly add_markers layout
#' @return Interactive 3D scatter object.
#' @export
morph.3d = function(object,
                              exp.list,
                              trans.list,
                              section.input.list,
                              max,
                              colorscale = "Jet",
                              background.col = 'rgb(0, 0, 0)',
                              grid.color = 'rgb(0, 0, 0)',
                              marker.size = 1.5,
                              min.size = 0,
                              max.size = 1.5,
                              range = NULL,
                              opacity = 0.1,
                              show.top = NULL,
                              scalefactor = 0,
                              showtickx = FALSE,
                              showticky = FALSE,
                              showtickz = FALSE,
                              add_sections = FALSE,
                              size.scale.factor = 0,
                              cam = NULL,
                              skip.every = 3,
                              return.df = FALSE,
                              subsample = NULL,
                              add.rotation = TRUE,
                              rot.angle = 0.01,
                              rot.offset = 20){

  object.list <- list()
  for (i in 1:length(exp.list)){
    norm.exp.values <- exp.list[[i]]
    trans <- trans.list[[i]]
    object.list[[i]] <- format_3D_data(norm.exp.values, trans, object)
  }

  df.list <- list()
  for (i in 1:length(object.list)){
    x <- as.integer(colnames(section.input.list[[i]])[1])
    y <- as.integer(colnames(section.input.list[[i]])[2])
    section <- section.input.list[[i]]
    section[is.na(section)] <- 0
    df.list[[i]] <- interpolate_2D_data(object.list[[i]], section, x, y, iterval = i)
  }
  if (is.null(range)){
    df.final <- na.omit(do.call(rbind, df.list))
  } else if (range == "zero_to_one"){
    df.final <- na.omit(do.call(rbind, df.list))
    df.final$val <- scale2range(df.final$val, a = 0, b = 1)
    max = 1
  } else if (!is.null(range)){
    print("Range not allowed")
    return(NULL)
  }

  if (!is.null(show.top)){
    df.final <- subset(df.final, (val >= (max(df.final$val)*show.top)))
    df.final$val <- scale2range(df.final$val, 0, 1)
    max <- 1
  }

  if (add_sections){
    df2 <- df.final
    df2$z <- df2$z + 0.33
    df3 <- df.final
    df3$z <- df3$z + 0.67
    df.final <- rbind(df.final, df2)
    df.final <- rbind(df.final, df3)
    df.final <- df.final[seq(1, nrow(df.final), skip.every), ]
    if (return.df) {
      return(df.final)
    }
  }

  if (return.df & !is.null(subsample)){
    df.final <- df.final[sample(1:nrow(df.final), subsample), ]
    return(df.final)
  } else if (return.df) {
    return(df.final)
  }

  if (!is.null(subsample)) {
    df.final <- df.final[sample(1:nrow(df.final), subsample), ]
  }

  if (!(size.scale.factor %in% c(-1, 0))){
    marker.size = ~val*size.scale.factor + min.size
  }

  if (size.scale.factor == -1){
    marker.size = ~scale2range(val, min.size, max.size)
  }

  if (is.null(cam)){
    cam = list(eye = list(x = -1.25, y = 1.25, z = 1))
  } else {
    cam = cam
  }

  p <- plot_ly(df.final,
               x = ~x, y = ~y, z = ~z,
               marker = list(color = ~scale2range(val, 0, 1),
                             cmin = 0,
                             cmax = max,
                             showscale = T,
                             colorscale = colorscale,
                             size = marker.size,
                             opacity = opacity)) %>%
    add_markers() %>%
    layout(scene = list(xaxis = list(range = c(min(min(df.final[, 1], na.rm = T), min(df.final[, 1], na.rm = T)), max(max(df.final[, 1], na.rm = T), max(df.final[, 2], na.rm = T))), title = 'x', showticklabels = showtickx, gridcolor = grid.color, zerolinecolor = grid.color),
                        yaxis = list(range = c(min(min(df.final[, 1], na.rm = T), min(df.final[, 1], na.rm = T)), max(max(df.final[, 1], na.rm = T), max(df.final[, 2], na.rm = T))), title = 'y', showticklabels = showticky, gridcolor = grid.color, zerolinecolor = grid.color),
                        zaxis = list(title = 'z', showticklabels = showtickz, range = c(-1 -scalefactor, 6 + scalefactor), gridcolor = grid.color, zerolinecolor = grid.color),
                        camera = cam)) %>%
    layout(paper_bgcolor = background.col)
  if (add.rotation) {
    javascript <- HTML(paste('<div id="LiveChat_1308239999" style="top:70%%;left:0px;position: absolute;">
                             <div class="btn-group" role="group" aria-label="Basic example">
                             <button id="playpause" class="btn btn-secondary">Rotate</button>
                             </div>
                             </div>
                             <script>
                             var id = document.getElementById("htmlwidget_container").firstChild.nextSibling.id
	gd = document.getElementById(id)
                             var rot_angle = ', rot.angle, '
                             i = 0
                             rotation = false
                             function timeout() {
                             if (rotation) {
                             rotate(rot_angle)
                             i += 1
                             }
                             if (i < 100000)
                             setTimeout(timeout, ', rot.offset, ')
                             }

                             function rotate(angle) {
                             var scene = gd._fullLayout["scene"];
                             var camera = scene.camera;

                             var rtz = xyz2rtz(camera.eye);
                             rtz.t += angle;
                             camera.eye = rtz2xyz(rtz);

                             Plotly.relayout(gd, "scene.camera", camera);
                             }

                             function rtz2xyz(rtz) {
                             return {
                             x: rtz.r * Math.cos(rtz.t),
                             y: rtz.r * Math.sin(rtz.t),
                             z: rtz.z
                             };
                             }

                             function xyz2rtz(xyz) {
                             return {
                             r: Math.sqrt(xyz.x * xyz.x + xyz.y * xyz.y),
                             t: Math.atan2(xyz.y, xyz.x),
                             z: xyz.z
                             };
                             }

                             $("#playpause").click(function() {rotation = !rotation})
                             $(document).ready(timeout)
                             </script>
                             '), sep = '')
    p <- appendContent(p, javascript)
  }
  return(p)
}

#' Create interactive scatterplot3js
#'
#' @description Given a data.frame with x, y, z coordinates and interpolates values, this function
#' can be used to generate a html widget visualizing 3D scatter data.
#' @param df data.frame with x, y, z coordinates and interpolated values.
#' @param colorscale specify which colorscale to use. Options: "viridis", "magma", "plasma", "inferno" or
#' an RColorBrewer palette
#' @param size Set marker size in rendered plot
#' @param offset Time offset between frames in milliseconds.
#' @param subsample Integer specifying number of data points to include. Downsampling will
#' decresase rendering time.
#' @return html widget
#' @importFrom viridis viridis magma inferno plasma
#' @importFrom threejs scatterplot3js
#' @importFrom htmlwidgets prependContent
#' @importFrom htmltools HTML
#' @importFrom RColorBrewer brewer.pal
#' @export
js3Dscatter <- function(df, filepath1 = NULL, colorscale = "viridis", size = 0.1, zlim = c(-30, 36), xlim = NULL, ylim = NULL, offset = 20, subsample = NULL, background.col = "#ffffff", zoom = 3, z = 2, neg = FALSE, alpha = 0.5, use.orbitcontrols = TRUE) {
  if (is.null(xlim)) {
    xlim = c(min(df[, 1]) - 5, max(df[, 1]) + 5)
    ylim = xlim
  }
  if (!is.null(subsample)) {
    df <- df[sample(1:nrow(df), subsample), ]
  }
  stopifnot(class(df) == "data.frame",
            class(colorscale) == "character")
  if (colorscale == "viridis") {
    colors <- rep(viridis::viridis(n = 100, alpha = alpha), each = round(length(df[, 4])/100))
  } else if (colorscale == "magma") {
    colors <- rep(viridis::magma(n = 100, alpha = alpha), each = round(length(df[, 4])/100))
  } else if (colorscale == "inferno") {
    colors <- rep(viridis::inferno(n = 100, alpha = alpha), each = round(length(df[, 4])/100))
  } else if (colorscale == "plasma") {
    colors <- rep(viridis::plasma(n = 100, alpha = alpha), each = round(length(df[, 4])/100))
  } else if (colorscale == "heat") {
    colors <- rep(heat.colors(n = 100, alpha = alpha), each = round(length(df[, 4])/100))
  } else if (colorscale == "jet") {
    jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    colors <- rep(jet.colors(n = 100), each = round(length(df[, 4])/100))
    colors <- addalpha(colors, alpha)
  } else {
    colors <- rep(colorRampPalette(brewer.pal(11, colorscale))(100), each = round(length(df[, 4])/100))
    colors <- addalpha(colors, alpha)
  }
  if (length(colors) > length(df[, 4])) {
    colors <- colors[1:length(df[, 4])]
  } else if (length(df[, 4]) > length(colors)) {
    ext <- length(df[, 4]) - length(colors)
    colors <- c(colors, rep(colors[length(colors)], ext))
  }
  colors = colors[(1:length(colors))[rank(df[, 4])]]
  merged.df <- data.frame(df, colors)

  # use alt. coloring
  # color.df <- data.frame(val = sort(unique(df[, 4])), col = viridis::viridis(n = length(unique(df[, 4])), alpha = 0.5))
  # merged.df <- merge(df, color.df, by = "val")
  if (neg) {
    y = -merged.df[, 2]
  } else {
    y = merged.df[, 2]
  }
  p <- threejs::scatterplot3js(x = merged.df[, 1],
                      y = y,
                      z = merged.df[, 3],
                      color = merged.df[, 5],
                      size = size,
                      axis = FALSE,
                      grid = FALSE,
                      zlim = zlim,
                      xlim = xlim,
                      ylim = ylim,
                      center = c(2, 2),
                      offset = 40,
                      bg = background.col,
                      use.orbitcontrols = use.orbitcontrols)
  if (!is.null(filepath1)) {
    txt1 <- base64encode(readBin(filepath1, "raw", file.info(filepath1)[1, "size"]), "txt")
    html_str1 <- sprintf('<img src="data:image/png;base64,%s height="20%%" width="20%%" style="top:0px;right:0px;position: absolute;">', txt1)
  }
 javascript <- HTML(paste('<div id="LiveChat_1308239999" style="top:70%%;left:0px;position: absolute;">
                           <div class="btn-group" role="group" aria-label="Basic example">
                           <button id="playpause" class="btn btn-secondary">Rotate</button>
                           <button id="reset" class="btn btn-secondary">Reset</button>
                           </div>
                           </div>
                           <div>
                           <script>
                           i = 1.575
                           zoom = ', zoom, '
                           rotate = false
                           function timeout() {
                           if (rotate) {
                           w = HTMLWidgets.find(".scatterplotThree").widget
                           w.camera.position.set(Math.cos(i)*zoom, ', z, ', Math.sin(i)*zoom)
                           w.controls.update()
                           i += 0.005
                           }
                           if (i < 1000)
                           setTimeout(timeout, ', offset, ')
                           }
                           function reset() {
                           w = HTMLWidgets.find(".scatterplotThree").widget
                           w.camera.position.set(Math.cos(1.575)*zoom, ', z, ', Math.sin(1.575)*zoom)
                           w.controls.update()
                           i = 1.575
                           }
                           $("#playpause").click(function() {rotate = !rotate})
                           $("#reset").click(reset)
                           $(document).ready(reset)
                           $(document).ready(timeout)
                           </script>
                           </div>
                           '), sep = '')
  if (!is.null(filepath1)) {
    javascript <- HTML(paste(javascript, html_str1), sep = '')
  }
  p <- prependContent(p,javascript)
  return(p)
}

#' Add alpha to channel 4 in rgba space
#'
#' @param colors Input colors which will be converted to rgb
#' @param alpha Set constant alpha level for rgb colors
#' @return rbga colors
addalpha <- function(colors, alpha = 1.0) {
  r <- col2rgb(colors, alpha = T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}
