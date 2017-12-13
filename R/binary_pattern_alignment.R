#' Scale, move and rotate binary dot pattern
#'
#' @description This function is used to manually align binary dot patterns.
#' @param data A list with a binary dot pattern and corresponding feature coordinates used to create the pattern.
#' @param df_ref A reference binary dot pattern that that you want to align the new pattern to.
#' @param sx Scale x coordinates with sx scalar, i.e. expand or contract along x axis.
#' @param sy Scale y coordinates with sy scalar, i.e. expand or contract along y axis.
#' @param shift_x Move pattern along x axis.
#' @param shift_y Move pattern along y axis.
#' @param flip_x Mirror along x axis.
#' @param fip_y Mirror along y axis.
#' @param alpha Set alpha level for data in plot.
#' @param size Set size of dots in plot.
#' @param keep.ref.fixed Fix the reference dot pattern if you only want to change the coordinates.
#' @return Dataframe or list with scaled data.
manipulate.grid <- function(data,
                            df_ref = NULL,
                            sx = 1,
                            sy = 1,
                            shift_x = 0,
                            shift_y = 0,
                            flip_x = F,
                            flip_y = F,
                            alpha = 0.5,
                            size = 0.5,
                            plot.coords = F,
                            return.df = F,
                            keep.ref.fixed = F)
{
  if (class(data) == "list"){
    df <- data[[1]]
    coords <- data[[2]]
    meanx <- mean(df[, 1])
    # Center x axis around 0
    if (meanx != 0) {
      if (keep.ref.fixed == F) {
        x <- df[, 1] - meanx
        x <- x*sx
      } else {
        x <- df[, 1] - meanx
      }
      if (keep.ref.fixed == T) {
        coords[, 1] <- coords[, 1] - mean(coords[, 1])
        coords[, 1] <- coords[, 1]*sx
      } else {
        coords[, 1] <- coords[, 1] - meanx
        coords[, 1] <- coords[, 1]*sx
      }
    }
    # Move pattern along x axis
    if (shift_x != 0){
      if (keep.ref.fixed == F) {
        x <- x + shift_x
      }
      coords[, 1] <- coords[, 1] + shift_x
    }
    # Flip x axis
    if (flip_x){
      if (keep.ref.fixed == F) {
        x <- -x
      }
      coords[, 1] <- -coords[, 1]
    }
    meany <- mean(df[, 2])
    # Center x axis around 0
    if (meany != 0) {
      if (keep.ref.fixed == F) {
        y <- df[, 2] - meany
        y <- y*sy
      } else {
        y <- df[, 2] - meany
      }
      if (keep.ref.fixed == T) {
        coords[, 2] <- coords[, 2] - mean(coords[, 2])
        coords[, 2] <- coords[, 2]*sy
      } else {
        coords[, 2] <- coords[, 2] - meany
        coords[, 2] <- coords[, 2]*sy
      }
    }
    # Move pattern along y axis
    if (shift_y != 0){
      if (keep.ref.fixed == F) {
        y <- y + shift_y
      }
      coords[, 2] <- coords[, 2] + shift_y
    }
    # Flip y axis
    if (flip_y){
      if (keep.ref.fixed == F) {
        y <- -y
      }
      coords[, 2] <- -coords[, 2]
    }
    # Define new data.frame
    df <- data.frame(x, y)
  }
  # CASE 1: If there's no reference data present and plot.coords = F, only dot pattern will be plotted
  # Any changes in df will be plotted
  if (is.null(df_ref)){
    # Define boundaries
    min.val <- min(min(df$x), min(df$y))
    max.val <- max(max(df$x), max(df$y))
    p <- ggplot2::ggplot(as.data.frame(df), ggplot2::aes(x, y)) +
      ggplot2::scale_x_continuous(limits = c(min.val, max.val)) +
      ggplot2::scale_y_continuous(limits = c(min.val, max.val)) +
      ggplot2::geom_point(stat = "identity", size = 0.5)
    # If plot.coords = T, coordinates are plotted as reactive data
    if (plot.coords){
      p <- p + ggplot2::geom_point(data = as.data.frame(coords), ggplot2::aes(x, y), stat = "identity", col = "red", alpha = alpha, size = size)
    }
    plot(p)
    # Return data if return.df = T
    if (return.df == T){
      return(list(df, coords))
    }
  }
  # CASE 2: If a reference dataset is present, the reference data will be fixed but the reactive df can change
  if (!is.null(df_ref)) {
    if (class(df_ref) == "list") {
      df_ref <- df_ref[[1]]
    }
    # Define boundaries
    min.val <- min(min(df$x), min(df$y), min(df_ref$x), min(df_ref$y))
    max.val <- max(max(df$x), max(df$y), max(df_ref$x), max(df_ref$y))
    colnames(df_ref) <- c("x", "y")
    p <- ggplot2::ggplot(as.data.frame(df_ref), ggplot2::aes(x, y)) +
      ggplot2::geom_point(stat = "identity", size = 0.5) +
      ggplot2::geom_point(data = as.data.frame(df), ggplot2::aes(x, y), stat = "identity", size = size, color = "red", alpha = alpha) +
      ggplot2::scale_x_continuous(limits = c(min.val, max.val)) +
      ggplot2::scale_y_continuous(limits = c(min.val, max.val))
    plot(p)
    if (class(data) == "list" & return.df == T){
      return(list(df, coords))
    }
  }
}

#' Align grid using shiny app
#'
#' @description Align binary dot pattern and calculate normalized coordinates. The alignment is doe using a local shiny app.
#' A reference dot pattern (black color) serves as a template for the alignment of a second dot pattern (red color). The second
#' dot pattern is automatically scaled to center around 0 and can then be rotated, scaled along x and y axis, flipped or moved.
#' Once the dot pattern is properly aligned with the reference dot pattern, the data can be downloaded as an .RData file.
#' @param scatter A dot pattern with x, y values or a list with dot pattern in as the first element and corresponding coordinates
#' as the second element.
#' @param df_ref list with dot pattern in as the first element and corresponding coordinates
#' as the second element used for alignment. The reference pattern can be obtained by running this function
#' without a reference dot pattern.
#' @param plot.coords Logical specifying if feature coordinates should be plotted on top of binary dot pattern.
#' For the first replicate, it is beneficial to set this to true to evaluate how well the features align with the dot pattern.
#' @param alpha Numeric value specifying alpha level of feature spots. Only required if plot.coords = TRUE.
#' @param size Numeric value specifying size of feature spots. Only required if plot.coords = TRUE.
#' @param return.df Logical specifying if the manipulated data should be returned.
#' @param output.name Character specifying the name of the output object.
#' by running the function without df_ref.
#' @return Download normalized data to an .RData file.
#' @export
align.grids <- function(scatter, df_ref = NULL, plot.coords = F, return.df = T, alpha = 0.5, size = 0.5, output.name = "trans.coord.list", keep.ref.fixed = F, scale.coords.ind = F) {
  ui <- shiny::fluidPage(
    shiny::fluidRow(
      shiny::column(3,
             shiny::sliderInput(
               inputId = "angle",
               label = "Rotation angle",
               value = 0, min = -180, max = 180, step = 0.5
             ),
             shiny::sliderInput(
               inputId = "scale_x",
               label = "Scale x axis",
               value = 1, min = 0.1, max = 2, step = 0.01
             ),
             shiny::sliderInput(
               inputId = "scale_y",
               label = "Scale y axis",
               value = 1, min = 0.1, max = 2, step = 0.01
             ),
             shiny::sliderInput(
               inputId = "shift_x",
               label = "Move along x axis",
               value = 0, min = -10, max = 10, step = 0.025
             ),
             shiny::sliderInput(
               inputId = "shift_y",
               label = "Move along y axis",
               value = 0, min = -10, max = 10, step = 0.025
             ),
             shiny::sliderInput(
               inputId = "size",
               label = "Change point size",
               value = 0.5, min = 0.1, max = 6, step = 0.1
             ),
             shiny::checkboxInput(inputId = "flip_x",
                           label = "Mirror along x axis",
                           value = FALSE),
             shiny::checkboxInput(inputId = "flip_y",
                           label = "Mirror along y axis",
                           value = FALSE),
             shiny::textInput(inputId = "output.var", label = "dataset name"),
             shiny::downloadButton('downloadModel', 'Download aligned data', class = "dlButton")
      ),

      shiny::column(8,
                    shiny::plotOutput("scatter")
      )
    )
  )

  server <- function(input, output) {
    getData <- shiny::reactive({
      if (class(scatter) == "list") {
        data <- list()
        if (keep.ref.fixed == F) {
          data[[1]] <- as.data.frame(spdep::Rotation(scatter[[1]], angle = as.numeric(input$angle*0.01745329)))
          colnames(data[[1]]) <- c("x", "y")
        } else {
          data[[1]] <- scatter[[1]]
        }
        data[[2]] <- as.data.frame(spdep::Rotation(scatter[[2]], angle = as.numeric(input$angle*0.01745329)))
        colnames(data[[2]]) <- c("x", "y")
      } else if (class(scatter) %in% c("data.frame", "matrix")) {
        data <- as.data.frame(spdep::Rotation(scatter, angle = as.numeric(input$angle*0.01745329)))
        colnames(data) <- c("x", "y")
      }
      manipulate.grid(data, df_ref = df_ref,
                sx = input$scale_x, sy = input$scale_y,
                shift_x = input$shift_x, shift_y = input$shift_y,
                flip_x = input$flip_x, flip_y = input$flip_y,
                alpha = alpha, size = input$size,
                plot.coords = plot.coords,
                return.df = return.df,
                keep.ref.fixed = keep.ref.fixed)
    })

    output$scatter <- shiny::renderPlot({

      getData()

    }, height = 800, width = 800)

    output$downloadModel <- shiny::downloadHandler(
      filename <- function() {
        paste("data-", Sys.Date(), ".RData", sep="")
      },
      content = function(file) {
        assign(output.name, getData())
        save(list = output.name, file = file)
      })
  }

  shiny::shinyApp(server = server, ui = ui)

}
