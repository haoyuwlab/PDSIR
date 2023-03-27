#' Draw trajectories of a SEM
#'
#' Line plot of the compartment sizes over time.
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
#' @inheritParams simulate_SEM
#'
#' @param SEM list corresponding to stochastic epidemic process
#' @param save_figure logical; whether to save the figure
#' @param path_figure path to save figure
#' @param id_figure logical; whether to save the figure
#' @param xlab label of x-axis
#' @param text_size size of text in figure
#'
#' @return save a figure of the compartments' trajectories
#' @export
#'

draw_trajectories <- function(
            SEM, t_end = 10,
            save_figure = FALSE, path_figure = NULL, id_figure = NULL,
            xlab = "Time", text_size = 30
            ) {

      t <- SEM[["t"]]
      S <- SEM[["S"]]
      I <- SEM[["I"]]
      X <- SEM[["X"]]
      N <- S[1] + I[1]

      tbl <- tibble::tibble(t = t, S = S,        I = I, R = N - S - I    , X = X)

      g <- tbl %>%
            tidyr::pivot_longer(cols = (.data$S) : (.data$R), names_to = "Compartments", values_to = "count") %>%
            dplyr::filter(.data$t < t_end) %>%
            dplyr::mutate(Compartments = factor(.data$Compartments, levels = c("S", "I", "R"))) %>%
            ggplot2::ggplot(ggplot2::aes(.data$t, .data$count, color = .data$Compartments)) +
            ggplot2::geom_line(size = 1.5) +
            ggplot2::theme(
                  text = ggplot2::element_text(size = text_size),
                  legend.position = "none"
            ) +
            ggplot2::labs(x = xlab, y = "Size of compartment")

      print(g)

      if(save_figure){
            ggplot2::ggsave(
                  paste0(id_figure, "_trajectories.jpeg"),
                  path = path_figure, width = 1.61803, height = 1, scale = 5
            )
      }

}

