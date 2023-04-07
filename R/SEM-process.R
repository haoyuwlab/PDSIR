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
#' @param xlab label of x-axis
#' @param text_size size of text in figure
#' @param include_legend logical; wehther to include a legend
#'
#'
#' @return save a figure of the compartments' trajectories
#' @export
#'

draw_trajectories <- function(
            SEM, t_end = 10,
            xlab = "time", include_legend = FALSE
            ) {

      legend.position <- if(include_legend){
            "right"
      }else{
            "none"
      }

      t <- SEM[["t"]]
      S <- SEM[["S"]]
      I <- SEM[["I"]]
      X <- SEM[["X"]]
      N <- S[1] + I[1]

      tbl <- tibble::tibble(t = t, Susceptible = S, Infectious = I, Removed = N - S - I, X = X)

      tbl %>%
            tidyr::pivot_longer(
               cols = (.data$Susceptible) : (.data$Removed),
               names_to = "Compartments", values_to = "count"
               ) %>%
            dplyr::filter(.data$t < t_end) %>%
            dplyr::mutate(Compartments = factor(.data$Compartments, levels = c("Susceptible", "Infectious", "Removed"))) %>%
            ggplot2::ggplot(ggplot2::aes(.data$t, .data$count, color = .data$Compartments)) +
            ggplot2::geom_line(size = 1.5) +
            ggplot2::theme(
                  legend.position = legend.position
            ) +
            ggplot2::labs(x = xlab, y = "size of\ncompartment")

}

