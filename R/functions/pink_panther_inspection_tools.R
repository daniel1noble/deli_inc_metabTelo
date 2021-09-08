#' Convert mice output to a long-formatted data frame
#'
#' @param df imputed data output from mice 
#'
#' @return
#' @export
#'
#' @examples
mice_imp_long_format<-function(df){
  mice::complete(df,"long", inc = TRUE) %>% 
    rename(imp=.imp, id=.id) %>% 
    mutate(imp=as.factor(imp))
  
}

#' Inspecting the distribution of original and imputed data: Density plot
#'
#' @param df imputed data output from mice in long format 
#' @param trait trait name
#' @param log_trans logical, log10 transformation
#' @param data_frac a decimal between 0 and 1
#'
#' @return
#' @export
#'
#' @examples imp_density_plot(df=completedData, trait=leaf_area, log_trans=T, data_frac = 0.5)
imp_density_plot <- function(df, trait, log_trans = F, data_frac=0.2) {
  df %<>% group_by(., imp) %>%  sample_frac(data_frac)
  
  dplot <- ggplot() +
    geom_density(
      data = df %>% filter(imp != 0),
      aes(x = {{trait}}, factor = imp),
      color = "lightblue",
      alpha = 0.8) +
    geom_density(
      data = df %>% filter(imp == 0),
      aes(x = {{trait}}),
      color = "red",
      size = 1) +
    theme_minimal() +
    geom_hline(color = "white", yintercept = 0)
  
  if (log_trans == T) {
    dplot + scale_x_log10()
  }
  return(dplot)
}

rphylo_imp_density_plot <- function(df, trait, log_trans = F, data_frac=0.2) {
  df %<>% group_by(., data_id) %>%  sample_frac(data_frac)
  
  dplot <- ggplot() +
    geom_density(
      data = df %>% filter(imp != 0),
      aes(x = {{trait}}, group = data_id),
      color = "lightblue",
      alpha = 0.1) +
    geom_density(
      data = df %>% filter(imp == 0),
      aes(x = {{trait}}),
      color = "red",
      size = 1) +
    theme_minimal() +
    geom_hline(color = "white", yintercept = 0)
  
  if (log_trans == T) {
    dplot + scale_x_log10()
  }
  return(dplot)
}


imp_scatter_plot <- function(df, trait_x, trait_y, log_trans = F, data_frac=0.2) {
  df %<>% group_by(., imp) %>%  sample_frac(data_frac)
  
  dplot <- ggplot() +
    geom_point(
      data = df %>% filter(imp != 0),
      aes(x = {{trait_x}}, y = {{trait_y}}),
      color = "lightblue",
      alpha = 1) +
    geom_line(
      data = df %>% filter(imp != 0),
      aes(x = {{trait_x}}, y = {{trait_y}}),
      color = "lightblue",
      stat="smooth", 
      method = "lm", 
      alpha = 0.8
    ) + 
  geom_point(
    data = df %>% filter(imp == 0),
    aes(x = {{trait_x}}, y = {{trait_y}}),
    color = "red",
    alpha = 0.5) +
    theme_minimal() + 
    geom_line(
      data = df %>% filter(imp == 0),
      aes(x = {{trait_x}}, y = {{trait_y}}),
      color = "red",
      stat="smooth", 
      method = "lm", 
      alpha = 0.8
    )
  
  if (log_trans == T) {
    dplot + scale_x_log10() + scale_y_log10()
  }
  return(dplot)
}
