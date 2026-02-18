# ==============================================================================
# TFG - ANÀLISI DE L'EQUILIBRI COMPETITIU (ATP)
# ==============================================================================
# 1. LLIBRERIES I CONFIGURACIÓ
# ==============================================================================
library(dplyr)
library(mgcv)
library(ggplot2)
library(plotly)
library(pROC)
library(tidyr)
library(ggtext)
library(lubridate)

# --- CÀRREGA DE DADES ---
# NOTA: Cal que el fitxer 'atp_tennis.csv' estiga a la carpeta 'dades', sinò ajusta la ruta. 
fitxer_dades <- "dades/atp_tennis.csv" 

if(file.exists(fitxer_dades)) {
  atp_tennis <- read.csv(fitxer_dades, sep = ",", stringsAsFactors = FALSE)
} else {
  atp_tennis <- read.csv("atp_tennis.csv", sep = ",", stringsAsFactors = FALSE)
}
atp_tennis$Date <- as.Date(atp_tennis$Date, format = "%Y-%m-%d")

# ==============================================================================
# 2. DEFINICIÓ DE FUNCIONS
# ==============================================================================

# --- Funció A: Calcular el model i graficació ---
cross_basis_bySurface <- function(atp_filtrat, k1, k2, period_name) {
  
  gam_6 <- gam(
    Y ~ s(millor_rank, k = k1, bs = "tp", by = Surface) + 
      s(ranking_diff, k = k2, bs = "tp", by = Surface) +
      millor_rank:ranking_diff,
    data   = atp_filtrat,
    family = binomial,
    method = "REML"
  )
  
  pred_probs <- predict(gam_6, type = "response")
  roc_obj <- pROC::roc(atp_filtrat$Y, pred_probs, quiet = TRUE)
  auc_value <- pROC::auc(roc_obj)
  aic_value <- AIC(gam_6)
  
  TOP <- 200
  nombres_superficies <- levels(atp_filtrat$Surface)
  n_rank <- 100 
  n_diff <- 100
  
  grid_df <- expand.grid(
    millor_rank   = seq(1, TOP, length.out = n_rank),
    ranking_diff  = seq(0, TOP, length.out = n_diff),
    Surface       = nombres_superficies
  )
  
  grid_df$Surface <- factor(grid_df$Surface, levels = c("Hard", "Clay", "Grass"))
  
  grid_df <- grid_df %>%
    dplyr::filter(millor_rank + ranking_diff <= TOP)
  
  pr_grid <- predict(gam_6, newdata = grid_df, type = "link", se.fit = TRUE)
  z <- qnorm(0.975)
  
  grid_df$pi_hat <- plogis(pr_grid$fit)
  grid_df$pi_lo  <- plogis(pr_grid$fit - z * pr_grid$se.fit)
  grid_df$pi_hi  <- plogis(pr_grid$fit + z * pr_grid$se.fit)
  
  p <- ggplot(grid_df, aes(x = millor_rank, y = ranking_diff, fill = pi_hat)) +
    geom_tile() +
    scale_fill_viridis_c(option = "viridis", limits = c(0, 1)) +
    facet_wrap(~Surface) + 
    coord_fixed() + 
    scale_x_continuous(limits = c(1, TOP), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, TOP), expand = c(0, 0)) +
    labs(
      x = "Millor rànquing", 
      y = "Diferència de rànquing"
    ) +
    theme_minimal()
  
  print(p)
  
  return(list(gam = gam_6, grid = grid_df, aic = aic_value, auc = auc_value, plot_abs = p))
}


# --- Funció B: Gràfic de diferències ---
plot_surface_diffs <- function(gam_obj, period_name, z_limits = c(-0.25, 0.25)) {
  
  TOP <- 200
  grid_data <- gam_obj$grid
  
  grid_data <- grid_data %>% 
    dplyr::filter(Surface != "Carpet") %>%
    dplyr::mutate(Surface = droplevels(as.factor(Surface)))
  
  grid_wide <- grid_data %>%
    tidyr::pivot_wider(
      id_cols = c(millor_rank, ranking_diff), 
      names_from = Surface, 
      values_from = pi_hat
    )
  
  red_hex  <- "#b2182b"
  blue_hex <- "#2166ac"
  
  label_HC <- paste0("<span style='color:", red_hex, "'>**Hard**</span> respecte <span style='color:", blue_hex, "'>**Clay**</span>")
  label_HG <- paste0("<span style='color:", red_hex, "'>**Hard**</span> respecte <span style='color:", blue_hex, "'>**Grass**</span>")
  label_CG <- paste0("<span style='color:", red_hex, "'>**Clay**</span> respecte <span style='color:", blue_hex, "'>**Grass**</span>")
  
  diff_list <- list()
  surfaces_avail <- names(grid_wide)
  
  if (all(c("Hard", "Clay") %in% surfaces_avail)) {
    diff_list[["HC"]] <- grid_wide %>%
      dplyr::mutate(Diff = Hard - Clay, Comparacio = label_HC) %>%
      dplyr::select(millor_rank, ranking_diff, Comparacio, Diff)
  }
  
  if (all(c("Hard", "Grass") %in% surfaces_avail)) {
    diff_list[["HG"]] <- grid_wide %>%
      dplyr::mutate(Diff = Hard - Grass, Comparacio = label_HG) %>%
      dplyr::select(millor_rank, ranking_diff, Comparacio, Diff)
  }
  
  if (all(c("Clay", "Grass") %in% surfaces_avail)) {
    diff_list[["CG"]] <- grid_wide %>%
      dplyr::mutate(Diff = Clay - Grass, Comparacio = label_CG) %>%
      dplyr::select(millor_rank, ranking_diff, Comparacio, Diff)
  }
  
  if (length(diff_list) == 0) return(NULL)
  
  plot_data <- dplyr::bind_rows(diff_list)
  
  plot_data$Comparacio <- factor(
    plot_data$Comparacio, 
    levels = c(label_HC, label_HG, label_CG)
  )
  
  # --- MÀSCARA TRIANGULAR: mateix domini TOP200 ---
  plot_data <- plot_data %>%
    dplyr::filter(
      ranking_diff >= 0,
      millor_rank >= 1, millor_rank <= TOP,
      millor_rank + ranking_diff <= TOP
    )
  
  p <- ggplot(plot_data, aes(x = millor_rank, y = ranking_diff, fill = Diff)) +
    geom_tile() +
    scale_fill_gradient2(
      low = blue_hex, mid = "#f7f7f7", high = red_hex,
      midpoint = 0, limits = z_limits, oob = scales::squish,
      name = "Diferència de \nProbabilitat"
    ) +
    facet_wrap(~Comparacio) +
    coord_fixed() + 
    scale_x_continuous(limits = c(1, TOP), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, TOP), expand = c(0, 0)) +
    labs(
      x = "Millor rànquing", 
      y = "Diferència de rànquing"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      strip.text = ggtext::element_markdown(size = 12)
    )
  
  return(p)
}

# --- Funció C: funcions auxiliars per extreure corbes 1D---
get_curve_diff_fixed <- function(models_list, val_diff) {
  do.call(rbind, lapply(names(models_list), function(nm) {
    obj <- models_list[[nm]]
    grid_local <- obj$grid
    idx <- which.min(abs(grid_local$ranking_diff - val_diff))
    sel_diff <- grid_local$ranking_diff[idx]
    
    grid_local %>%
      filter(abs(ranking_diff - sel_diff) < 1e-9) %>%
      group_by(millor_rank) %>%
      summarise(pi_hat = mean(pi_hat, na.rm=T), pi_lo = mean(pi_lo, na.rm=T), pi_hi = mean(pi_hi, na.rm=T), .groups="drop") %>%
      mutate(Subconjunt = nm)
  }))
}

get_curve_rank_fixed <- function(models_list, val_rank) {
  do.call(rbind, lapply(names(models_list), function(nm) {
    obj <- models_list[[nm]]
    grid_local <- obj$grid
    idx <- which.min(abs(grid_local$millor_rank - val_rank))
    sel_rank <- grid_local$millor_rank[idx]
    
    grid_local %>%
      filter(abs(millor_rank - sel_rank) < 1e-9) %>%
      group_by(ranking_diff) %>%
      summarise(pi_hat = mean(pi_hat, na.rm=T), pi_lo = mean(pi_lo, na.rm=T), pi_hi = mean(pi_hi, na.rm=T), .groups="drop") %>%
      mutate(Subconjunt = nm)
  }))
}

# Funció per a Interactius
plot_comparison_interactive <- function(data, x_var, title_prefix, val_fixed, x_lab) {
  p <- ggplot(data, aes_string(x = x_var, y = "pi_hat", color = "Subconjunt")) +
    geom_line(size = 1) + 
    geom_ribbon(aes(ymin = pi_lo, ymax = pi_hi, fill = Subconjunt), alpha = 0.15, linewidth = 0) +
    labs(title = paste(title_prefix, val_fixed), x = x_lab, y = "Probabilitat", color = "Subconjunt", fill = "Subconjunt") +
    theme_minimal(base_size = 14) + 
    scale_y_continuous(breaks = seq(0.3, 1.0, by = 0.1)) +
    coord_cartesian(ylim = c(0.25, 1), xlim = c(0, 200), expand = FALSE)
  ggplotly(p, tooltip = c("x", "y", "Subconjunt"))
}

# ==============================================================================
# 3. 'ANÀLISI PER PERÍODES
# ==============================================================================

# --- BLOC 1: 2000-2005 ---
atp_00_05 <- atp_tennis %>%
  filter(year(Date) >= 2000, year(Date) <= 2005, Rank_1 <= 200, Rank_2 <= 200) %>%
  mutate(millor_rank = pmin(Rank_1, Rank_2), Y = if_else(Winner == if_else(Rank_1 < Rank_2, Player_1, Player_2), 1, 0), ranking_diff = abs(Rank_1 - Rank_2)) %>%
  filter(Surface != "Carpet") %>% mutate(Surface = droplevels(as.factor(Surface)))

gam_2000_2005 <- cross_basis_bySurface(atp_00_05, k1=5, k2=5, "2000-2005")
p_diff_00 <- plot_surface_diffs(gam_2000_2005, "2000-2005", c(-0.4, 0.4))
print(p_diff_00)

# --- BLOC 2: 2006-2010 ---
atp_06_10 <- atp_tennis %>%
  filter(year(Date) >= 2006, year(Date) <= 2010, Rank_1 <= 200, Rank_2 <= 200) %>%
  mutate(millor_rank = pmin(Rank_1, Rank_2), Y = if_else(Winner == if_else(Rank_1 < Rank_2, Player_1, Player_2), 1, 0), ranking_diff = abs(Rank_1 - Rank_2)) %>%
  filter(Surface != "Carpet") %>% mutate(Surface = droplevels(as.factor(Surface)))

gam_2006_2010 <- cross_basis_bySurface(atp_06_10, k1=5, k2=5, "2006-2010")
p_diff_06 <- plot_surface_diffs(gam_2006_2010, "2006-2010", c(-0.4, 0.4))
print(p_diff_06)

# --- BLOC 3: 2011-2015 ---
atp_11_15 <- atp_tennis %>%
  filter(year(Date) >= 2011, year(Date) <= 2015, Rank_1 <= 200, Rank_2 <= 200) %>%
  mutate(millor_rank = pmin(Rank_1, Rank_2), Y = if_else(Winner == if_else(Rank_1 < Rank_2, Player_1, Player_2), 1, 0), ranking_diff = abs(Rank_1 - Rank_2)) %>%
  filter(Surface != "Carpet") %>% mutate(Surface = droplevels(as.factor(Surface)))

gam_2011_2015 <- cross_basis_bySurface(atp_11_15, k1=5, k2=5, "2011-2015")
p_diff_11 <- plot_surface_diffs(gam_2011_2015, "2011-2015", c(-0.4, 0.4))
print(p_diff_11)

# --- BLOC 4: 2016-2020 ---
atp_16_20 <- atp_tennis %>%
  filter(year(Date) >= 2016, year(Date) <= 2020, Rank_1 <= 200, Rank_2 <= 200) %>%
  mutate(millor_rank = pmin(Rank_1, Rank_2), Y = if_else(Winner == if_else(Rank_1 < Rank_2, Player_1, Player_2), 1, 0), ranking_diff = abs(Rank_1 - Rank_2)) %>%
  filter(Surface != "Carpet") %>% mutate(Surface = droplevels(as.factor(Surface)))

gam_2016_2020 <- cross_basis_bySurface(atp_16_20, k1=5, k2=5, "2016-2020")
p_diff_16 <- plot_surface_diffs(gam_2016_2020, "2016-2020", c(-0.4, 0.4))
print(p_diff_16)

# --- BLOC 5: 2021-2025 ---
atp_21_25 <- atp_tennis %>%
  filter(year(Date) >= 2021, year(Date) <= 2025, Rank_1 <= 200, Rank_2 <= 200) %>%
  mutate(millor_rank = pmin(Rank_1, Rank_2), Y = if_else(Winner == if_else(Rank_1 < Rank_2, Player_1, Player_2), 1, 0), ranking_diff = abs(Rank_1 - Rank_2)) %>%
  filter(Surface != "Carpet") %>% mutate(Surface = droplevels(as.factor(Surface)))

gam_2021_2025 <- cross_basis_bySurface(atp_21_25, k1=5, k2=5, "2021-2025")
p_diff_21 <- plot_surface_diffs(gam_2021_2025, "2021-2025", c(-0.4, 0.4))
print(p_diff_21)

# ==============================================================================
# 4. COMPARATIVA EVOLUTIVA (GRÀFICS INTERACTIUS)
# ==============================================================================

all_gams <- list(
  "2000-2005" = gam_2000_2005,
  "2006-2010" = gam_2006_2010,
  "2011-2015" = gam_2011_2015,
  "2016-2020" = gam_2016_2020,
  "2021-2025" = gam_2021_2025
)
valors_test <- c(1, 10, 20, 30, 40, 50)

for (val in valors_test) {
  dades <- get_curve_diff_fixed(all_gams, val)
  pl <- plot_comparison_interactive(dades, "millor_rank", "Diff_Rank =", val, "Millor rànquing")
  print(pl)
}
for (val in valors_test) {
  dades <- get_curve_rank_fixed(all_gams, val)
  pl <- plot_comparison_interactive(dades, "ranking_diff", "Millor_Rank =", val, "Diferència de rànquing")
  print(pl)
}