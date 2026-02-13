library(tidyverse)
library(lme4)
library(lmerTest)     
library(performance)  
library(broom.mixed)  
library(ggplot2)



library(tidyverse)
library(sf)
library(ggplot2)
library(scales)
library(patchwork)

library(tidyverse)
library(sf)
library(tigris)
library(ggplot2)

# read data
setwd("~/Desktop/CU/BIS9185/BISTP9185-Proj2")
df <- read_csv("tract_cbsa.csv", show_col_types = FALSE)

# Quick checks for NA 
glimpse(df)
colSums(is.na(df)) %>% sort(decreasing = TRUE) %>% head(20)

df <- df %>%
  filter(!is.na(pm25), !is.na(bc))

# Basic cleaning / transforms
df <- df %>%
  mutate(
    cbsa_id    = as.character(cbsa_id),
    ruca_agg   = factor(ruca_agg, levels = c("urban", "suburban", "rural"))
  )




######################
#     plots!!!!!!    #
######################
# state level summary of the data
states_sf <- tigris::states(cb = TRUE, class = "sf") %>%
  filter(!STUSPS %in% c("AK","HI","PR")) %>%
  st_transform(5070)

state_summary <- df %>%
  group_by(state_fips, state_abb) %>%
  summarise(
    pm25_mean = mean(pm25, na.rm = TRUE),
    dens_mean = mean(dens, na.rm = TRUE),
    pct_white_mean = mean(pct_white, na.rm = TRUE),
    pct_black_mean = mean(pct_black, na.rm = TRUE),
    pct_hispanic_mean = mean(pct_hispanic, na.rm = TRUE),
    pct_asian_mean = mean(pct_asian, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(state_fips = str_pad(as.character(state_fips), 2, pad = "0"))

state_map_df <- states_sf %>%
  mutate(state_fips = STATEFP, state_abb = STUSPS) %>%
  left_join(state_summary, by = c("state_fips","state_abb"))

state_map_df2 <- state_map_df %>%
  mutate(
    white_p = pct_white_mean / 100,
    black_p = pct_black_mean / 100,
    asian_p = pct_asian_mean / 100,
    hisp_p  = pct_hispanic_mean / 100
  ) %>% filter(!is.na(pm25_mean))



# Theme helper: make maps fill more area + consistent titles
theme_map <- theme_void() +
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text  = element_text(size = 10),
    legend.position = "right"
  )

# ------------------------------------------------------------
# (A) PM2.5 + Density as ONE faceted plot (aligned sizing)

top_long <- state_map_df2 %>%
  st_as_sf() %>%
  select(geometry, pm25_mean, dens_mean) %>%
  pivot_longer(
    cols = c(pm25_mean, dens_mean),
    names_to = "measure",
    values_to = "value"
  ) %>%
  mutate(
    measure = recode(
      measure,
      pm25_mean = "Mean PM2.5 (2010)",
      dens_mean = "Population Density"
    )
  )

# PM2.5 plot (own scale)
p_pm25 <- ggplot(state_map_df2) +
  geom_sf(aes(fill = pm25_mean), color = "white", linewidth = 0.2) +
  coord_sf(datum = NA, expand = FALSE) +
  scale_fill_gradient(
    low = "#E3F2FD", high = "#0D47A1",
    na.value = "grey90",
    name = "PM2.5"
  ) +
  theme_map +
  labs(title = "Mean PM2.5 (2010)") +
  theme(
    plot.title = element_text(
      hjust = 0.5,      # center
      face = "bold",    # bold
      size = 12        # optional: bigger text
    )
  )

# Density plot (own scale)
p_dens <- ggplot(state_map_df2) +
  geom_sf(aes(fill = dens_mean), color = "white", linewidth = 0.2) +
  coord_sf(datum = NA, expand = FALSE) +
  scale_fill_gradient(
    low = "#E3F2FD", high = "#0D47A1",
    na.value = "grey90",
    name = "Density"
  ) +
  labs(title = "Population Density") +
  theme_map +
  theme(
    plot.title = element_text(
      hjust = 0.5,      # center
      face = "bold",    # bold
      size = 12        # optional: bigger text
    )
  )

# Combine side-by-side
p_top <- p_pm25 | p_dens
p_top


# Race as ONE faceted plot (0–1 scale, aligned sizing)

race_long <- state_map_df2 %>%
  st_as_sf() %>%
  select(geometry, white_p, black_p, asian_p, hisp_p) %>%
  pivot_longer(
    cols = c(white_p, black_p, asian_p, hisp_p),
    names_to = "race",
    values_to = "value"
  ) %>%
  mutate(
    race = recode(
      race,
      white_p = "Proportion White",
      black_p = "Proportion Black",
      asian_p = "Proportion Asian",
      hisp_p  = "Proportion Hispanic"
    )
  )

p_race <- ggplot(race_long) +
  geom_sf(aes(fill = value), color = "white", linewidth = 0.2) +
  coord_sf(datum = NA) +
  facet_wrap(~ race, ncol = 2) +
  scale_fill_gradient(
    low = "#E3F2FD",   # light blue
    high = "#0D47A1",  # dark blue
    limits = c(0, 1),
    oob = squish,
    labels = percent_format(accuracy = 1),
    na.value = "grey90"
  )+
  theme_map

# ------------------------------------------------------------
# Combine: make bottom row taller so race maps aren't tiny
# ------------------------------------------------------------
combined_plot <- p_top / p_race +
  plot_layout(heights = c(1, 2.3)) 

combined_plot


ggsave(
  "us_maps_large.png",
  plot = combined_plot,
  width = 18,   # increase size
  height = 14,  # increase size
  dpi = 300
)



# EDA

# PM2.5 histogram
ggplot(df, aes(x = pm25)) +
  geom_histogram(bins = 40) +
  labs(
    title = "Distribution of PM2.5 across census tracts",
    x = "PM2.5 (2010 modeled, tract-level)",
    y = "Number of tracts"
  ) +
  theme_minimal()

# PM2.5 by RUCA category
ggplot(df, aes(x = ruca_agg, y = pm25)) +
  geom_boxplot(outlier.alpha = 0.2) +
  labs(
    title = "PM2.5 by rural-urban category (RUCA_agg)",
    x = "RUCA category",
    y = "PM2.5"
  ) +
  theme_minimal()

# by race percentile
race_vars <- c(
  pct_white    = "White",
  pct_black    = "Black",
  pct_asian    = "Asian",
  pct_hispanic = "Hispanic"
)

df %>%
  # keep needed columns and drop missing pm25
  select(pm25, all_of(names(race_vars))) %>%
  filter(!is.na(pm25)) %>%
  # long format: one row per (tract, race_group)
  pivot_longer(
    cols = all_of(names(race_vars)),
    names_to = "race_var",
    values_to = "pct"
  ) %>%
  mutate(
    race_group = recode(race_var, !!!race_vars),
    
    # decile bins, include 100% reliably by going to 101 and using right=FALSE
    pct_bin = cut(
      pct,
      breaks = seq(0, 101, by = 10),
      right = FALSE,
      include.lowest = TRUE
    ),
    
    # nicer labels like "0–10", "10–20", ..., "90–100"
    pct_bin_label = case_when(
      is.na(pct_bin) ~ NA_character_,
      TRUE ~ str_replace_all(as.character(pct_bin), "\\[|\\)|\\]", "") %>%
        str_replace(",", "–") %>%
        str_replace("100–101", "90–100")   # just in case it prints oddly
    )
  ) %>%
  filter(!is.na(pct_bin)) %>%
  group_by(race_group, pct_bin, pct_bin_label) %>%
  summarise(
    mean_pm25 = mean(pm25, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  # order x-axis bins correctly
  mutate(
    pct_bin_label = factor(
      pct_bin_label,
      levels = c("0–10","10–20","20–30","30–40","40–50","50–60","60–70","70–80","80–90","90–100")
    )
  ) %>%
  ggplot(aes(x = pct_bin_label, y = mean_pm25, color = race_group, group = race_group)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    title = "Mean PM2.5 by Race/Ethnicity Share Deciles",
    x = "Percent of tract population (binned into deciles)",
    y = "Mean PM2.5",
    color = "Group"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# State-level summaries (useful even without tract geometries)
state_pm <- df %>%
  group_by(state_abb) %>%
  summarise(
    n_tracts = n(),
    pm25_mean = mean(pm25, na.rm = TRUE),
    pm25_p95  = quantile(pm25, 0.95, na.rm = TRUE),
    pct_white_mean = mean(pct_white, na.rm = TRUE),
    dens_mean = mean(dens, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(pm25_mean))

print(head(state_pm, 20))

# Top-20 states by mean PM2.5 (bar)
state_pm %>%
  slice_head(n = 20) %>%
  ggplot(aes(x = reorder(state_abb, pm25_mean), y = pm25_mean)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Top 20 states by mean tract-level PM2.5 (2010)",
    x = "State",
    y = "Mean PM2.5"
  ) +
  theme_minimal()

# Correlations (selected vars)
corr_vars <- c(
  "pm25","bc","hyads","dens",
  "pct_white","pct_black",
  "pct_hispanic","pct_asian","pct_pov"
)

# readable labels
var_labels <- c(
  pm25         = "PM2.5",
  bc           = "Black Carbon",
  hyads        = "Coal Plant Impact",
  dens         = "Population Density",
  pct_white    = "% White",
  pct_black    = "% Black",
  pct_hispanic = "% Hispanic",
  pct_asian    = "% Asian",
  pct_pov      = "% Poverty"
)

# correlation matrix
cor_mat <- df %>%
  select(all_of(corr_vars)) %>%
  cor(use = "pairwise.complete.obs")

# reshape + relabel
cor_long <- cor_mat %>%
  as.data.frame() %>%
  rownames_to_column("var1") %>%
  pivot_longer(-var1, names_to = "var2", values_to = "correlation") %>%
  mutate(
    var1 = recode(var1, !!!var_labels),
    var2 = recode(var2, !!!var_labels)
  )

# heatmap
ggplot(cor_long, aes(var1, var2, fill = correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", correlation)),
            size = 3) +
  
  scale_fill_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Correlation"
  ) +
  
  coord_equal() +
  
  labs(
    title = "Correlation Heatmap: Environmental & Demographic Variables",
    x = NULL,
    y = NULL
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )




















# -------------------------------
# 2) Hierarchical model: lmer
# -------------------------------
# Primary mixed model:
#   - outcome: pm25
#   - exposure: pct_white_10 (per 10 percentage points)
#   - confounders: log_dens, pct_pov, ruca_agg, bc, hyads
#   - hierarchy: (1|state_fips) + (1|state_fips:cbsa_id)
#     (CBSA random intercept nested within state)

df_model <- df %>%
  mutate(
    # nested cluster ID
    cbsa_in_state = interaction(state_fips, cbsa_id, drop = TRUE),
  )

# -------------------------
# 1) Primary multilevel model
#    Random intercepts: state + CBSA-within-state
# -------------------------
# m1 <- lmer(
#   pm25 ~ pct_white + dens + pct_pov + ruca_agg + bc + hyads +
#     (1 | state_fips) + (1 | cbsa_in_state),
#   data = df_model,
#   REML = TRUE
# )
# 
# summary(m1)
# confint(m1, method = "Wald")  # quick CI
# 
# # Tidy fixed effects
# broom.mixed::tidy(m1, effects = "fixed", conf.int = TRUE)
# 
# # Diagnostics
# performance::check_model(m1)
# 
# # Interpretation helper: effect per +10 percentage points White
# beta_white <- fixef(m1)["pct_white_10"]
# cat("\nEstimated change in PM2.5 for +10 percentage points White:",
#     round(beta_white, 4), "ug/m^3\n")

# ------------------------------------------------
# 2) Sensitivity A: random slope for % White by CBSA-in-state
#    (more flexible, can be slower / may warn if near-singular)
# ------------------------------------------------
m2 <- lmer(
  pm25 ~ pct_white + log(dens) + pct_pov + ruca_agg + bc + hyads +
    (1 | state_fips) + (1 + pct_white | cbsa_in_state),
  data = df_model,
  REML = TRUE
)

summary(m2)
performance::check_model(m2)

# Compare m1 vs m2 (REML=TRUE isn’t ideal for LRT; refit with REML=FALSE if you want formal LRT)
m1_ml <- update(m1, REML = FALSE)
m2_ml <- update(m2, REML = FALSE)
anova(m1_ml, m2_ml)



