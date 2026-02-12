# load packages
library(tidyverse)
library(lme4)
library(lmerTest)     # p-values for lmer (Satterthwaite)
library(performance)  # model diagnostics helpers
library(broom.mixed)  # tidy() for lmer
library(ggplot2)

# -------------
# 0) Load data
# -------------
setwd("~/Desktop/CU/BIS9185/BISTP9185-Proj2")
df <- read_csv("tract_cbsa.csv", show_col_types = FALSE)

# Quick checks
glimpse(df)
colSums(is.na(df)) %>% sort(decreasing = TRUE) %>% head(20)

# Basic cleaning / transforms
df <- df %>%
  mutate(
    state_fips = str_pad(as.character(state_fips), 2, pad = "0"),
    cbsa_id    = as.character(cbsa_id),
    ruca_agg   = factor(ruca_agg, levels = c("urban", "suburban", "rural")),
    pct_white_10 = pct_white / 10,                   # per 10 percentage points
    log_dens   = log(dens + 1e-6)                    # avoid log(0)
  )

# Drop rows missing key outcomes/exposures
df_ana <- df %>%
  filter(!is.na(pm25), !is.na(bc))

# -----------------------------
# 1) EDA: distributions & bivar
# -----------------------------

# PM2.5 histogram
ggplot(df_ana, aes(x = pm25)) +
  geom_histogram(bins = 40) +
  labs(
    title = "Distribution of PM2.5 across census tracts",
    x = "PM2.5 (2010 modeled, tract-level)",
    y = "Number of tracts"
  ) +
  theme_minimal()

# PM2.5 by RUCA category
ggplot(df_ana, aes(x = ruca_agg, y = pm25)) +
  geom_boxplot(outlier.alpha = 0.2) +
  labs(
    title = "PM2.5 by rural-urban category (RUCA_agg)",
    x = "RUCA category",
    y = "PM2.5"
  ) +
  theme_minimal()

# Mean PM2.5 by % White deciles
df_ana %>%
  mutate(pct_white_bin = cut(pct_white, breaks = seq(0, 100, by = 10), right = FALSE)) %>%
  group_by(pct_white_bin) %>%
  summarise(mean_pm25 = mean(pm25, na.rm = TRUE), n = n(), .groups = "drop") %>%
  ggplot(aes(x = pct_white_bin, y = mean_pm25, group = 1)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Mean PM2.5 by % White deciles",
    x = "% White (binned)",
    y = "Mean PM2.5"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# State-level summaries (useful even without tract geometries)
state_pm <- df_ana %>%
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
corr_vars <- c("pm25","bc","hyads","dens","pct_white","pct_black","pct_hispanic","pct_asian","pct_pov")
cor_mat <- df_ana %>%
  select(all_of(corr_vars)) %>%
  cor(use = "pairwise.complete.obs")
print(round(cor_mat, 3))

# -------------------------------
# 2) Hierarchical model: lmer
# -------------------------------
# Primary mixed model:
#   - outcome: pm25
#   - exposure: pct_white_10 (per 10 percentage points)
#   - confounders: log_dens, pct_pov, ruca_agg, bc, hyads
#   - hierarchy: (1|state_fips) + (1|state_fips:cbsa_id)
#     (CBSA random intercept nested within state)

df_model <- df_ana %>%
  filter(
    !is.na(pm25),
    !is.na(pct_white_10),
    !is.na(log_dens),
    !is.na(pct_pov),
    !is.na(ruca_agg),
    !is.na(bc),
    !is.na(hyads),
    !is.na(state_fips),
    !is.na(cbsa_id)
  ) %>%
  mutate(cbsa_in_state = interaction(state_fips, cbsa_id, drop = TRUE))

m1 <- lmer(
  pm25 ~ pct_white_10 + log_dens + pct_pov + ruca_agg + bc + hyads +
    (1 | state_fips) + (1 | cbsa_in_state),
  data = df_model,
  REML = TRUE
)

summary(m1)
confint(m1, method = "Wald")  # fast CI

# Tidy coefficients
broom.mixed::tidy(m1, effects = "fixed", conf.int = TRUE)

# Model diagnostics (residual checks)
performance::check_model(m1)

# Interpretation helper:
# effect per +10 percentage points white = beta for pct_white_10
beta_white <- fixef(m1)["pct_white_10"]
cat("\nEstimated change in PM2.5 for +10 percentage points White:", beta_white, "ug/m^3\n")

# ---------------------------------------------
# 3) Sensitivity checks you may want to include
# ---------------------------------------------

# (A) Random slopes for % White by CBSA (if you expect heterogeneity across metro areas)
# NOTE: may be heavier / slower.
m2 <- lmer(
  pm25 ~ pct_white_10 + log_dens + pct_pov + ruca_agg + bc + hyads +
    (1 | state_fips) + (1 + pct_white_10 | cbsa_in_state),
  data = df_model,
  REML = TRUE
)
summary(m2)

# (B) CBSA fixed effects (within-CBSA comparisons), plus state fixed effects
# This is NOT a random-effects model; it controls for all CBSA-level confounding.
m_fe <- lm(
  pm25 ~ pct_white_10 + log_dens + pct_pov + ruca_agg + bc + hyads +
    factor(state_fips) + factor(cbsa_id),
  data = df_model
)
summary(m_fe)

# If you want cluster-robust SE by CBSA for the FE model:
# install.packages(c("sandwich","lmtest"))
library(sandwich)
library(lmtest)
coeftest(m_fe, vcov = vcovCL(m_fe, cluster = ~ cbsa_id))

# -----------------------------------------
# 4) (Optional) Map / choropleth scaffolding
# -----------------------------------------
# You do NOT have tract geometries in tract_cbsa.csv.
# If your assignment allows downloading TIGER/Line tract shapefiles,
# you can do a tract choropleth by joining geometries on GEOID.
#
# Minimal example (requires internet + sf + tigris):
# install.packages(c("sf","tigris"))
# library(sf); library(tigris)
# options(tigris_use_cache = TRUE)
# tracts <- tigris::tracts(state = "NY", year = 2010, cb = TRUE)  # example state
# tracts <- tracts %>% mutate(geoid = GEOID)
# map_df <- tracts %>% left_join(df, by = "geoid")
# ggplot(map_df) + geom_sf(aes(fill = pm25), color = NA) + theme_void()