



library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# 
normal_data <- data.frame(
  p_N = c(0.1, 0.3, 0.6, 0.9, 0.1, 0.3, 0.6, 0.9, 0.1, 0.3),
  N = c(50, 50, 50, 50, 100, 100, 100, 100, 200, 200),
  pT_CsF_c1 = c(0.023, 0.061, 0.121, 0.31, 0.032, 0.071, 0.084, 0.283, 0.05, 0.056),
  pT_CsF_c2 = c(0.041, 0.062, 0.067, 0.115, 0.054, 0.071, 0.049, 0.138, 0.047, 0.068),
  pT_CsF_cr = c(0.041, 0.062, 0.056, 0.025, 0.054, 0.07, 0.035, 0.017, 0.047, 0.066),
  pEBA4_rls = c(0.047, 0.040, 0.085, 0.154, 0.043, 0.035, 0.077, 0.184, 0.028, 0.039),
  pols_2_rls = c(0.046, 0.041, 0.051, 0.044, 0.043, 0.027, 0.008, 0.008, 0.030, 0.008)
)

# 
test_names <- c(
  "pT_CsF_c1" = "CsF-c1",
  "pT_CsF_c2" = "CsF-c2",
  "pT_CsF_cr" = "CsF-cr",
  "pEBA4_rls" = "EBA4_rls",
  "pols_2_rls" = "ols_2_rls"
)

# 
data_long <- normal_data %>%
  pivot_longer(
    cols = -c(p_N, N),
    names_to = "Test",
    values_to = "Rejection_Rate"
  ) %>%
  mutate(
    Test = factor(Test, levels = names(test_names), labels = test_names),
    N = factor(N, levels = c(50, 100, 200)),
    Test_Group = case_when(
      grepl("F", Test) ~ "F-variants",
      grepl("CsF", Test) ~ "CsF-variants",
      grepl("EBA", Test) ~ "EBA-variants",
      TRUE ~ "Other"
    )
  )

# 
ggplot(data_long, aes(x = p_N, y = Rejection_Rate, color = Test, group = Test)) +
  geom_line(linewidth = 0.8, alpha = 0.9) +
  geom_point(size = 1.5, alpha = 0.9) +
  facet_wrap(~ N, nrow = 1, labeller = labeller(N = ~ paste("N =", .x))) +
  scale_y_continuous(
    limits = c(0, 0.5),                     # 
    breaks = seq(0, 0.5, 0.1),              # 
    expand = c(0, 0.02)
  ) +
  scale_x_continuous(
    breaks = unique(normal_data$p_N),
    labels = scales::number_format(accuracy = 0.1),
    expand = c(0.05, 0.05)
  ) +
  geom_hline(
    yintercept = 0.05, 
    linetype = "dashed", 
    color = "red", 
    linewidth = 0.5,
    alpha = 0.6
  ) +
  labs(
    title = "Rejection Rates by p/N Ratio (non-normality (skewness=1, kurtosis=3) Distribution)",
    x = "p/N Ratio",
    y = "Rejection Rate",
    color = "Test Method"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",        # 
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),      # 
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    axis.title = element_text(face = "bold", size = 10),
    axis.text = element_text(color = "black", size = 8),
    panel.spacing = unit(0.8, "lines")
  ) +
  scale_color_manual(values = scales::hue_pal()(length(test_names))) +
  guides(color = guide_legend(
    nrow = 1,                                # 
    override.aes = list(linewidth = 1.5, alpha = 1)
  ))

# 
ggsave("all_tests_rejection_rates.png", width = 14, height = 8, dpi = 300)



