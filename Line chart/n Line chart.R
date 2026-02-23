



library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# 
normal_data <- data.frame(
  p_N = c(0.1, 0.3, 0.6, 0.9, 0.1, 0.3, 0.6, 0.9, 0.1, 0.3),
  N = c(50, 50, 50, 50, 100, 100, 100, 100, 200, 200),
  pT_CsF_c1 = c(0.028, 0.066, 0.102, 0.311, 0.039, 0.052, 0.12, 0.438, 0.046, 0.07),
  pT_CsF_c2 = c(0.05, 0.067, 0.049, 0.133, 0.052, 0.051, 0.074, 0.27, 0.043, 0.08),
  pT_CsF_cr = c(0.05, 0.067, 0.041, 0.057, 0.052, 0.048, 0.057, 0.072, 0.042, 0.079),
  pEBA4_rls = c(0.057, 0.048, 0.055, 0.104, 0.039, 0.028, 0.072, 0.141, 0.036, 0.025),
  pols_2_rls = c(0.057, 0.051, 0.036, 0.022, 0.038, 0.026, 0.007, 0.003, 0.036, 0.007)
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
    title = "Rejection Rates by p/N Ratio (Normal Distribution)",
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



