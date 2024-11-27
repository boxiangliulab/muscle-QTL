# Figure S14A
ggplot(expr_ANK1, aes(x = Time, y = Expression, group = Sample)) +
  geom_line(color = "grey") +
  geom_point(aes(color = Time), size = 3, position = position_dodge(0.2)) +
  stat_summary(fun = mean, geom = "line", aes(group = 1), color = "blue", size = 1) +
  labs(title = "ANKT Expression (TPM)",
       subtitle = "Paired t-test p-value: 0.0106",
       x = NULL,
       y = "ANKT Exp (TPM)") +
  scale_color_manual(values = c("Pre" = "red", "Post" = "blue")) +
  theme_minimal() +
  theme(legend.position = "none")

  # Figure S14B
ggplot(isoform_comparison, aes(x = Variant, y = Value, color = Group)) +
  geom_point(position = position_dodge(0.2), size = 3) +
  geom_errorbar(aes(ymin = Value - se, ymax = Value + se), width = 0.2, position = position_dodge(0.2)) +
  scale_color_manual(values = c("Chinese obese" = "blue", "Chinese lean" = "red")) +
  labs(x = NULL, y = "-dCt Value") +
  theme_minimal() +
  theme(legend.position = "none")

  # Figure S14C
  ggplot(expression_by_ANK1, aes(x = Genotype, y = Expression, fill = Time)) +
  geom_boxplot(outlier.color = "black", position = position_dodge(0.8)) +
  scale_fill_manual(values = c("Pre" = "red", "Post" = "blue")) +
  labs(x = "rs508419", y = "ANK1 Expression") +
  theme_minimal() +
  theme(legend.position = "right", legend.title = element_blank())
