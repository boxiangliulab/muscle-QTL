# Figure S5A
ggplot(A_bw, aes(x = Time, y = BodyWeight)) +
  geom_jitter(width = 0.2, alpha = 0.6, color = "steelblue") +
  stat_summary(fun = mean, geom = "point", color = "black", size = 3) +
  stat_summary(fun = mean, geom = "errorbar", width = 0.2, color = "black") +
  labs(title = "Body Weight Changes Over 12 Months", y = "Body Weight (kg)", x = "") +
  theme_minimal()

# Add pvalue manually

# Figure S5B
ggplot(B_bmi, aes(x = Time, y = BMI)) +
  geom_jitter(width = 0.2, alpha = 0.6, color = "steelblue") +
  stat_summary(fun = mean, geom = "point", color = "black", size = 3) +
  stat_summary(fun = mean, geom = "errorbar", width = 0.2, color = "black") +
  labs(title = "BMI Changes Over 12 Months", y = "BMI (kg/m²)", x = "") +
  theme_minimal()

# Figure S5C
ggplot(C_bw, aes(x = Time, y = BodyWeight, group = ID)) +
  geom_line(alpha = 0.5, color = "gray") +
  stat_summary(aes(group = 1), fun = mean, geom = "line", color = "red", size = 1) +
  stat_summary(aes(group = 1), fun = median, geom = "line", color = "blue", size = 1) +
  labs(title = "Body Weight Changes Across Key Time Points", y = "Body Weight (kg)", x = "") +
  theme_minimal()

# Figure S5D
ggplot(D_bmi, aes(x = Time, y = BMI, group = ID)) +
  geom_line(alpha = 0.5, color = "gray") +
  stat_summary(aes(group = 1), fun = mean, geom = "line", color = "red", size = 1) +
  stat_summary(aes(group = 1), fun = median, geom = "line", color = "blue", size = 1) +
  labs(title = "BMI Changes Across Key Time Points", y = "BMI (kg/m²)", x = "") +
  theme_minimal()
