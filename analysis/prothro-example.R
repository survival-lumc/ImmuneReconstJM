?prothro

tdCox.pro <- coxph(Surv(start, stop, event) ~ pro + treat,
                   data = prothro)
summary(tdCox.pro)

dt <- as.data.table(prothro)
last <- dt[, .SD[.N], by = "id"]
prothro |>
  ggplot(aes(time, pro, group = id)) +
  #geom_point() +
  geom_line(alpha = 0.75) +
  geom_point(data = last, size = 2, col = "red") +
  facet_wrap(~ death)

last |>
  ggplot(aes(x = factor(death), y = pro, fill = factor(death))) +
  geom_violin() +
  #geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  geom_point(position = position_dodge2(width = 0.1), alpha = 0.5)

last[, start_new := 0]
coxph(Surv(stop, event) ~ pro + treat, data = last) |>  coef()
