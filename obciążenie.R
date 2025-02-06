
library(ggplot2)

h <- 0.1

L_hom <- c(18, 127)  
L_het <- c(45, 75)
osobniki <- c("C_pyg_22", "C_ruf_08")


genetic_load <- L_hom * 1 + L_het * 0.5
realized_load <- L_hom * 1 + L_het * h

df <- data.frame(
  Osobnik = rep(osobniki, 2),
  Typ = rep(c("Całkowite", "Zrealizowane"), each = 2),
  Obciążenie = c(genetic_load, realized_load)
)

ggplot(df, aes(x = Osobnik, y = Obciążenie, fill = Typ)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Obciążenie genetyczne",
       x = "Osobnik", y = "Obciążenie") +
  theme_minimal() +
  scale_fill_manual(values = c("blue", "green"))
