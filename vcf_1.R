knitr::opts_chunk$set(echo = TRUE)

library(ggplot2)
library(tidyr)
library(dplyr)

vcf_1 <- read.table("C:/Users/margo/OneDrive - Uniwersytet im. Adama Mickiewicza w Poznaniu/gen pop/vcf_C_pyg_22_forward_indelRm_intersect_annotated_SNP_is_alt.bed")[,c(1,3,6,20)]

vcf_2 <- read.table("C:/Users/margo/OneDrive - Uniwersytet im. Adama Mickiewicza w Poznaniu/gen pop/vcf_C_ruf_08_forward_indelRm_intersect_annotated_SNP_is_alt.bed")[,c(1,3,6,20)]

colnames(vcf_1) <- c("chr", "pos", "chCADD", "FORMAT")
colnames(vcf_2) <- c("chr", "pos", "chCADD", "FORMAT")

vcf_1$gat <- "C_pyg"
vcf_2$gat <- "C_ruf"

vcf <- rbind(vcf_1, vcf_2)
# policz warianty
n_variants <- nrow(vcf_1)

# wyciągnij informacje o homozygotach/heterozygotach. Użyj funkcji separate z paczki tidyr. 
vcf_1 %>%
  separate(FORMAT, into = "genotype", sep = ":")

# jeśli chcesz pozbyć się ostrzeżenia, dodaj paramter extra:
vcf_1 %>%
  separate(FORMAT, into = "genotype", sep = ":", extra = "drop") 


# zmiana 0/0 na hom_ref etc. 
vcf_1 <- vcf_1 %>%
  separate(FORMAT, into = "genotype", sep = ":", extra = "drop") %>%
  mutate(genotype = ifelse(genotype == "0/1", "HET", 
                           ifelse(genotype == "1/1", "HOM_ALT", 
                                  ifelse(genotype == "0/0", "HOM_REF", "error"))))


# policz heterozygoty i homozygoty:
table(vcf_1$genotype)


# policz ponownie:
vcf_1 %>%
  group_by(genotype) %>%
  tally()

n_HET <- vcf_1 %>%
  group_by(genotype) %>%
  tally() %>%
  filter(genotype == "HET") %>%
  pull()

n_HOM <- vcf_1 %>%
  group_by(genotype) %>%
  tally() %>%
  filter(genotype == "HOM") %>%
  pull()

paste0("W danych zidentyfikowano ", n_HET, " heterozygot i ", n_HOM, " homozygot.")
```


```{r wizualizacja}
ggplot(vcf_1, aes(pos, chCADD)) +
  geom_point() +
  ggtitle("chCADD, scaffold X, individual Y", 
          subtitle = paste0(n_variants, " variants"))

# zwizualizuj chCADD wzdłuż scaffoldu z podziałem na genotyp.
```


```{r pyg_26 testowanie}
head(vcf_1)
ggplot(vcf_1, aes(chCADD)) +
  geom_histogram() +
  facet_wrap(~genotype)


vcf_1 %>%
  mutate(chCADD_log = log10(chCADD)) %>%
  ggplot(aes(chCADD_log)) +
  geom_histogram() +
  facet_wrap(~genotype)

vcf_1 %>%
  mutate(chCADD_log = log10(chCADD)) %>%
  ggplot(aes(chCADD_log)) +
  geom_density() +
  facet_wrap(~genotype)
```
###


```{r chCADD ~ genotyp}
# Wykonamy to testem- najpierw warto sprawdzić rozkład w grupach:
vcf_1 %>%
  ggplot(aes(chCADD)) +
  geom_histogram() + 
  facet_grid(row = "genotype")

# brak rozkładu normalnego --> transformacja danych
# logartym (dodaj 1 żeby uniknąć log(0)) --> źle
vcf_1 %>%
  ggplot(aes(log(chCADD) + 1)) +
  geom_histogram() + 
  facet_grid(row = "genotype")
# pierwiastek --> wygląda OK
vcf_1 %>%
  ggplot(aes(sqrt(chCADD))) +
  geom_histogram() + 
  facet_grid(row = "genotype")
# sprawdź po rozkład po pierwiastkowaniu testem Shapiro-Wilka (sprawdza rozkład normalny)
x <- vcf_1 %>%
  filter(genotype == "HET") %>%
  mutate(chCADD_sqrt = sqrt(chCADD))

shapiro.test(x$chCADD_sqrt) 


Swilcox.test(chCADD ~ genotype, data = vcf_1)

vcf %>%
  ggplot(aes(chCADD)) +
  geom_density() +
  facet_wrap(~gat)


vcf %>%
  ggplot(aes(x = gat, y = chCADD)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = .2)


wilcox.test(chCADD ~ gat, data = vcf)


```{r}
# utworzenie losowego zestawu danych
df <- data.frame(wartość = c(rgamma(100, shape = 2, scale = 2),
                             rgamma(100, shape = 2, scale = 3.9)),
                 grupa = c(rep("A", 500), 
                           rep("B", 500)))

print(df[c(1:5, 500:505),])

```

Obliczamy średnią wartość dla grupy A i B i odejmujemy A - B.

```{r}
obs_diff <- mean(df$wartość[df$grupa == "A"]) - mean(df$wartość[df$grupa == "B"])

paste0("Różnica średniej wartości wynosi ", round(obs_diff, 3))
```

```{r}
# 1000 powtórzeń
perm_diffs <- vector("numeric", 1000)

for(i in 1:1000) {
  # Losowo zmieniaj przyporządkowanie wartości do grupy:
  perm_groups <- sample(df$grupa)
  
  # Oblicz różnice średniej wartości pomiędzy grupami A i B dla każdego powtórzenia:
  perm_diffs[i] <- mean(df$wartość[perm_groups == "A"]) - mean(df$wartość[perm_groups == "B"])
}

p_value <- sum(perm_diffs <= obs_diff) / length(perm_diffs)

print(p_value)

# wizualizacja:
perm_diffs_df <- data.frame(perm_diffs = perm_diffs)

ggplot(perm_diffs_df, aes(perm_diffs)) +
  geom_histogram() +
  geom_vline(xintercept = obs_diff, linetype = "dashed")


```