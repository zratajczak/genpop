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
n_variants <- nrow(vcf_2)

# wyciągnij informacje o homozygotach/heterozygotach. Użyj funkcji separate z paczki tidyr. 
vcf_2 %>%
  separate(FORMAT, into = "genotyp", sep = ":")

# jeśli chcesz pozbyć się ostrzeżenia, dodaj paramter extra:
vcf_2 %>%
  separate(FORMAT, into = "genotyp", sep = ":", extra = "drop") 


# zmiana 0/0 na hom_ref etc. 
vcf_2 <- vcf_2 %>%
  separate(FORMAT, into = "genotyp", sep = ":", extra = "drop") %>%
  mutate(genotyp = ifelse(genotyp == "0/1", "HET", 
                           ifelse(genotyp == "1/1", "HOM_ALT", 
                                  ifelse(genotyp == "0/0", "HOM_REF", "error"))))


# policz heterozygoty i homozygoty:
table(vcf_2$genotyp)


# policz ponownie:
vcf_2 %>%
  group_by(genotyp) %>%
  tally()

n_HET <- vcf_2 %>%
  group_by(genotyp) %>%
  tally() %>%
  filter(genotyp == "HET") %>%
  pull()

n_HOM <- vcf_2 %>%
  group_by(genotyp) %>%
  tally() %>%
  filter(genotyp == "HOM") %>%
  pull()

paste0("W danych zidentyfikowano ", n_HET, " heterozygot i ", n_HOM, " homozygot.")
```


```{r wizualizacja}
ggplot(vcf_2, aes(pos, chCADD)) +
  geom_point() +
  ggtitle("chCADD, scaffold X, individual Y", 
          subtitle = paste0(n_variants, " variants"))

# zwizualizuj chCADD wzdłuż scaffoldu z podziałem na genotyp.
```


```{r pyg_26 testowanie}
head(vcf_2)
ggplot(vcf_2, aes(chCADD)) +
  geom_histogram() +
  facet_wrap(~genotyp)


vcf_2 %>%
  mutate(chCADD_log = log10(chCADD)) %>%
  ggplot(aes(chCADD_log)) +
  geom_histogram() +
  facet_wrap(~genotyp)

vcf_2 %>%
  mutate(chCADD_log = log10(chCADD)) %>%
  ggplot(aes(chCADD_log)) +
  geom_density() +
  facet_wrap(~genotyp)


```{r pyg_26 testowanie}
head(vcf_2)
ggplot(vcf_2, aes(chCADD)) +
  geom_histogram() +
  facet_wrap(~genotyp)

library(dplyr)

vcf_2 %>%
  mutate(chCADD_log = log10(chCADD)) %>%
  ggplot(aes(chCADD_log)) +
  geom_histogram() +
  facet_wrap(~genotyp)

```

vcf_2 %>%
  ggplot(aes(chCADD)) +
  geom_density() +
  facet_wrap(~genotyp)

vcf_2 %>%
  ggplot(aes(x = genotyp, y = chCADD)) +
  geom_boxplot() +
  geom_jitter(width = .2)


