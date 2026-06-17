# ============================================================
# Monte Carlo metodų taikymas statistinių hipotezių tikrinimui
# genų raiškos duomenyse
# Margiris Antanas Malakauskas
# ============================================================


# ------------------------------------------------------------
# 1. Paketų įkėlimas
# -> GEOquery leidžia parsisiųsti duomenis tiesiai iš NCBI GEO duomenų bazės
# -> BiocManager yra Bioconductor paketų diegimo įrankis
# ------------------------------------------------------------

library(GEOquery)


# ------------------------------------------------------------
# 2. Duomenų atsiuntimas ir patikrinimas
# -> Parsisiunčiamas GSE1297 iš GEO duomenų bazės
# -> Patikrinamos dimensijos: eilutės = genai, stulpeliai = samplai
# -> Peržiūrima mėginių informaciją (pavadinimas, grupės)
# -> Suskaičiuojamas samplus pagal grupę
# ------------------------------------------------------------

gse <- getGEO("GSE1297", GSEMatrix = TRUE)
gse <- gse[[1]]

dim(exprs(gse))
pData(gse)[, 1:3]
table(pData(gse)$title)


# ------------------------------------------------------------
# 3. Kontrolinių ir sunkių AD samplų išskyrimas
# -> grepl() paieška grupių pavadinimuose
# -> Mėginių patikrinimas - tikiuosi 9 kontroliniai ir 7 sunkūs AD 
# -> Išskiriama ekspresijos matricą tik šioms dviem grupėms
# ------------------------------------------------------------

pheno <- pData(gse)

control_idx <- grepl("Control", pheno$title)
severe_idx  <- grepl("Severe",  pheno$title)

sum(control_idx)
sum(severe_idx)

expr         <- exprs(gse)
expr_control <- expr[, control_idx]
expr_severe  <- expr[, severe_idx]

dim(expr_control)  # 22283 x 9
dim(expr_severe)   # 22283 x 7


# ------------------------------------------------------------
# 4. Neapdorotų duomenų peržiūra
# -> Tikrinama reikšmių skalė - didelės reikšmės rodo log transformacijos būtinumą
# -> Vidurkis >> mediana rodo "right skew"
# ------------------------------------------------------------

expr_control[1:5, ]
expr_severe[1:5, ]
summary(expr_control[, 1])


# ------------------------------------------------------------
# 5. Log2 transformacija
# -> Transformacija stabilizuoja dispersiją ir sumažina asimetriją
# -> Po transformacijos vidurkis ≈ mediana - skirstinys simetriškesnis
# -> Vizualiai palyginame pasiskirstymą prieš ir po transformacijos
# ------------------------------------------------------------

expr_control_log <- log2(expr_control + 1)
expr_severe_log  <- log2(expr_severe  + 1)

summary(expr_control_log[, 1])

par(mfrow = c(1, 2))

hist(expr_control[, 1],
     main   = "Prieš log transformaciją",
     xlab   = "Ekspresijos lygis",
     col    = "lightblue",
     breaks = 50)

hist(expr_control_log[, 1],
     main   = "Po log2 transformacijos",
     xlab   = "log2 Ekspresijos lygis",
     col    = "lightgreen",
     breaks = 50)

par(mfrow = c(1, 1))


# ------------------------------------------------------------
# 6. Vieno geno analizė (Genas 1: 1007_s_at)
# -> Analizuojamas vienas genas, patikrinti metodo logiką
# -> Rankiniu būdu suskaičiuota t-statistiką ir palyginta su t.test()
# -> Monte Carlo: maišo grupes 10000 kartų ir skaičiuoja empirinę p-reikšmę
# -> Vizualizuojamas permutacijų pasiskirstymas su stebima statistika
# ------------------------------------------------------------

gene1_control <- expr_control_log[1, ]
gene1_severe  <- expr_severe_log[1, ]

mean(gene1_control)
mean(gene1_severe)
mean(gene1_severe) - mean(gene1_control)

n1        <- length(gene1_control)
n2        <- length(gene1_severe)
mean_diff <- mean(gene1_severe) - mean(gene1_control)
var1      <- var(gene1_control)
var2      <- var(gene1_severe)
se        <- sqrt(var1/n1 + var2/n2)
t_stat    <- mean_diff / se

cat("Vidurkių skirtumas:", mean_diff, "\n")
cat("Standartinė paklaida:", se, "\n")
cat("t-statistika:", t_stat, "\n")

t_result <- t.test(gene1_severe, gene1_control)
t_result

gene1_all <- c(gene1_control, gene1_severe)
obs_diff  <- mean(gene1_severe) - mean(gene1_control)
B         <- 10000

set.seed(42)
perm_diffs <- numeric(B)

for (i in 1:B) {
  shuffled      <- sample(gene1_all)
  perm_severe   <- shuffled[1:7]
  perm_control  <- shuffled[8:16]
  perm_diffs[i] <- mean(perm_severe) - mean(perm_control)
}

p_monte_carlo <- (1 + sum(abs(perm_diffs) >= abs(obs_diff))) / (B + 1)

cat("Stebimas skirtumas:", obs_diff, "\n")
cat("Klasikinis t-testo p:", t_result$p.value, "\n")
cat("Monte Carlo p:", p_monte_carlo, "\n")

hist(perm_diffs,
     breaks = 50,
     main   = "Permutacijų pasiskirstymas: Genas 1007_s_at",
     xlab   = "Vidurkių skirtumas (Sunkus AD - Kontrolė)",
     ylab   = "Dažnis",
     col    = "lightblue",
     border = "white")
abline(v =  obs_diff, col = "red", lwd = 2, lty = 1)
abline(v = -obs_diff, col = "red", lwd = 2, lty = 2)
legend("topright",
       legend = c("Stebimas skirtumas", "Veidrodinis (dvipusis)"),
       col    = "red",
       lty    = c(1, 2),
       lwd    = 2)


# ------------------------------------------------------------
# 7. Visų 22 283 genų analizė
# -> Kiekvienam genui atliekamas t-testas ir Monte Carlo permutacinis testas
# -> B = 5000 permutacijų - kompromisas tarp tikslumo ir skaičiavimo laiko
# -> Rezultatai saugomi duomenų lentelėje
# ------------------------------------------------------------

n_genes     <- nrow(expr_control_log)
ttest_pvals <- numeric(n_genes)
mc_pvals    <- numeric(n_genes)
mean_diffs  <- numeric(n_genes)

B <- 5000
set.seed(42)

for (g in 1:n_genes) {
  ctrl <- expr_control_log[g, ]
  sev  <- expr_severe_log[g, ]
  
  tt             <- t.test(sev, ctrl)
  ttest_pvals[g] <- tt$p.value
  mean_diffs[g]  <- tt$estimate[1] - tt$estimate[2]
  
  all_vals <- c(ctrl, sev)
  obs      <- mean(sev) - mean(ctrl)
  perm_d   <- numeric(B)
  
  for (i in 1:B) {
    s         <- sample(all_vals)
    perm_d[i] <- mean(s[1:7]) - mean(s[8:16])
  }
  
  mc_pvals[g] <- (1 + sum(abs(perm_d) >= abs(obs))) / (B + 1)
}

gene_names <- rownames(expr_control_log)

results <- data.frame(
  gene      = gene_names,
  mean_diff = mean_diffs,
  p_ttest   = ttest_pvals,
  p_mc      = mc_pvals
)

cat("Rezultatai:", nrow(results), "genų\n")
head(results)


# ------------------------------------------------------------
# 8. Daugybinio tikrinimo korekcija (Benjamini-Hochberg FDR)
# -> Vienu metu tikrinamos 22 283 hipotezes
# -> Benjamini-Hochberg metodas kontroliuoja klaidingo atradimo dažnį
# -> Palyginamas reikšmingų genų skaičius prieš ir po korekcijos
# ------------------------------------------------------------

results$p_ttest_adj <- p.adjust(results$p_ttest, method = "BH")
results$p_mc_adj    <- p.adjust(results$p_mc,    method = "BH")

cat("| Nekoreguota p-reikšmė < 0.05 |\n")
cat("t-testas - reikšmingi genai:       ", sum(results$p_ttest < 0.05), "\n")
cat("Monte Carlo - reikšmingi genai:    ", sum(results$p_mc    < 0.05), "\n\n")

cat("| FDR koreguota p-reikšmė < 0.05 |\n")
cat("t-testas - reikšmingi genai:       ", sum(results$p_ttest_adj < 0.05), "\n")
cat("Monte Carlo - reikšmingi genai:    ", sum(results$p_mc_adj    < 0.05), "\n")


# ------------------------------------------------------------
# 9. Metodų palyginimas
# -> Spearman koreliacija matuoja metodų sutarimą genų reitingavime
# -> Sutapimo tarp metodų "Dot diagram"
# -> Identifikuojamas genas, kur metodai labiausiai skiriasi
# ------------------------------------------------------------

correlation <- cor(results$p_ttest, results$p_mc, method = "spearman")
cat("Spearman koreliacija tarp t-testo ir MC p-reikšmių:", round(correlation, 4), "\n")

plot(results$p_ttest, results$p_mc,
     main = "t-testo ir Monte Carlo p-reikšmių palyginimas",
     xlab = "t-testo p-reikšmė",
     ylab = "Monte Carlo p-reikšmė",
     pch  = 16,
     cex  = 0.3,
     col  = rgb(0, 0, 1, 0.2))
abline(a = 0, b = 1, col = "red", lwd = 2)

top_gene <- results[results$p_ttest_adj < 0.05, ]
cat("\nGeriausias genas pagal t-testą (FDR koreguota):\n")
print(top_gene)

results_sorted_mc <- results[order(results$p_mc), ]
cat("\nTop 10 genų pagal Monte Carlo p-reikšmę:\n")
print(head(results_sorted_mc[, c("gene", "mean_diff", "p_ttest", "p_mc")], 10))


# ------------------------------------------------------------
# 10. Specifinio (tarp meodų labiausiai nesutarto) geno analizė: Genas 203894_at
# -> t-testas: p = 6.18e-07 (reikšmingas), Monte Carlo: p = 0.0002 (nereikšmingas)
# -> Permutacijų diagrama parodokad joks permutavimas 
#    neviršijo stebimo skirtumo, tad MC negali tiksliau įvertinti p-reikšmės
# ------------------------------------------------------------

gene_idx  <- which(results$gene == "203894_at")
ctrl_star <- expr_control_log[gene_idx, ]
sev_star  <- expr_severe_log[gene_idx, ]

cat("Kontrolės reikšmės:\n"); print(round(ctrl_star, 3))
cat("Sunkaus AD reikšmės:\n"); print(round(sev_star, 3))
cat("Vidurkių skirtumas:", round(mean(sev_star) - mean(ctrl_star), 4), "\n")

obs_star <- mean(sev_star) - mean(ctrl_star)
all_star  <- c(ctrl_star, sev_star)

set.seed(42)
perm_star <- numeric(5000)
for (i in 1:5000) {
  s            <- sample(all_star)
  perm_star[i] <- mean(s[1:7]) - mean(s[8:16])
}

hist(perm_star,
     breaks = 50,
     main   = "Permutacijų pasiskirstymas: Genas 203894_at",
     xlab   = "Vidurkių skirtumas (Sunkus AD - Kontrolė)",
     ylab   = "Dažnis",
     col    = "lightblue",
     border = "white")
abline(v =  obs_star, col = "red", lwd = 2, lty = 1)
abline(v = -obs_star, col = "red", lwd = 2, lty = 2)
legend("topright",
       legend = c("Stebimas skirtumas", "Veidrodinis (dvipusis)"),
       col    = "red",
       lty    = c(1, 2),
       lwd    = 2)


# ------------------------------------------------------------
# 11. Veikimo laiko palyginimas
# -> t-testas naudoja paprastą matematinę formulę - labai greitas
# -> Monte Carlo kartoja 5000 atsitiktinių maišymų kiekvienam genui, dideliame
#    genų rinkynyję
# -> Matuojamas MC laikas 500 genų imtyje ir ekstrapoliuojamas visiems
# -> Skirtumas rodo vieną iš Monte Carl trūkumū - skaičiavimo laiką
# ------------------------------------------------------------

time_ttest <- system.time({
  ttest_pvals_timed <- numeric(n_genes)
  for (g in 1:n_genes) {
    ctrl <- expr_control_log[g, ]
    sev  <- expr_severe_log[g, ]
    ttest_pvals_timed[g] <- t.test(sev, ctrl)$p.value
  }
})

time_mc_500 <- system.time({
  mc_pvals_timed <- numeric(500)
  set.seed(42)
  for (g in 1:500) {
    ctrl     <- expr_control_log[g, ]
    sev      <- expr_severe_log[g, ]
    all_vals <- c(ctrl, sev)
    obs      <- mean(sev) - mean(ctrl)
    perm_d   <- numeric(5000)
    for (i in 1:5000) {
      s         <- sample(all_vals)
      perm_d[i] <- mean(s[1:7]) - mean(s[8:16])
    }
    mc_pvals_timed[g] <- (1 + sum(abs(perm_d) >= abs(obs))) / (5001)
  }
})

cat("| Skaičiavimo laiko palyginimas |\n\n")
cat("t-testas (visi 22 283 genai):\n")
cat("  Laikas:", round(time_ttest["elapsed"], 2), "sekundžių\n\n")
cat("Monte Carlo B=5000 (500 genų):\n")
cat("  Laikas:", round(time_mc_500["elapsed"], 2), "sekundžių\n")
cat("  Įvertintas laikas visiems 22 283 genams:",
    round(time_mc_500["elapsed"] / 500 * 22283, 1), "sekundžių\n")
cat("  Tai yra apytiksliai",
    round(time_mc_500["elapsed"] / 500 * 22283 / 60, 1), "minučių\n\n")
cat("Monte Carlo yra apytiksliai",
    round((time_mc_500["elapsed"] / 500 * 22283) / time_ttest["elapsed"]),
    "kartų lėtesnis už t-testą\n")


# ------------------------------------------------------------
# 12. Permutacijų skaičiaus stabilumo analizė
# -> Tikrinama kaip p-reikšmė keičiasi didėjant permutacijų skaičiui B
# -> Kiekvienas B pakartotas 20 kartų - matuojamas vidurkis ir standartinįs nuokrypis
# -> Genas 1007_s_at (ribinis): didelė variacija ties B=100, stabilizuojasi ~5000
# -> Genas 203894_at (stiprus): pasiskirsto ties grindų reikšme 1/(B+1)
# ------------------------------------------------------------

B_values <- c(100, 500, 1000, 2000, 5000, 10000)
n_reps   <- 20

g1_all <- c(expr_control_log[1, ], expr_severe_log[1, ])
g1_obs <- mean(expr_severe_log[1, ]) - mean(expr_control_log[1, ])

g2_idx <- which(rownames(expr_control_log) == "203894_at")
g2_all <- c(expr_control_log[g2_idx, ], expr_severe_log[g2_idx, ])
g2_obs <- mean(expr_severe_log[g2_idx, ]) - mean(expr_control_log[g2_idx, ])

g1_pvals <- matrix(NA, nrow = length(B_values), ncol = n_reps)
g2_pvals <- matrix(NA, nrow = length(B_values), ncol = n_reps)

set.seed(42)

for (b_idx in 1:length(B_values)) {
  B_curr <- B_values[b_idx]
  
  for (rep in 1:n_reps) {
    
    perm_d <- numeric(B_curr)
    for (i in 1:B_curr) {
      s         <- sample(g1_all)
      perm_d[i] <- mean(s[1:7]) - mean(s[8:16])
    }
    g1_pvals[b_idx, rep] <- (1 + sum(abs(perm_d) >= abs(g1_obs))) / (B_curr + 1)
    
    perm_d <- numeric(B_curr)
    for (i in 1:B_curr) {
      s         <- sample(g2_all)
      perm_d[i] <- mean(s[1:7]) - mean(s[8:16])
    }
    g2_pvals[b_idx, rep] <- (1 + sum(abs(perm_d) >= abs(g2_obs))) / (B_curr + 1)
  }
  
  cat("Atlikta B =", B_curr, "\n")
}

g1_mean <- apply(g1_pvals, 1, mean)
g1_sd   <- apply(g1_pvals, 1, sd)
g2_mean <- apply(g2_pvals, 1, mean)
g2_sd   <- apply(g2_pvals, 1, sd)

cat("| Genas 1007_s_at (ribinis reikšmingas) |\n")
cat(sprintf("%-10s %-14s %-14s\n", "B", "Vidurkio p", "Stand. nuokr."))
for (i in 1:length(B_values)) {
  cat(sprintf("%-10d %-14.4f %-14.4f\n", B_values[i], g1_mean[i], g1_sd[i]))
}

cat("\n| Genas 203894_at (stipriai reikšmingas) |\n")
cat(sprintf("%-10s %-14s %-14s\n", "B", "Vidurkio p", "Stand. nuokr."))
for (i in 1:length(B_values)) {
  cat(sprintf("%-10d %-14.4f %-14.4f\n", B_values[i], g2_mean[i], g2_sd[i]))
}

par(mfrow = c(1, 2))

plot(B_values, g1_mean,
     type  = "b",
     pch   = 16,
     col   = "steelblue",
     ylim  = c(0, max(g1_mean + g1_sd) * 1.2),
     main  = "Genas 1007_s_at\n(ribinis reikšmingas)",
     xlab  = "Permutacijų skaičius (B)",
     ylab  = "Monte Carlo p-reikšmė",
     log   = "x")
arrows(B_values, g1_mean - g1_sd,
       B_values, g1_mean + g1_sd,
       angle = 90, code = 3, length = 0.05, col = "steelblue")
abline(h = results$p_ttest[1], col = "red", lty = 2, lwd = 2)
legend("topright", legend = "t-testo p-reikšmė", col = "red", lty = 2, lwd = 2)

plot(B_values, g2_mean,
     type  = "b",
     pch   = 16,
     col   = "darkgreen",
     ylim  = c(0, max(g2_mean + g2_sd) * 1.2),
     main  = "Genas 203894_at\n(stipriai reikšmingas)",
     xlab  = "Permutacijų skaičius (B)",
     ylab  = "Monte Carlo p-reikšmė",
     log   = "x")
arrows(B_values, g2_mean - g2_sd,
       B_values, g2_mean + g2_sd,
       angle = 90, code = 3, length = 0.05, col = "darkgreen")
abline(h = results$p_ttest[g2_idx], col = "red", lty = 2, lwd = 2)
legend("topright", legend = "t-testo p-reikšmė", col = "red", lty = 2, lwd = 2)

par(mfrow = c(1, 1))