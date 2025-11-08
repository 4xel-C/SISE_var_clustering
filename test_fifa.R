# =============================================================================
# TEST FIFA COMPLET AVEC VISUALISATIONS
# =============================================================================

library(devtools)
load_all()

# Charger données
cat("\n=== CHARGEMENT FIFA ===\n")
fifa <- read.csv("data/players_22_cleaned.csv")
skill_cols <- grep("pace|shooting|passing|dribbling|defending|physic|skill_|attacking_|movement_|power_|mentality_|goalkeeping_", 
                   names(fifa), value = TRUE)
fifa_skills <- fifa[, skill_cols]
cat("✓", nrow(fifa), "joueurs,", ncol(fifa_skills), "variables\n")

# Fit K-means
cat("\n=== FIT K-MEANS (K=5) ===\n")
km <- KmeansVariables$new(n_clusters = 5, n_init = 10, random_state = 123)
km$fit(fifa_skills)
km$summary()

# Créer dossier plots
if (!dir.exists("plots")) dir.create("plots")

# Métriques
cat("\n=== CALCUL MÉTRIQUES ===\n")
sil <- kmeans_silhouette(fifa_skills, km$clusters, km$centroids)
cat("✓ Silhouette moyenne:", round(mean(sil$silhouette), 4), "\n")

ch <- kmeans_calinski_harabasz(fifa_skills, km$clusters, km$centroids)
cat("✓ CH Index:", round(ch, 2), "\n")

contrib <- kmeans_contributions(fifa_skills, km$clusters, km$centroids)

# Tableau de corrélations variable-cluster
cat("\n=== TABLEAU CORRÉLATIONS VARIABLE-CLUSTER ===\n")
cor_table <- kmeans_correlation_table(fifa_skills, km$clusters, km$centroids)
cat("✓ Tableau créé avec", nrow(cor_table), "variables\n")
cat("✓ Variables ambiguës (separation < 0.2):",
    sum(cor_table$separation < 0.2), "/", nrow(cor_table), "\n")

# VISUALISATIONS
cat("\n=== GÉNÉRATION PLOTS ===\n")

# Plot 1: Silhouette avec noms de variables (horizontal)
png("plots/fifa_silhouette.png", width = 1200, height = 1400)
plot_kmeans_silhouette(sil, main = "Silhouette par Variable (FIFA)", cex_names = 0.7)
dev.off()
cat("✓ Silhouette plot (avec noms) sauvegardé\n")

# Plot 2: Contributions par cluster
png("plots/fifa_contributions.png", width = 1400, height = 600)
plot_kmeans_contributions(contrib, top_n = 5)
dev.off()
cat("✓ Contributions sauvegardées\n")

# Plot 3: Projection PCA des variables
png("plots/fifa_projection_pca.png", width = 1000, height = 1000)
plot_kmeans_projection(fifa_skills, km$clusters,
                       main = "Projection des variables FIFA (PCA)",
                       show_labels = TRUE,
                       label_threshold = 0.4)
dev.off()
cat("✓ Projection PCA sauvegardée\n")

# Plot 4: Summary panel complet
png("plots/fifa_summary_panel.png", width = 1600, height = 1200)
plot_kmeans_summary(km, fifa_skills, top_n = 3)
dev.off()
cat("✓ Summary panel sauvegardé\n")

# Plot 5-6: Choix de K
cat("\n=== ANALYSE K OPTIMAL ===\n")
elbow_data <- kmeans_elbow(fifa_skills, k_range = 2:8, n_init = 5, random_state = 123)

png("plots/fifa_elbow.png", width = 800, height = 600)
plot_kmeans_elbow(elbow_data)
dev.off()
cat("✓ Elbow plot sauvegardé\n")

sil_range <- kmeans_silhouette_range(fifa_skills, k_range = 2:8, n_init = 5, random_state = 123)

png("plots/fifa_silhouette_range.png", width = 800, height = 600)
plot_kmeans_silhouette_range(sil_range)
dev.off()
cat("✓ Silhouette range sauvegardée\n")

# Plot 7: Diagnostics complet
png("plots/fifa_diagnostics.png", width = 1200, height = 1000)
plot_kmeans_diagnostics(fifa_skills, k_range = 2:8, fitted_k = 5, n_init = 5, random_state = 123)
dev.off()
cat("✓ Diagnostics panel sauvegardé\n")

# Afficher métriques
cat("\n=== MÉTRIQUES PAR K ===\n")
ch_range <- kmeans_calinski_harabasz_range(fifa_skills, k_range = 2:8, n_init = 5, random_state = 123)
metrics <- merge(elbow_data, sil_range, by = "k")
metrics <- merge(metrics, ch_range, by = "k")
print(metrics)

cat("\n✅ TERMINÉ ! 7 plots dans plots/\n")

# =============================================================================
# TABLEAU RÉCAPITULATIF POUR LE RAPPORT
# =============================================================================

cat("\n=== CRÉATION TABLEAU RÉCAPITULATIF ===\n")

# Créer dossier results
if (!dir.exists("results")) dir.create("results")

# Interpréter les clusters (basé sur les top variables)
cluster_names <- c(
  "1" = "Offensif",
  "2" = "Vitesse",
  "3" = "Défensif", 
  "4" = "Physique",
  "5" = "Technique"
)

# Extraire top 3 variables par cluster
top_vars_list <- lapply(1:5, function(k) {
  cluster_contrib <- contrib[contrib$cluster == k, ]
  top3 <- head(cluster_contrib[order(-cluster_contrib$contribution), ], 3)
  paste(top3$variable, collapse = ", ")
})

# Créer le tableau
cluster_summary <- data.frame(
  Cluster = 1:5,
  Interprétation = unname(cluster_names[as.character(1:5)]),
  N_variables = as.numeric(table(km$clusters)),
  Inertie_intra = round(km$cluster_inertias, 4),
  Corrélation_moyenne = sapply(1:5, function(k) {
    var_indices <- which(km$clusters == k)
    cors <- sapply(var_indices, function(j) {
      abs(cor(fifa_skills[, j], km$centroids[, k]))
    })
    round(mean(cors), 4)
  }),
  Top_3_variables = unlist(top_vars_list),
  stringsAsFactors = FALSE
)

# Sauvegarder CSV
write.csv(cluster_summary, "results/cluster_summary.csv", row.names = FALSE)
cat("✓ Tableau sauvegardé : results/cluster_summary.csv\n")

# Afficher dans la console
cat("\n=== RÉSUMÉ DES CLUSTERS ===\n")
print(cluster_summary)

# Créer version LaTeX pour le rapport
latex_table <- knitr::kable(cluster_summary, format = "latex",
                            caption = "Résumé des 5 clusters de variables FIFA",
                            booktabs = TRUE)
writeLines(latex_table, "results/cluster_table.tex")
cat("✓ Tableau LaTeX sauvegardé : results/cluster_table.tex\n")

# Sauvegarder le tableau de corrélations
write.csv(cor_table, "results/correlation_table.csv", row.names = FALSE)
cat("✓ Tableau corrélations sauvegardé : results/correlation_table.csv\n")

# Afficher les variables les plus ambiguës
cat("\n=== TOP 5 VARIABLES AMBIGUËS ===\n")
ambiguous <- cor_table[order(cor_table$separation), ]
print(head(ambiguous[, c("variable", "assigned_cluster", "max_cor", "separation")], 5))

# Afficher les variables les mieux assignées
cat("\n=== TOP 5 VARIABLES MIEUX ASSIGNÉES ===\n")
strong <- cor_table[order(-cor_table$separation), ]
print(head(strong[, c("variable", "assigned_cluster", "max_cor", "separation")], 5))


