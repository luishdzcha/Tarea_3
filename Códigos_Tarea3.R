###Códigos de toda la tarea 3

#Cargamos nuestralibrería 
library(phyloseq)
library(microbiome)
library(ggplot2)
library(vegan)
library(dplyr)


#Para cargar el objeto phyloseq
data("dietswap", package = "microbiome")
ps <- dietswap
ps

#Curvas de Rarefacción 

# Extraer matriz de OTUs y transponer
otu_mat <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)  # Asegurar que taxones sean filas

# Generar curvas de Rarefacción
rarecurve(otu_mat, step = 50, label = FALSE, col = "blue", ylab = "Especies", xlab ="Número de lecturas")

# Guardamos la gráfica en pdf
pdf("Figuras_Tarea3/Rarefacción_dietswap.pdf", width = 8, height = 6)  
rarecurve(otu_mat, step = 50, label = FALSE, col = "blue", 
          main = "Curvas de Rarefacción (dietswap)",
          xlab = "Número de lecturas", ylab = "Especies")
dev.off()  # para cerra bien el pdf

#Diversidad alfa

# Calcular y graficar índices alfa (Este es de los grupos)
plot_richness(ps, x = "group", measures = c("Observed", "Shannon", "Simpson")) + geom_boxplot(aes(fill = group)) + theme_bw()       


#Guardarlo en la carpeta
pdf("Figuras_Tarea3/Índicesalfa_dietswap.pdf", width = 8, height = 5)  
plot_richness(ps, x = "group", measures = c("Observed", "Shannon", "Simpson")) + geom_boxplot(aes(fill = group)) + theme_bw()  
dev.off()  

#Filtrado y transformación


ps_relativa <- transform_sample_counts(ps, function(x) x / sum(x) * 100)  # Abundancia relativa %

taxa_abundantes <- genefilter_sample(ps_relativa, filterfun_sample(function(x) x >= 0.1), A = 0.1 * nsamples(ps_relativa))  
ps_filtrado <- prune_taxa(taxa_abundantes, ps) 
#Ahora tenemos los datos filtrados de las taxas con .1% de abundancia relativa en el 10% de la muestra

#Diversidad Beta

# Calcular distancia Bray-Curtis
dist_bray <- phyloseq::distance(ps_filtrado, method = "bray")

# Realizar PCoA (Ordination)
pcoa <- ordinate(ps_filtrado, method = "PCoA", distance = dist_bray)

# Graficar PCoA
plot_ordination(ps_filtrado, ordinate(ps_filtrado, "PCoA", distance = "bray"), color = "group", title = "PCoA (Bray-Curtis)") + geom_point(size = 3) + stat_ellipse(level = 0.95) + theme_minimal()

#Lo guardamos en la carpeta
pdf("Figuras_Tarea3/DiversidadBeta_dietswap.pdf", width = 8, height = 5)
plot_ordination(ps_filtrado, ordinate(ps_filtrado, "PCoA", distance = "bray"), color = "group", title = "PCoA (Bray-Curtis)") + geom_point(size = 3) + stat_ellipse(level = 0.95) + theme_minimal()
dev.off()

#Gráficas Rank Abundance
otu_df <- as.data.frame(t(otu_table(ps_relativa)))  # Transponer para taxones en columnas
colnames(otu_df) <- tax_table(ps_relativa)[, "Genus"]  # Asignar nombres de géneros

#Calcular abundancias promedio por taxón
taxa_avg <- sort(colMeans(otu_df), decreasing = TRUE)
rank_df <- data.frame(
  Rank = seq_along(taxa_avg),
  Taxon = names(taxa_avg),
  Abundance = taxa_avg
)

head(rank_df, 10) # Estos son los géneros más abundantes
#Gráfica
ggplot(rank_df, aes(x = Rank, y = Abundance)) + geom_point(size = 3) + geom_line() + scale_y_log10() + labs(title = "Curva de Abundancia-Rango", x = "Rango", y = "Abundancia relativa (log10)") + theme_minimal()

#Guardamos la gráfica
pdf("Figuras_Tarea3/RankAbundance_dietswap.pdf", width = 8, height = 5)
ggplot(rank_df, aes(x = Rank, y = Abundance)) + geom_point(size = 3) + geom_line() + scale_y_log10() + labs(title = "Curva de Abundancia-Rango", x = "Rango", y = "Abundancia relativa (log10)") + theme_minimal()
dev.off()

#GRáfica apiladas de abundancia por taxón

#Gráfica por Phylum
plot_bar(ps_filtrado, fill = "Phylum") +
  geom_bar(stat = "identity") +
  labs(title = "Composición por Phylum", 
       x = "Muestras", 
       y = "Abundancia relativa") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
#Guardamos la gráfica
pdf("Figuras_Tarea3/GráficaPhylum_dietswap.pdf", width = 8, height = 5)
plot_bar(ps_filtrado, fill = "Phylum") +
  geom_bar(stat = "identity") +
  labs(title = "Composición por Phylum", 
       x = "Muestras", 
       y = "Abundancia relativa") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


#Gráfica por Género 
plot_bar(ps_filtrado, fill = "Genus") +
  geom_bar(stat = "identity") +
  labs(title = "Composición por Género", 
       x = "Muestras", 
       y = "Abundancia relativa") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
#Guardamos la gráfica
pdf("Figuras_Tarea3/GráficaGénero_dietswap.pdf", width = 8, height = 5)
plot_bar(ps_filtrado, fill = "Genus") +
  geom_bar(stat = "identity") +
  labs(title = "Composición por Género", 
       x = "Muestras", 
       y = "Abundancia relativa") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
# Exportar resultados

####GlobalPatterns

data("GlobalPatterns")
gp <- GlobalPatterns
#Preprocesamiento 

# 1. Filtrar taxones con <5 lecturas en ≥20% de muestras
gp_filtrado <- filter_taxa(gp, function(x) sum(x >= 5) >= 0.2 * nsamples(gp), prune = TRUE)

# 2. Aglomerar a nivel de Familia
gp_familia <- tax_glom(gp_filtrado, taxrank = "Family")

# 3. Transformar a abundancias relativas (%)
gp_relativo <- transform_sample_counts(gp_familia, function(x) x / sum(x) * 100)

# 4. Subset para Soil, Feces, Skin
gp_final <- subset_samples(gp_relativo, SampleType %in% c("Soil", "Feces", "Skin"))

gp_final

#Diversidad alfa
# 1. Calcular índices (usando gp_familia, no gp_final)
alpha_div <- estimate_richness(gp_familia, measures = c("Shannon", "Simpson", "Observed"))

# 2. Unir con metadatos y graficar
alpha_meta <- cbind(sample_data(gp_final), alpha_div)  # Usar gp_final para los metadatos

# Verificar muestras en alpha_div vs gp_final
rownames(alpha_div) %in% rownames(sample_data(gp_final))

# Filtrar alpha_div para incluir solo muestras presentes en gp_final
alpha_div <- alpha_div[rownames(alpha_div) %in% rownames(sample_data(gp_final)), ]
#Alpha_meta corregido para que coincida
alpha_meta <- cbind(
  sample_data(gp_final)[rownames(alpha_div), ],  # Metadatos en el mismo orden
  alpha_div
)


# 3. Boxplots comparativos por tipo de muestra (Soil, Feces, Skin)
plot_shannon <- ggplot(alpha_meta, aes(x = SampleType, y = Shannon, fill = SampleType)) +
  geom_boxplot() +
  labs(title = "Índice de Shannon por Tipo de Muestra", x = "", y = "Shannon") +
  theme_minimal()
plot_shannon
plot_simpson <- ggplot(alpha_meta, aes(x = SampleType, y = Simpson, fill = SampleType)) +
  geom_boxplot() +
  labs(title = "Índice de Simpson por Tipo de Muestra", x = "", y = "Simpson") +
  theme_minimal()
plot_simpson
plot_observed <- ggplot(alpha_meta, aes(x = SampleType, y = Observed, fill = SampleType)) +
  geom_boxplot() +
  labs(title = "Riqueza Observada por Tipo de Muestra", x = "", y = "Observed") +
  theme_minimal()
plot_observed
#Guardamos los boxplots de cada uno
#El de Shannon
pdf("Figuras_Tarea3/Gráficashannon_GP.pdf", width = 8, height = 5)
plot_shannon <- ggplot(alpha_meta, aes(x = SampleType, y = Shannon, fill = SampleType)) +
  geom_boxplot() +
  labs(title = "Índice de Shannon por Tipo de Muestra", x = "", y = "Shannon") +
  theme_minimal()
plot_shannon
dev.off()
#EL de Simpson
pdf("Figuras_Tarea3/GráficaSimpson_GP.pdf", width = 8, height = 5)
plot_simpson <- ggplot(alpha_meta, aes(x = SampleType, y = Simpson, fill = SampleType)) +
  geom_boxplot() +
  labs(title = "Índice de Simpson por Tipo de Muestra", x = "", y = "Simpson") +
  theme_minimal()
plot_simpson
dev.off()
#El de Observed
pdf("Figuras_Tarea3/GráficaObserved_GP.pdf", width = 8, height = 5)
plot_observed <- ggplot(alpha_meta, aes(x = SampleType, y = Observed, fill = SampleType)) +
  geom_boxplot() +
  labs(title = "Riqueza Observada por Tipo de Muestra", x = "", y = "Observed") +
  theme_minimal()
plot_observed
dev.off()

# 4. Prueba de Kruskal-Wallis para cada índice
kruskal_shannon <- kruskal.test(Shannon ~ SampleType, data = alpha_meta)
kruskal_simpson <- kruskal.test(Simpson ~ SampleType, data = alpha_meta)
kruskal_observed <- kruskal.test(Observed ~ SampleType, data = alpha_meta)
# Mostrar resultados estadísticos
list(
  Shannon = kruskal_shannon,
  Simpson = kruskal_simpson,
  Observed = kruskal_observed
)

#Curvas de Rango-Abundancia

# 1. Preparar datos (abundancia relativa por muestra)
ps_rel <- transform_sample_counts(gp_final, function(x) x / sum(x))  # Abundancia relativa (0-1)

# 2. Extraer matriz de abundancia y metadatos
otu_df <- as.data.frame(otu_table(ps_rel))
meta_df <- as.data.frame(sample_data(ps_rel))

# 3. Función para generar datos de rango-abundancia por muestra
prepare_rank_data <- function(sample_id) {
  sample_abund <- sort(otu_df[, sample_id], decreasing = TRUE)
  sample_abund <- sample_abund[sample_abund > 0]  # Filtrar ceros
  data.frame(
    SampleID = sample_id,
    Rank = seq_along(sample_abund),
    Abundance = sample_abund,
    SampleType = meta_df[sample_id, "SampleType"]
  )
}

# 4. Aplicar a todas las muestras y combinar
rank_data <- bind_rows(lapply(colnames(otu_df), prepare_rank_data))

# 5. Gráfico (log10 en Y, colorear por SampleType)
ggplot(rank_data, aes(x = Rank, y = Abundance, color = SampleType)) +
  geom_line(aes(group = SampleID), alpha = 0.6) +  # Líneas por muestra
  scale_y_log10() +  # Escala logarítmica en Y
  labs(
    title = "Curvas de Rango-Abundancia por Tipo de Muestra",
    x = "Rango (de más a menos abundante)",
    y = "Abundancia relativa (log10)",
    color = "Ambiente"
  ) +
  theme_minimal()
#Guardamos la gráfica
pdf("Figuras_Tarea3/GráficaRA_GP.pdf", width = 8, height = 5)
ggplot(rank_data, aes(x = Rank, y = Abundance, color = SampleType)) +
  geom_line(aes(group = SampleID), alpha = 0.6) +  # Líneas por muestra
  scale_y_log10() +  # Escala logarítmica en Y
  labs(
    title = "Curvas de Rango-Abundancia por Tipo de Muestra",
    x = "Rango (de más a menos abundante)",
    y = "Abundancia relativa (log10)",
    color = "Ambiente"
  ) +
  theme_minimal()
dev.off()
#Perfil taxonómico 

# 1. Agrupar a nivel de Phylum y calcular abundancia relativa
gp_phylum <- gp_final %>%
  tax_glom(taxrank = "Phylum") %>%               # Aglomerar por Phylum
  transform_sample_counts(function(x) x/sum(x))  # Abundancia relativa (0-1)

# 2. Identificar los 5 phyla más abundantes (promedio global)
top_phyla <- names(sort(taxa_sums(gp_phylum), decreasing = TRUE)[1:5])
gp_top5 <- prune_taxa(top_phyla, gp_phylum)      # Filtrar solo los 5 principales

# 3. Preparar datos para ggplot
df_plot <- psmelt(gp_top5) %>%                   # Convertir a dataframe
  group_by(SampleType, Phylum) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop")  # Promedio por ambiente

# 4. Gráfico de barras apiladas + facet_wrap
ggplot(df_plot, aes(x = SampleType, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill") +          # Barras al 100%
  scale_y_continuous(labels = scales::percent) +            # Eje Y en porcentaje
  facet_wrap(~ SampleType, scales = "free_x") +             # Separar por ambiente
  labs(
    title = "Perfil Taxonómico: Top 5 Phyla por Ambiente",
    x = "",
    y = "Abundancia relativa",
    fill = "Phylum"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_blank())                      # Ocultar etiquetas redundantes
#Guardamos la gráica
pdf("Figuras_Tarea3/Gráficaperfiltax_GP.pdf", width = 8, height = 5)
ggplot(df_plot, aes(x = SampleType, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill") +          # Barras al 100%
  scale_y_continuous(labels = scales::percent) +            # Eje Y en porcentaje
  facet_wrap(~ SampleType, scales = "free_x") +             # Separar por ambiente
  labs(
    title = "Perfil Taxonómico: Top 5 Phyla por Ambiente",
    x = "",
    y = "Abundancia relativa",
    fill = "Phylum"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_blank())
dev.off()
#Diversidad Beta

# 1. Calcular distancia Bray-Curtis
dist_bray <- phyloseq::distance(gp_final, method = "bray")

# 2. Realizar PCoA
pcoa <- ordinate(gp_final, method = "PCoA", distance = dist_bray)

# 3. Visualización (colores por SampleType + elipses)
plot_pcoa <- plot_ordination(gp_final, pcoa, color = "SampleType") +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.95, linetype = 2) +  # Elipses del 95%
  labs(
    title = "PCoA basado en Bray-Curtis",
    caption = paste("Stress:", round(pcoa$values$Relative_eig[1] * 100, 1), "% explicado por PCo1")
  ) +
  theme_minimal()
#Guardamos la gráfica
pdf("Figuras_Tarea3/GráficaBray-Curtis_GP.pdf", width = 8, height = 5)
plot_pcoa <- plot_ordination(gp_final, pcoa, color = "SampleType") +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.95, linetype = 2) +  # Elipses del 95%
  labs(
    title = "PCoA basado en Bray-Curtis",
    caption = paste("Stress:", round(pcoa$values$Relative_eig[1] * 100, 1), "% explicado por PCo1")
  ) +
  theme_minimal()
plot_pcoa
dev.off()

# 4. Stress plot (para evaluar calidad del PCoA)
stressplot(cmdscale(dist_bray, k = 2), dist_bray)  # Me sale que vegan espera un objeto de un análisis NMDS no PCoA
#Y que en ese caso se utilizaría scree plot, ya que stressplot solo se con objetos de tipo metaMDS
# 5. PERMANOVA para diferencias entre grupos
permanova <- adonis2(dist_bray ~ SampleType, data = as(sample_data(gp_final), "data.frame"))
print(permanova)

# Mostrar gráfico PCoA
plot_pcoa
