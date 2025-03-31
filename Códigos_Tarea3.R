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
pdf("Figuras_Tarea3/DiversidadBeta_dietswap", width = 8, height = 5)
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
pdf("Figuras_Tarea3/RankAbundance_dietswap", width = 8, height = 5)
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

#Gráfica por Género 
plot_bar(ps_filtrado, fill = "Genus") +
  geom_bar(stat = "identity") +
  labs(title = "Composición por Género", 
       x = "Muestras", 
       y = "Abundancia relativa") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

####GlobalPatterns

data("GlobalPatterns")
gp <- GlobalPatterns
#Preprocesamiento 


#Diversidad alfa

#Curvas de Rango-Abundancia

#Perfil taxonómico 


#Diversidad Beta


