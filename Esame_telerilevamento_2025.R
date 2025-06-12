# PROGETTO DI ANALISI NDVI E CLASSIFICAZIONE DEL SUOLO – NEW DELHI (2000–2023)

# Nicolò Garbo
# Esame telerilevamento geo ecologico 2025

# Obiettivo:
# Analizzare l’evoluzione della copertura vegetale e dell’urbanizzazione nella città 
# di New Delhi e nelle aree periurbane circostanti, attraverso immagini satellitari 
# multispettrali Landsat. L’analisi si concentra su tre anni: 2000, 2013, 2023.

# Contesto:
# New Delhi ha vissuto una rapida espansione urbana: in circa 25 anni la popolazione 
# è quasi triplicata. Questo ha generato un'enorme pressione sul territorio, con perdita 
# di vegetazione, consumo di suolo e crescita del cemento. L’NDVI (Normalized Difference 
# Vegetation Index) è utilizzato per quantificare la vegetazione sana, mentre la classificazione 
# spettrale permette di suddividere il territorio in aree tematiche (verde, suolo, cemento).

# Dati:
# - Immagini Landsat 7 (ETM+) e Landsat 8 (OLI)
# - Immagini degli anni 2000, 2013, 2023
# - Path 140, Row 051
# - 4 bande usate: B1 (Blue), B2 (Green), B3 (Red), B4 (NIR)
# - Cloud cover < 5%, mese di novembre (stagione secca stabile)
# - Download da USGS Earth Explorer




# Caricamento dei pacchetti
library(terra)         #  per gestire raster SpatRaster, crop(), ext(), rast()
library(ggplot2)       #  per grafici e mappe ggplot
library(viridis)       #  palette colori viridis (adatta anche per daltonici)
library(patchwork)     #  unisce più grafici ggplot in layout multipli
library(imageRy)       #  funzioni personalizzate per NDVI, RGB, esportazioni

# Directory di lavoro

setwd("C:/Users/Utente/Desktop/Esame_nuova_delhi_2025")

# Definizione del bounding box per l'area di interesse (New Delhi e fascia periurbana)
# EPSG:32643 – UTM Zone 43N (coordinate in metri)
# ext() proviene dal pacchetto terra

delhi_ext = ext(745000, 805000, 3125000, 3180000)  

### Importazione delle immagini satellitari dei tre anni 

# 2000 

files2000 = sort(list.files("2000", pattern = "SR_B[1-4]\\.TIF$", full.names = TRUE))
lsat2000 = rast(files2000)       # crea uno SpatRaster
lsat2000 = crop(lsat2000, delhi_ext)

# 2013 

files2013 = sort(list.files("2013", pattern = "SR_B[1-4]\\.TIF$", full.names = TRUE))
lsat2013 = rast(files2013)
lsat2013 = crop(lsat2013, delhi_ext)

# 2023 

files2023 = sort(list.files("2023", pattern = "SR_B[1-4]\\.TIF$", full.names = TRUE))
lsat2023 = rast(files2023)
lsat2023 = crop(lsat2023, delhi_ext)

# Verifica compatibilità tra raster

all.equal(res(lsat2000), res(lsat2013))  # TRUE
all.equal(res(lsat2013), res(lsat2023))  # TRUE


### Visualizazzione RGB
# Utilizziamo la funzione imageRy::im.plotRGB 
# con bande NIR (4), GREEN (2), BLUE (1) al fine di evidenziare la vegetazione.
# Più la vegetazione è attiva, più risulterà chiara.
# Bande: NIR (4), GREEN (2), BLUE (1)

# Esportazione delle immagini RGB come TIFF utilizzando imageRy::im.export

im.export(lsat2000, filename = "output_esame/RGB_2000.tif")
im.export(lsat2013, filename = "output_esame/RGB_2013.tif")
im.export(lsat2023, filename = "output_esame/RGB_2023.tif")

# utilizzo im.multiframe dal pacchetto imageRy per creare un pdf che mostri i tre plot affiancati

pdf("output_esame/RGB_comparazione_2000_2013_2023.pdf", width = 18, height = 6)
im.multiframe(1, 3)  # layout: 1 riga, 3 colonne
im.plotRGB(lsat2000, r = 4, g = 2, b = 1, title = "RGB 2000")
im.plotRGB(lsat2013, r = 4, g = 2, b = 1, title = "RGB 2013")
im.plotRGB(lsat2023, r = 4, g = 2, b = 1, title = "RGB 2023")
dev.off()


## Analisi immagini RGB

# Le immagini RGB a falsi colori (NIR, Green, Blue) di Nuova Delhi rivelano una progressiva diminuzione 
# delle aree rosse intense (vegetazione sana) tra il 2000, 2013 e 2023. Questo cambiamento cromatico,
# con l'emergere di tonalità grigio-bluastre e marroni, indica una significativa espansione dell'urbanizzazione, 
# con la conversione di spazi verdi in infrastrutture edificate. In sintesi, 
# la variazione dei colori mostra visivamente la perdita di copertura vegetale a favore dell'espansione urbana nel tempo.



### Calcolo dell'NDVI per ciascun anno tramite la funzione imageRy::im.ndvi
# NDVI = (NIR - RED) / (NIR + RED)
# NIR = banda 4, RED = banda 3


# NDVI 2000

ndvi_2000 = im.ndvi(lsat2000, nir = 4, red = 3)

# NDVI 2013

ndvi_2013 = im.ndvi(lsat2013, nir = 4, red = 3)

# NDVI 2023

ndvi_2023 = im.ndvi(lsat2023, nir = 4, red = 3)


# Esportazione NDVI come raster GeoTIFF

im.export(ndvi_2000, filename = "output_esame/NDVI_2000.tif")
im.export(ndvi_2013, filename = "output_esame/NDVI_2013.tif")
im.export(ndvi_2023, filename = "output_esame/NDVI_2023.tif")


# Salvataggio NDVI 

pdf("output_esame/NDVI_2000_2013_2023.pdf", width = 18, height = 6)
im.multiframe(1, 3)
plot(ndvi_2000, main = "NDVI 2000", col = viridis::magma(100))
plot(ndvi_2013, main = "NDVI 2013", col = viridis::magma(100))
plot(ndvi_2023, main = "NDVI 2023", col = viridis::magma(100))
dev.off()



## Analisi NDVI per 2000, 2013, 2023

# Le mappe NDVI permettono di confrontare la copertura vegetale a New Delhi 
# in tre anni chiave. Nel 2000 si osservano aree con NDVI elevato (>0.2), 
# indicative di vegetazione densa, particolarmente nella parte nord-orientale.
# Nel 2013 il valore medio NDVI si abbassa nettamente e le zone vegetate 
# si frammentano. Nel 2023, pur osservandosi una leggera ripresa visiva del contrasto, 
# l'NDVI rimane su valori molto bassi, la maggior parte <0,05, con assenza quasi totale di vegetazione vigorosa.
# Questa evoluzione è coerente con un'espansione urbana che ha ridotto 
# sia la quantità sia la qualità della vegetazione residua.


### Calcolo differenze NDVI

ndvi_diff_2000_2013 = ndvi_2000 - ndvi_2013
ndvi_diff_2013_2023 = ndvi_2013 - ndvi_2023

ndvi_diff_2000_2013 #massima diminuizione dell' NDVI: 0,5

ndvi_diff_2013_2023 #massima diminuizione dell' NDVI: 0,2

# Salvataggio Raster

im.export(ndvi_diff_2000_2013, filename = "output_esame/NDVI_diff_2000_2013.tif")
im.export(ndvi_diff_2013_2023, filename = "output_esame/NDVI_diff_2013_2023.tif")


# Prepara i dataframe per ggplot2

df_diff_1 = as.data.frame(ndvi_diff_2000_2013, xy = TRUE)
names(df_diff_1)[3] <- "diff"

df_diff_2 = as.data.frame(ndvi_diff_2013_2023, xy = TRUE)
names(df_diff_2)[3] <- "diff"

# Plot ggplot

plot_diff_1 = ggplot(df_diff_1, aes(x = x, y = y, fill = diff)) +
  geom_raster() +
  scale_fill_gradientn(colours = c("blue", "white", "red"), limits = c(-1, 1)) +
  ggtitle("Delta_NDVI 2000 -> 2013") +
  labs(x = "Easting", y = "Northing") +
  theme_minimal()

plot_diff_2 = ggplot(df_diff_2, aes(x = x, y = y, fill = diff)) +
  geom_raster() +
  scale_fill_gradientn(colours = c("blue", "white", "red"), limits = c(-1, 1)) +
  ggtitle("Delta_NDVI 2013 -> 2023") +
  labs(x = "Easting", y = "Northing") +
  theme_minimal()

# Composizione e salvataggio

diff_patch = plot_diff_1 + plot_diff_2 +
  plot_annotation(
    title = "Differenze NDVI – New Delhi",
    subtitle = "Delta_NDVI: 2000–2013 e 2013–2023",
    caption = "Rosso = perdita, Blu = guadagno"
  )

pdf("output_esame/NDVI_diff_patch.pdf", width = 10, height = 6)
print(diff_patch)
dev.off()


## ANALISI:
# La perdita tra 2000 e 2013 è massiccia, soprattutto in aree che erano vegetate (visibili in rosso).
# Tra 2013 e 2023 non si registrano guadagni, solo una leggera continuità della perdita.


### Classificazione pixel 
# Quantificare la perdita di aree verdi fra il 2000 e il 2013, dove si è registrata la maggiore differenza di NDVI
# im.classify distribuisce i pixel dell'immagine in modo non supervisionato fra il numero di classi indicato.
# Le etichette delle classi sono state assegnate confrontando visivamente le immagini NDVI con le mappe classificate:
# in questo modo è possibile associare ai gruppi spettrali le categorie ecologiche (vegetazione, cemento, suolo).

# Classificazione

class_2000 = im.classify(lsat2000, seed = 3, num_clusters = 3)
class_2013 = im.classify(lsat2013, seed = 3, num_clusters = 3)


# Salvataggio in PDF

pdf("output_esame/Classificazione_2000_2013.pdf", width = 8, height = 6)
plot(class_2000, main = "Classificazione non supervisionata – 2000")
plot(class_2013, main = "Classificazione non supervisionata – 2013")
dev.off()

# L'output della classificazione spettrale con im.classify() assegna numeri alle classi (1, 2, 3) in modo arbitrario. 
# Per interpretare correttamente i risultati, abbiamo confrontato visivamente le mappe classificate con l'NDVI:
# questo ci ha permesso di attribuire a ciascuna classe un significato ecologico coerente (vegetazione, suolo, cemento).
# Le percentuali riportate nel barplot sono state riordinate in base al contenuto informativo delle classi, 
# non al numero automatico assegnato dall’algoritmo.

# Totale pixel

tot = ncell(class_2000)
tot #3668000 coincide con lsat2000 e lsat2013

# Frequenze per classe utilizzando fre() dal pacchetto raster

library(raster)

freq_2000 = freq(raster(class_2000))
freq_2013 = freq(raster(class_2013))


# Percentuali

perc_2000 = freq_2000[ , "count"] * 100 / tot
perc_2013 = freq_2013[ , "count"] * 100 / tot

perc_2000 #vegetazione: 37.4, suolo: 40.6, cemento: 21.9
perc_2013 #vegetazione: 26,9, suolo: 31,3, cemento: 41.7


# Definizione etichette dopo confronto con NDVI (ordine: Classe 1 = Vegetazione, 2 = Suolo, 3 = Cemento)

class = c("Vegetazione", "Suolo", "Cemento")
y2000 = perc_2000
y2013 = perc_2013

tab = data.frame(class, y2000, y2013)

# Colori personalizzati per le classi

colori_classi = c("Vegetazione" = "forestgreen",
                  "Cemento" = "grey50",
                  "Suolo" = "saddlebrown")

# Barplot 2000 con colori

p1 = ggplot(tab, aes(x = class, y = y2000, fill = class, color = class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colori_classi) +
  scale_color_manual(values = colori_classi) +
  ylim(c(0, 100))

# Barplot 2013 con colori

p2 = ggplot(tab, aes(x = class, y = y2013, fill = class, color = class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colori_classi) +
  scale_color_manual(values = colori_classi) +
  ylim(c(0, 100))

# Creazione del patchwork

barplot_patch = p1 + p2 + 
  plot_annotation(
    title = "Distribuzione percentuale delle classi",
    subtitle = "New Delhi – Classificazione non supervisionata"
  )

# Visualizzazione 

print(barplot_patch)

# Salvataggio in PDF con ggsave

pdf("output_esame/class_percentuali_barplot.pdf", width = 10, height = 6)
print(barplot_patch)
dev.off()



## Analisi Classificazione
# La classificazione non supervisionata messa a confronto con l’ndvi ha evidenziato un forte aumento delle superfici cementificate
# tra il 2000 e il 2013. Contestualmente si osserva una significativa riduzione della vegetazione, 
# sia in estensione sia in vigore. I risultati evidenti nel barplot confermano la progressiva urbanizzazione 
# e la perdita di copertura vegetale già emersa dalle analisi NDVI.


### Calcolo della deviazione standard locale su finestra 3x3 pixel.

# La funzione terra::focal() applica il calcolo della SD (fun = sd) su ogni cella,
# considerando i 9 pixel (3x3) adiacenti. Questo approccio consente di misurare
# l'eterogeneità locale della vegetazione. Valori alti di SD indicano variazioni NDVI
# più marcate nello spazio -> maggiore frammentazione del paesaggio vegetale.
# Utilizziamo come raster direttamente il vettore ndvi (per il 2000 e il 2023), che è il più significativo per la vegetazione

ndvi_sd_2000 = focal(ndvi_2000, w = matrix(1/9, 3, 3), fun = sd, na.rm = TRUE)

ndvi_sd_2023 = focal(ndvi_2023, w = matrix(1/9, 3, 3), fun = sd, na.rm = TRUE)


# Esportazione dei raster SD come GeoTIFF

im.export(ndvi_sd_2000, filename = "output_esame/NDVI_SD_2000.tif")
im.export(ndvi_sd_2023, filename = "output_esame/NDVI_SD_2023.tif")


# Salvataggio SD 

pdf("output_esame/NDVI_SD_comparison.pdf", width = 12, height = 6)
im.multiframe(1, 2)
plot(ndvi_sd_2000, main = "SD NDVI – 2000", col = viridis::rocket(100))
plot(ndvi_sd_2023, main = "SD NDVI – 2023", col = viridis::rocket(100))
dev.off()


# Calcolo del valore medio di SD su tutta l'immagine (funzione terra::global)
# Questo valore sintetizza l'eterogeneità vegetativa per ogni anno in un unico indicatore.

mean_sd_2000 = global(ndvi_sd_2000, fun = mean, na.rm = TRUE)[1]
mean_sd_2023 = global(ndvi_sd_2023, fun = mean, na.rm = TRUE)[1]



cat("📊 SD media NDVI 2000:", round(mean_sd_2000[[1]], 4), "\n") #0,0031
cat("📊 SD media NDVI 2023:", round(mean_sd_2023[[1]], 4), "\n") #0,0008

## Analisi deviazione standard.

# L'SD media dell'NDVI nel 2000 era pari a 0.0031, mentre nel 2023 si è ridotta drasticamente a 0.0008.
# Questa differenza indica una forte riduzione dell'eterogeneità vegetazionale.
# In termini ecologici, questo suggerisce che nel 2000 la vegetazione urbana e periurbana era distribuita in modo più variegato,
# con la presenza di aree ad alta e bassa densità vegetativa.
# Nel 2023, invece, il paesaggio risulta molto più omogeneo, con NDVI tendenzialmente basso ovunque.
# Questo conferma una transizione urbana in cui aree verdi e vegetazione sparsa sono state rimosse o degradate,
# lasciando una copertura più uniforme, tipicamente associata a superfici artificiali come cemento o suolo nudo.


#Sintesi dei risultati – New Delhi (2000–2023)

#Visualizzazione RGB (NIR, Green, Blue)
#Le immagini a falsi colori mostrano un’evidente perdita visiva di aree vegetate tra il 2000 e il 2023.
#I toni chiari (vegetazione sana) si riducono progressivamente, sostituiti da colori grigio-bluastri o marroni, tipici di suolo o superfici artificiali.

#Vegetazione NDVI > 0.15 in forte calo
#Presente nel 2000, quasi assente nel 2023.
#La vegetazione “vigorosa” è scomparsa dalle aree urbane e periurbane.

#Delta NDVI (differenze)
#2000–2013: perdita massiccia (fino a -0.5 NDVI).
#2013–2023: perdita più contenuta (fino a -0.2 NDVI), ma ancora attiva.

#Classificazione spettrale 
#Classe con NDVI medio più alto (0.21) nel 2000 → nel 2013 scesa a 0.02.
#Apparente crescita dell’area verde, ma con valori NDVI molto più bassi → vegetazione degradata o non funzionale.

#Eterogeneità NDVI 
#SD media 2000 = 0.0031 → nel 2023 = 0.0008.
#Il paesaggio diventa più omogeneo, segno di semplificazione ecologica (povera varietà vegetativa).

#Conclusione
#La copertura vegetale è diminuita sia in quantità che in qualità.
#L’urbanizzazione ha portato a un paesaggio più uniforme, con meno biodiversità e meno contrasto vegetativo.
