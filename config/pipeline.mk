# Paramètres de la chaîne de traitement (versionnés ; les modifier relance les produits
# qui en dépendent au prochain « make »).

# Dossier de compilation (Release)
BUILD ?= build/Release
# Degré maximal des harmoniques sphériques
LMAX = 30
# |latitude| maximale des pixels ajustés (°) ; 90 = carte entière
MAX_LAT = 90
# Rayon de la surface source (rayons solaires)
RSS = 2.5
# En dessous de ce nombre d'heures valides, la fenêtre est marquée insufficient_data
MIN_VALID_HOURS = 48
