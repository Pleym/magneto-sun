# Chaîne de traitement magneto-sun : exploitation d'une campagne de fenêtres d'observation.
#
#   make build          compile (Release) et lance les tests
#   make fetch          télécharge les entrées L1 (réseau : frontale ou poste local)
#   make campaign       produits L2 et L3 de chaque fenêtre, puis synthèse L4
#   make replay         vitrine : rejeu animé de la campagne (figures/replay.gif, .mp4)
#   make status         état de chaque fenêtre
#   make clean-products supprime les produits (les entrées L1 restent)
#
# Chaque produit dépend de ses entrées, des paramètres (config/pipeline.mk) et du code
# compilé : seul ce qui est périmé est recalculé. « make -j N » traite N étapes en
# parallèle ; « make -k » poursuit après l'échec d'une fenêtre et le signale.
# Niveaux : L1 entrées brutes (data/), L2 standardisés, L3 scientifiques, L4 campagne.

include config/pipeline.mk

CAMPAIGN ?= config/campaign_solo_2020_2022.txt
PRODUCTS ?= products
STAMP := $(BUILD)/.stamp
PROV := scripts/provenance.sh
SOURCES := CMakeLists.txt $(wildcard src/cpp/*.cpp src/cpp/*.hpp src/fortran/*.f90)

WINDOWS := $(shell awk '!/^\#/ && NF {print $$1}' $(CAMPAIGN))

# Valeurs de la table de campagne pour une fenêtre : $(call start,w01), etc.
field = $(shell awk '$$1 == "$(1)" {print $$$(2)}' $(CAMPAIGN))
start = $(call field,$(1),2)
stop = $(call field,$(1),3)
hours = $(call field,$(1),4)
gong = data/$(call field,$(1),5)
tag = $(firstword $(subst T, ,$(call start,$(1))))_$(firstword $(subst T, ,$(call stop,$(1))))
mag = data/solo_mag_rtn_$(call tag,$(1)).txt
pas = data/solo_pas_v_$(call tag,$(1)).txt
eph = data/solo_horizons_$(call tag,$(1)).txt
log = $(PRODUCTS)/$(1)/logs/$(2).log

SCORES := $(foreach w,$(WINDOWS),$(PRODUCTS)/$(w)/L3/polarity_score.txt)
SURFACES := $(foreach w,$(WINDOWS),$(PRODUCTS)/$(w)/L3/pfss_ss.bin)

.PHONY: all build fetch campaign replay status clean-products
.DELETE_ON_ERROR:
# Les produits L2 et L3 sont des résultats, pas des fichiers intermédiaires à effacer
.SECONDARY:
.SECONDEXPANSION:

all: campaign

# --- Logiciel ------------------------------------------------------------------------
$(STAMP): $(SOURCES)
	cmake -S . -B $(BUILD) -DCMAKE_BUILD_TYPE=Release
	cmake --build $(BUILD) -j
	touch $@

build: $(STAMP)
	ctest --test-dir $(BUILD) --output-on-failure

# --- L1 : entrées brutes (réseau) ------------------------------------------------------
fetch: $(foreach w,$(WINDOWS),$(call gong,$(w)) $(call mag,$(w)) $(call pas,$(w)) $(call eph,$(w)))

data/mrzqs%.fits.gz:
	scripts/fetch_gong.sh mrzqs$*.fits.gz

# Les trois fichiers in situ d'une fenêtre viennent d'un même téléchargement
data/solo_mag_rtn_%.txt data/solo_pas_v_%.txt data/solo_horizons_%.txt:
	scripts/fetch_solo.sh $(word 1,$(subst _, ,$*))T00:00:00 $(word 2,$(subst _, ,$*))T00:00:00

# --- L2 : série horaire in situ standardisée ------------------------------------------
$(PRODUCTS)/%/L2/insitu_hourly.txt: $$(call mag,$$*) $$(call pas,$$*) $$(call eph,$$*) \
                                    $(STAMP) $(CAMPAIGN)
	@mkdir -p $(@D) $(PRODUCTS)/$*/logs
	$(BUILD)/insitu_l2 $(call mag,$*) $(call pas,$*) $(call eph,$*) $(call start,$*) \
	  $(call hours,$*) $@ > $(call log,$*,insitu_l2)
	$(PROV) $@ L2 $(BUILD)/insitu_l2 "start=$(call start,$*) hours=$(call hours,$*)" \
	  $(call mag,$*) $(call pas,$*) $(call eph,$*)

# --- L3 : produits scientifiques ------------------------------------------------------
$(PRODUCTS)/%/L3/fit_coeffs.txt: $$(call gong,$$*) $(STAMP) config/pipeline.mk
	@mkdir -p $(@D) $(PRODUCTS)/$*/logs
	$(BUILD)/fit_map $< $(LMAX) $(MAX_LAT) $(PRODUCTS)/$*/L3/fit > $(call log,$*,fit_map)
	$(PROV) $@ L3 $(BUILD)/fit_map "lmax=$(LMAX) max_lat=$(MAX_LAT)" $<

$(PRODUCTS)/%/L3/pfss_ss.bin: $(PRODUCTS)/%/L3/fit_coeffs.txt $(STAMP) config/pipeline.mk
	$(BUILD)/pfss_map $< $(RSS) $(PRODUCTS)/$*/L3/pfss > $(call log,$*,pfss_map)
	$(PROV) $@ L3 $(BUILD)/pfss_map "rss=$(RSS)" $<

$(PRODUCTS)/%/L3/polarity_score.txt: $(PRODUCTS)/%/L3/fit_coeffs.txt \
                                     $(PRODUCTS)/%/L2/insitu_hourly.txt $(STAMP) config/pipeline.mk
	$(BUILD)/polarity $(word 1,$^) $(RSS) $(word 2,$^) $(MIN_VALID_HOURS) \
	  $(PRODUCTS)/$*/L3/polarity > $(call log,$*,polarity)
	$(PROV) $@ L3 $(BUILD)/polarity "rss=$(RSS) min_valid_hours=$(MIN_VALID_HOURS)" \
	  $(word 1,$^) $(word 2,$^)

# --- L4 : synthèse de campagne ---------------------------------------------------------
$(PRODUCTS)/campaign_summary.txt: $(SCORES) scripts/campaign_summary.sh
	scripts/campaign_summary.sh $(SCORES) > $@
	$(PROV) $@ L4 scripts/campaign_summary.sh "campaign=$(CAMPAIGN)" $(SCORES)

$(PRODUCTS)/campaign_summary.png: $(PRODUCTS)/campaign_summary.txt plots/campaign.gp
	gnuplot -e "in='$<'; out='$@'" plots/campaign.gp

campaign: $(PRODUCTS)/campaign_summary.png $(SURFACES)

# --- Vitrine : rejeu de la campagne (GIF + MP4) et figure de synthèse versionnée -------
figures/replay.gif: $(PRODUCTS)/campaign_summary.txt $(SURFACES) scripts/make_replay.sh \
                    plots/replay_frame.gp
	scripts/make_replay.sh $(CAMPAIGN) $(PRODUCTS) figures/replay

figures/campaign_summary.png: $(PRODUCTS)/campaign_summary.png
	cp $< $@

replay: figures/replay.gif figures/campaign_summary.png

# --- Exploitation ----------------------------------------------------------------------
status:
	@for w in $(WINDOWS); do \
	  score=$(PRODUCTS)/$$w/L3/polarity_score.txt; \
	  if [ -f $$score ]; then \
	    awk -v w=$$w '{v[$$1] = $$2} END {printf "%s  %-18s %-17s accord %s %%\n", w, v["start"], v["status"], v["sector_agreement"]}' $$score; \
	  else \
	    echo "$$w  (pas encore traitée)"; \
	  fi; \
	done

clean-products:
	rm -rf $(PRODUCTS)
