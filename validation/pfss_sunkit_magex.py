"""Validation externe du modèle PFSS (Phase 3).

Compare notre B_r sur la surface source (sortie *_ss.bin de pfss_map) à celui de
sunkit-magex, solveur numérique indépendant, calculé sur la même carte GONG. La carte
est relue par sunpy : la remise en ordre des longitudes de notre lecteur est donc
vérifiée elle aussi. Seul code Python du projet ; il ne sert qu'à cette comparaison.

Usage :
  validation/.venv/bin/python validation/pfss_sunkit_magex.py \
      data/mrzqs200727t0004c2233_145.fits.gz data/cr2233_full_ss.bin 2.5 data/cr2233_sunkit
Écrit <préfixe>_ss.bin (B_r de sunkit-magex sur notre grille) et <préfixe>_neutral.txt.
"""
import sys

import astropy.units as u
import numpy as np
import sunpy.map
from astropy.coordinates import SkyCoord
from sunkit_magex import pfss

N_LON, N_LAT = 360, 180
N_R = 50  # cellules radiales du solveur de sunkit-magex


def canonical_grid():
    """Centres des pixels de notre grille : longitude (°), latitude (°), forme (N_LAT, N_LON)."""
    lon = (np.arange(N_LON) + 0.5) * 360 / N_LON
    sinlat = -1 + (np.arange(N_LAT) + 0.5) * 2 / N_LAT
    return np.meshgrid(lon, np.degrees(np.arcsin(sinlat)))


def sample_on_canonical_grid(m):
    """Valeurs de la carte sunpy m aux centres de notre grille (pixel le plus proche)."""
    lon, lat = canonical_grid()
    x, y = m.world_to_pixel(SkyCoord(lon * u.deg, lat * u.deg, frame=m.coordinate_frame))
    i = np.round(x.value).astype(int) % m.data.shape[1]
    j = np.clip(np.round(y.value).astype(int), 0, m.data.shape[0] - 1)
    return m.data[j, i]


def neutral_line(br):
    """Zéros de B_r entre pixels voisins, comme pfss_map : (longitude °, sin(latitude))."""
    lon = (np.arange(N_LON) + 0.5) * 360 / N_LON
    sinlat = -1 + (np.arange(N_LAT) + 0.5) * 2 / N_LAT
    east = np.roll(br, -1, axis=1)
    j, i = np.nonzero(br * east < 0)
    points = [np.column_stack([lon[i] + br[j, i] / (br[j, i] - east[j, i]) * 360 / N_LON, sinlat[j]])]
    j, i = np.nonzero(br[:-1] * br[1:] < 0)
    t = br[j, i] / (br[j, i] - br[j + 1, i])
    points.append(np.column_stack([lon[i], sinlat[j] + t * 2 / N_LAT]))
    return np.vstack(points)


def main(gong_path, ours_path, rss, prefix):
    gong = sunpy.map.Map(gong_path)
    theirs = sample_on_canonical_grid(pfss.pfss(pfss.Input(gong, N_R, rss)).source_surface_br)
    ours = np.fromfile(ours_path, dtype="<f8").reshape(N_LAT, N_LON)

    print(f"sunkit-magex, nr = {N_R}, Rss = {rss} R_sun")
    print(f"même signe de B_r(Rss)     : {100 * np.mean(np.sign(ours) == np.sign(theirs)):.1f} % de l'aire")
    print(f"corrélation                : {np.corrcoef(ours.ravel(), theirs.ravel())[0, 1]:.4f}")
    print(f"rapport des écarts-types   : {np.std(ours) / np.std(theirs):.3f} (nous / sunkit-magex)")

    points = neutral_line(theirs)
    points[:, 0] %= 360
    theirs.astype("<f8").tofile(prefix + "_ss.bin")
    np.savetxt(prefix + "_neutral.txt", points,
               header="ligne neutre sunkit-magex : longitude (°) sin(latitude)")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4])
