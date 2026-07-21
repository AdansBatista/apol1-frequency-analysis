"""Draw a schematic West Africa-to-Americas historical context map."""

import json
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Polygon


THIS = Path(__file__).resolve()
ROOT = THIS.parent
REFERENCE_DIR = ROOT / "data" / "reference"
GEOJSON_PATH = REFERENCE_DIR / "ne_110m_admin_0_countries.geojson"
OUTPUT = ROOT / "figures" / "west_africa_to_americas_map.png"
NATURAL_EARTH_URL = (
    "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/"
    "geojson/ne_110m_admin_0_countries.geojson"
)


SOURCE_SITES = [
    ("GWD", "The Gambia", -15.31, 13.45),
    ("MSL", "Sierra Leone", -11.78, 8.46),
    ("ESN / YRI", "Nigeria", 7.49, 9.06),
]

AMERICAN_SITES = [
    ("ACB", "Barbados", -59.54, 13.19),
    ("ASW", "United States", -99.0, 31.0),
    ("ABraOM", "Sao Paulo", -46.63, -23.55),
]


def load_boundaries():
    """Download and cache the Natural Earth country boundaries."""
    REFERENCE_DIR.mkdir(parents=True, exist_ok=True)
    if not GEOJSON_PATH.exists():
        request = urllib.request.Request(
            NATURAL_EARTH_URL,
            headers={"User-Agent": "BIOL625-APOL1-map/1.0"},
        )
        with urllib.request.urlopen(request, timeout=60) as response:
            GEOJSON_PATH.write_bytes(response.read())
    return json.loads(GEOJSON_PATH.read_text(encoding="utf-8"))


def iter_outer_rings(geometry):
    """Yield outer polygon rings from Polygon and MultiPolygon features."""
    if geometry["type"] == "Polygon":
        yield geometry["coordinates"][0]
    elif geometry["type"] == "MultiPolygon":
        for polygon in geometry["coordinates"]:
            yield polygon[0]


def country_color(name):
    """Color countries that contain study reference or comparison sites."""
    if name in {"Nigeria", "Gambia", "Sierra Leone"}:
        return "#efc66a"
    if name == "Brazil":
        return "#72aa8a"
    if name == "United States of America":
        return "#8fb7cc"
    return "#f2f0e9"


def add_route(axis, start, end, curvature):
    """Add one broad, schematic transatlantic route arrow."""
    arrow = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        connectionstyle=f"arc3,rad={curvature}",
        mutation_scale=15,
        linewidth=2.0,
        color="#9f4b3f",
        alpha=0.76,
        zorder=4,
    )
    axis.add_patch(arrow)


def add_site(axis, code, label, longitude, latitude, color, offset):
    """Add a study site marker and compact label."""
    axis.scatter(
        longitude,
        latitude,
        s=62,
        color=color,
        edgecolor="white",
        linewidth=1.2,
        zorder=7,
    )
    axis.annotate(
        f"{code}\n{label}",
        (longitude, latitude),
        xytext=offset,
        textcoords="offset points",
        fontsize=8,
        fontweight="bold",
        ha="left" if offset[0] >= 0 else "right",
        va="center",
        zorder=8,
        bbox={
            "boxstyle": "round,pad=0.22",
            "facecolor": "white",
            "edgecolor": "none",
            "alpha": 0.84,
        },
    )


def main():
    boundaries = load_boundaries()
    figure, axis = plt.subplots(figsize=(11, 6.5))
    figure.patch.set_facecolor("white")
    axis.set_facecolor("#dcebf0")

    for feature in boundaries["features"]:
        name = feature["properties"].get("ADMIN", "")
        for ring in iter_outer_rings(feature["geometry"]):
            axis.add_patch(
                Polygon(
                    ring,
                    closed=True,
                    facecolor=country_color(name),
                    edgecolor="#aaa89f",
                    linewidth=0.45,
                    zorder=1,
                )
            )

    route_origin = (-7.0, 8.0)
    add_route(axis, route_origin, (-38.50, -12.97), 0.18)
    add_route(axis, route_origin, (-59.54, 13.19), -0.12)
    add_route(axis, route_origin, (-79.0, 27.0), -0.22)

    source_offsets = {
        "GWD": (-10, 18),
        "MSL": (-8, -20),
        "ESN / YRI": (8, 14),
    }
    for code, label, longitude, latitude in SOURCE_SITES:
        add_site(
            axis,
            code,
            label,
            longitude,
            latitude,
            "#c98716",
            source_offsets[code],
        )

    comparison_offsets = {
        "ACB": (8, 13),
        "ASW": (9, 3),
        "ABraOM": (9, -2),
    }
    for code, label, longitude, latitude in AMERICAN_SITES:
        add_site(
            axis,
            code,
            label,
            longitude,
            latitude,
            "#176b87",
            comparison_offsets[code],
        )

    axis.scatter(
        -38.50,
        -12.97,
        s=62,
        facecolor="white",
        edgecolor="#176b87",
        linewidth=2,
        zorder=7,
    )
    axis.annotate(
        "Salvador / Bahia\nregional context",
        (-38.50, -12.97),
        xytext=(9, 9),
        textcoords="offset points",
        fontsize=8,
        ha="left",
        va="bottom",
        bbox={
            "boxstyle": "round,pad=0.22",
            "facecolor": "white",
            "edgecolor": "none",
            "alpha": 0.84,
        },
        zorder=8,
    )

    axis.text(
        -42,
        25,
        "Schematic historical\nforced-migration routes",
        color="#84382f",
        fontsize=10,
        fontweight="bold",
        ha="center",
        zorder=8,
    )
    axis.set_xlim(-112, 25)
    axis.set_ylim(-39, 39)
    axis.set_aspect("equal", adjustable="box")
    axis.axis("off")
    axis.set_title(
        "West African reference sites and historical context in the Americas",
        fontsize=16,
        fontweight="bold",
        pad=13,
    )

    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="#c98716",
            markeredgecolor="white",
            markersize=8,
            label="1000 Genomes West African site",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="#176b87",
            markeredgecolor="white",
            markersize=8,
            label="American comparison site",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="white",
            markeredgecolor="#176b87",
            markeredgewidth=1.5,
            markersize=8,
            label="Regional context only",
        ),
        Line2D(
            [0],
            [0],
            color="#9f4b3f",
            linewidth=2,
            label="Broad historical route",
        ),
    ]
    axis.legend(
        handles=legend_handles,
        loc="lower left",
        frameon=True,
        framealpha=0.92,
        fontsize=8,
        ncol=2,
    )
    figure.text(
        0.5,
        0.02,
        "Arrows provide historical context only. They do not represent individual ancestry or exact voyage counts.",
        ha="center",
        fontsize=9,
        color="#444444",
    )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(figure)
    print(f"[data] {GEOJSON_PATH}")
    print(f"[figure] {OUTPUT}")


if __name__ == "__main__":
    main()