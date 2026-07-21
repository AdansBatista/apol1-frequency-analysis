"""Draw APOL1 G1/G2 changes in the canonical UniProt O14791 sequence."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch


THIS = Path(__file__).resolve()
ROOT = THIS.parent
OUTPUT = ROOT / "figures" / "apol1_g1_g2_alignment.png"

REGION_START = 338
REFERENCE_REGION = "VAPVSFFLVLDVVYLVYESKHLHEGAKSETAEELKKVAQELEEKLNILNNNYKILQ"


def residue(position):
    """Return one canonical residue by its one-based APOL1 position."""
    return REFERENCE_REGION[position - REGION_START]


def draw_sequence(axis, positions, rows, changes, title):
    """Draw aligned amino-acid boxes and emphasize variant positions."""
    axis.set_xlim(-3.5, len(positions) + 0.2)
    axis.set_ylim(-0.4, len(rows) + 1.15)
    axis.axis("off")
    axis.set_title(title, loc="left", fontsize=12, fontweight="bold", pad=10)

    for column, position in enumerate(positions):
        axis.text(
            column + 0.5,
            len(rows) + 0.55,
            str(position),
            ha="center",
            va="center",
            fontsize=7,
            color="#555555",
            rotation=45,
        )

    for row_index, (row_label, sequence) in enumerate(rows.items()):
        y_position = len(rows) - row_index - 0.35
        axis.text(
            -0.25,
            y_position + 0.42,
            row_label,
            ha="right",
            va="center",
            fontsize=10,
            fontweight="bold" if row_label != "Reference" else "normal",
        )
        for column, amino_acid in enumerate(sequence):
            change_type = changes.get((row_label, positions[column]))
            if change_type == "substitution":
                face_color, text_color, edge_color = "#c94c36", "white", "#9c2f20"
            elif change_type == "deletion":
                face_color, text_color, edge_color = "#fff1d6", "#9b5d00", "#d3942f"
            elif row_label == "Reference" and positions[column] in {342, 384, 388, 389}:
                face_color, text_color, edge_color = "#dfe9ec", "#17363f", "#78929a"
            else:
                face_color, text_color, edge_color = "#f7f7f5", "#222222", "#c5c5c0"
            box = FancyBboxPatch(
                (column + 0.08, y_position),
                0.84,
                0.84,
                boxstyle="round,pad=0.02,rounding_size=0.04",
                facecolor=face_color,
                edgecolor=edge_color,
                linewidth=1,
            )
            axis.add_patch(box)
            axis.text(
                column + 0.5,
                y_position + 0.42,
                amino_acid,
                ha="center",
                va="center",
                fontsize=11,
                family="DejaVu Sans Mono",
                fontweight="bold",
                color=text_color,
            )


def main():
    assert len(REFERENCE_REGION) == 56
    assert residue(342) == "S"
    assert residue(384) == "I"
    assert residue(388) == "N"
    assert residue(389) == "Y"

    first_positions = list(range(338, 349))
    first_reference = "".join(residue(position) for position in first_positions)
    first_g1 = list(first_reference)
    first_g1[first_positions.index(342)] = "G"

    second_positions = list(range(380, 394))
    second_reference = "".join(residue(position) for position in second_positions)
    second_g1 = list(second_reference)
    second_g1[second_positions.index(384)] = "M"
    second_g2 = list(second_reference)
    second_g2[second_positions.index(388)] = "-"
    second_g2[second_positions.index(389)] = "-"

    figure, axes = plt.subplots(2, 1, figsize=(11, 6.3), gridspec_kw={"hspace": 0.62})
    figure.suptitle(
        "APOL1 G1 substitutions and G2 in-frame deletion",
        fontsize=17,
        fontweight="bold",
        y=0.99,
    )

    draw_sequence(
        axes[0],
        first_positions,
        {"Reference": first_reference, "G1": "".join(first_g1)},
        {("G1", 342): "substitution"},
        "A. G1 first substitution: S342G (serine to glycine)",
    )
    draw_sequence(
        axes[1],
        second_positions,
        {
            "Reference": second_reference,
            "G1": "".join(second_g1),
            "G2": "".join(second_g2),
        },
        {
            ("G1", 384): "substitution",
            ("G2", 388): "deletion",
            ("G2", 389): "deletion",
        },
        "B. G1 second substitution: I384M; G2 deletion: N388-Y389del",
    )

    figure.text(
        0.5,
        0.015,
        "G1 contains both substitutions on the same chromosome. G2 removes six DNA bases and two amino acids.",
        ha="center",
        fontsize=10,
        color="#444444",
    )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(figure)
    print(f"[figure] {OUTPUT}")


if __name__ == "__main__":
    main()