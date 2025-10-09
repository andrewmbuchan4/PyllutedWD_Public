#!/usr/bin/env python

# import csv
# import numpy as np
# import os

import graph_factory as gf
import pwd_utils as pu

dump_config1 = {
    "y_dimension": 3,
    "x_dimension": 3,
    "filenames": ["hollands_ksp_multiplots.pdf", "hollands_ksp_multiplots.png"],
    "fig_height": 10,
    "fig_width": 15,
    "gridspec_wspace": 0,
    "gridspec_hspace": 0,
    "sharey_axes": True,
    "sharex_axes": True,
    "list_of_dumps": [
        pu.get_path_to_pipeline_base_dir()
        + "deltavarerr/pipeline_ksp_comparison_SyntheticHollandsTidal_Hollands"
        + "_Ca_Fe_fragment_core_frac.pdf.txt",
        pu.get_path_to_pipeline_base_dir()
        + "deltavarerr/pipeline_ksp_comparison_SyntheticHollandsTidal_Hollands"
        + "_Mg_Fe_fragment_core_frac.pdf.txt",
        pu.get_path_to_pipeline_base_dir()
        + "deltavarerr/pipeline_ksp_comparison_SyntheticHollandsTidal_Hollands"
        + "_Ca_Mg_fragment_core_frac.pdf.txt",
        pu.get_path_to_pipeline_base_dir()
        + "deltavarerr/pipeline_ksp_comparison_SyntheticHollandsCollisional_Hollands"
        + "_Ca_Fe_fragment_core_frac.pdf.txt",
        pu.get_path_to_pipeline_base_dir()
        + "deltavarerr/pipeline_ksp_comparison_SyntheticHollandsCollisional_Hollands"
        + "_Mg_Fe_fragment_core_frac.pdf.txt",
        pu.get_path_to_pipeline_base_dir()
        + "deltavarerr/pipeline_ksp_comparison_SyntheticHollandsCollisional_Hollands"
        + "_Ca_Mg_fragment_core_frac.pdf.txt",
        pu.get_path_to_pipeline_base_dir()
        + "deltavarerr/pipeline_ksp_comparison_SyntheticHollandsDeltaPop_Hollands"
        + "_Ca_Fe_fragment_core_frac.pdf.txt",
        pu.get_path_to_pipeline_base_dir()
        + "deltavarerr/pipeline_ksp_comparison_SyntheticHollandsDeltaPop_Hollands"
        + "_Mg_Fe_fragment_core_frac.pdf.txt",
        pu.get_path_to_pipeline_base_dir()
        + "deltavarerr/pipeline_ksp_comparison_SyntheticHollandsDeltaPop_Hollands"
        + "_Ca_Mg_fragment_core_frac.pdf.txt",
    ],
}

dump_config2 = {
    "y_dimension": 1,
    "x_dimension": 3,
    "filenames": [
        "HollandsDelta0p1_CaFe_MgFe_CaMg.pdf",
        "HollandsDelta0p1_CaFe_MgFe_CaMg.png",
    ],
    "fig_height": 10,
    "fig_width": 15,
    "gridspec_wspace": 0,
    "gridspec_hspace": 0,
    "sharey_axes": True,
    "sharex_axes": False,
    "list_of_dumps": [
        pu.get_path_to_pipeline_base_dir()
        + "deltavarerr"
        + "/SyntheticHollandsDeltaPop_RealisticObservererr0p1_StandardModeller"
        + "/fcf_Observed_pdf_Ca_Fe.pdf.txt",
        pu.get_path_to_pipeline_base_dir()
        + "deltavarerr"
        + "/SyntheticHollandsDeltaPop_RealisticObservererr0p1_StandardModeller"
        + "/fcf_Observed_pdf_Mg_Fe.pdf.txt",
        pu.get_path_to_pipeline_base_dir()
        + "deltavarerr"
        + "/SyntheticHollandsDeltaPop_RealisticObservererr0p1_StandardModeller"
        + "/fcf_Observed_pdf_Ca_Mg.pdf.txt",
    ],
}

dump_config3 = {
    "y_dimension": 3,
    "x_dimension": 3,
    "filenames": ["ternaryplots.pdf", "ternaryplots.png"],
    "fig_height": 10,
    "fig_width": 15,
    "gridspec_wspace": 0.2,
    "gridspec_hspace": 0.2,
    "sharey_axes": False,
    "sharex_axes": False,
    "list_of_dumps": [
        pu.get_path_to_pipeline_base_dir()
        + "dzcomparison/SyntheticHollandsDeltaPop_HollandsObservererr0_NullModeller"
        + "/ternary_plot_sampled_SyntheticHollandsDeltaPop_HollandsObservererr0.pdf"
        + ".txt",
        pu.get_path_to_pipeline_base_dir()
        + "dzcomparison/SyntheticHollandsDeltaPop_HollandsObservererr0p2_NullModeller"
        + "/ternary_plot_sampled_SyntheticHollandsDeltaPop_HollandsObservererr0p2.pdf"
        + ".txt",
        pu.get_path_to_pipeline_base_dir()
        + "dzcomparison/SyntheticHollandsDeltaPop_HollandsObservererr0p4_NullModeller"
        + "/ternary_plot_sampled_SyntheticHollandsDeltaPop_HollandsObservererr0p4.pdf"
        + ".txt",
        pu.get_path_to_pipeline_base_dir()
        + "dzcomparison/SyntheticHollandsCollisional_HollandsObservererr0_NullModeller"
        + "/ternary_plot_sampled_SyntheticHollandsCollisional_HollandsObservererr0.pdf"
        + ".txt",
        pu.get_path_to_pipeline_base_dir()
        + "dzcomparison/SyntheticHollandsCollisional_HollandsObservererr0p2_NullModeller"
        + "/ternary_plot_sampled_SyntheticHollandsCollisional_HollandsObservererr0p2"
        + ".pdf.txt",
        pu.get_path_to_pipeline_base_dir()
        + "dzcomparison/SyntheticHollandsCollisional_HollandsObservererr0p4_NullModeller"
        + "/ternary_plot_sampled_SyntheticHollandsCollisional_HollandsObservererr0p4"
        + ".pdf.txt",
        pu.get_path_to_pipeline_base_dir()
        + "dzcomparison/SyntheticHollandsTidal_HollandsObservererr0_NullModeller"
        + "/ternary_plot_sampled_SyntheticHollandsTidal_HollandsObservererr0.pdf.txt",
        pu.get_path_to_pipeline_base_dir()
        + "dzcomparison/SyntheticHollandsTidal_HollandsObservererr0p2_NullModeller"
        + "/ternary_plot_sampled_SyntheticHollandsTidal_HollandsObservererr0p2.pdf.txt",
        pu.get_path_to_pipeline_base_dir()
        + "dzcomparison/SyntheticHollandsTidal_HollandsObservererr0p4_NullModeller"
        + "/ternary_plot_sampled_SyntheticHollandsTidal_HollandsObservererr0p4.pdf.txt",
    ],
}

dump_config4 = {
    "y_dimension": 4,
    "x_dimension": 3,
    "filenames": [
        "da_overshoot_timescale_comparisons_mk1.pdf",
        "da_overshoot_timescale_comparisons_mk1.png",
    ],
    "fig_height": 10,
    "fig_width": 15,
    "gridspec_wspace": 0.2,
    "gridspec_hspace": 0,
    "sharey_axes": False,
    "sharex_axes": True,
    "list_of_dumps": [
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_C_Ca_bo_bn_H_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_C_Ca_vo_bn_H_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_C_Ca_3o_bn_H_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_O_Ca_bo_bn_H_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_O_Ca_vo_bn_H_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_O_Ca_3o_bn_H_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_Fe_Ca_bo_bn_H_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_Fe_Ca_vo_bn_H_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_Fe_Ca_3o_bn_H_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_Mg_Ca_bo_bn_H_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_Mg_Ca_vo_bn_H_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_Mg_Ca_3o_bn_H_SS.pdf.txt",
    ],
}

dump_config5 = {
    "y_dimension": 4,
    "x_dimension": 3,
    "filenames": [
        "da_overshoot_timescale_comparisons.pdf",
        "da_overshoot_timescale_comparisons.png",
    ],
    "fig_height": 10,
    "fig_width": 15,
    "gridspec_wspace": 0.2,
    "gridspec_hspace": 0,
    "sharey_axes": False,
    "sharex_axes": True,
    "list_of_dumps": [
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_Ca_Mg_bo_bn_H_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_Ca_Mg_vo_bn_H_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_Ca_Mg_3o_bn_H_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_Fe_Mg_bo_bn_H_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_Fe_Mg_vo_bn_H_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_Fe_Mg_3o_bn_H_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_Si_Mg_bo_bn_H_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_Si_Mg_vo_bn_H_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_Si_Mg_3o_bn_H_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_O_Mg_bo_bn_H_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_O_Mg_vo_bn_H_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_O_Mg_3o_bn_H_SS.pdf.txt",
    ],
}

dump_config6 = {
    "y_dimension": 4,
    "x_dimension": 2,
    "filenames": [
        "db_overshoot_timescale_comparisons.pdf",
        "db_overshoot_timescale_comparisons.png",
    ],
    "fig_height": 13,
    "fig_width": 11,
    "gridspec_wspace": 0.2,
    "gridspec_hspace": 0,
    "sharey_axes": False,
    "sharex_axes": True,
    "list_of_dumps": [
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_Ca_Mg_bo_bn_He_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_Ca_Mg_vo_bn_He_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_Fe_Mg_bo_bn_He_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_Fe_Mg_vo_bn_He_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_Si_Mg_bo_bn_He_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_Si_Mg_vo_bn_He_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_O_Mg_bo_bn_He_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_O_Mg_vo_bn_He_SS.pdf.txt",
    ],
}

dump_config7 = {
    "y_dimension": 4,
    "x_dimension": 2,
    "filenames": [
        "da_bvk_timescale_comparisons.pdf",
        "da_bvk_timescale_comparisons.png",
    ],
    "fig_height": 13,
    "fig_width": 11,
    "gridspec_wspace": 0.2,
    "gridspec_hspace": 0,
    "sharey_axes": False,
    "sharex_axes": True,
    "list_of_dumps": [
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_Ca_Mg_ko_bo_H_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_Ca_Mg_kn_bn_H_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_Fe_Mg_ko_bo_H_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_Fe_Mg_kn_bn_H_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_Si_Mg_ko_bo_H_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_Si_Mg_kn_bn_H_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_O_Mg_ko_bo_H_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_O_Mg_kn_bn_H_SS.pdf.txt",
    ],
}

dump_config8 = {
    "y_dimension": 4,
    "x_dimension": 2,
    "filenames": [
        "db_bvk_timescale_comparisons.pdf",
        "db_bvk_timescale_comparisons.png",
    ],
    "fig_height": 13,
    "fig_width": 11,
    "gridspec_wspace": 0.2,
    "gridspec_hspace": 0,
    "sharey_axes": False,
    "sharex_axes": True,
    "list_of_dumps": [
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_Ca_Mg_ko_bo_He_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_Ca_Mg_kn_bn_He_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_Fe_Mg_ko_bo_He_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_Fe_Mg_kn_bn_He_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_Si_Mg_ko_bo_He_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_Si_Mg_kn_bn_He_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_O_Mg_ko_bo_He_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_O_Mg_kn_bn_He_SS.pdf.txt",
    ],
}

dump_config9 = {
    "y_dimension": 1,
    "x_dimension": 2,
    "filenames": [
        "da_overshoot_timescale_comparison.pdf",
        "da_overshoot_timescale_comparison.png",
    ],
    "fig_height": 13,
    "fig_width": 11,
    "gridspec_wspace": 0.2,
    "gridspec_hspace": 0,
    "sharey_axes": False,
    "sharex_axes": True,
    "list_of_dumps": [
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_Fe_Mg_3p_bn_H_SS.pdf.txt",
        pu.get_path_to_default_graphs()
        + "timescale_type_comparison_O_Mg_3p_bn_H_SS.pdf.txt",
    ],
}

dump_config10 = {
    "y_dimension": 2,
    "x_dimension": 2,
    "filenames": ["robust_systems.pdf", "robust_systems.png"],
    "fig_height": 11,
    "fig_width": 14,
    "gridspec_wspace": 0,
    "gridspec_hspace": 0,
    "sharey_axes": True,
    "sharex_axes": True,
    "list_of_dumps": [
        pu.get_path_to_default_graphs()
        + "../../../"
        + "Graphs of Robust Systems"
        + "/WDJ183352.68+321757.25"
        + "_TremblayGR_kn_n_p2000_HD0123_NEL_D_composition_rel_Mg.pdf.txt",
        pu.get_path_to_default_graphs()
        + "../../../"
        + "Graphs of Robust Systems"
        + "/HE0106-3253_XuGR_3p_t_p2000_HD013_NEL_D_composition_rel_Mg.pdf.txt",
        pu.get_path_to_default_graphs()
        + "../../../"
        + "Graphs of Robust Systems"
        + "/PG1225-079Model2_KleinGR_kn_n_p2000_HD01_NEL_D_composition_rel_Mg.pdf.txt",
        pu.get_path_to_default_graphs()
        + "../../../"
        + "Graphs of Robust Systems"
        + "/GD133_XuGR_3p_t_p2000_HD01_NEL_D_composition_rel_Mg.pdf.txt",
    ],
}

dump_config11 = {
    "y_dimension": 3,
    "x_dimension": 2,
    "filenames": ["timescale_crosssections.pdf", "timescale_crosssections.png"],
    "fig_height": 19,
    "fig_width": 14,
    "gridspec_wspace": 0.1,
    "gridspec_hspace": 0.1,
    "sharey_axes": False,
    "sharex_axes": False,
    "list_of_dumps": [
        pu.get_path_to_default_graphs() + "timescales_v_Teff_H_3p_bn.pdf.txt",
        pu.get_path_to_default_graphs() + "timescales_v_Teff_He_vo_bn.pdf.txt",
        pu.get_path_to_default_graphs() + "timescales_v_Teff_H_relMg_3p.pdf.txt",
        pu.get_path_to_default_graphs() + "timescales_v_Teff_H_bn_kn.pdf.txt",
        pu.get_path_to_default_graphs() + "timescales_v_Teff_He_bn_kn.pdf.txt",
    ],
}


def multipanelise_from_dumps(list_of_dump_configs):
    graph_fac = gf.GraphFactory()
    for dump_config in list_of_dump_configs:
        graph_fac.multipanelise_from_dumps(
            dump_config["list_of_dumps"],
            dump_config["y_dimension"],
            dump_config["x_dimension"],
            dump_config["filenames"],
            dump_config["fig_height"],
            dump_config["fig_width"],
            dump_config["gridspec_wspace"],
            dump_config["gridspec_hspace"],
            dump_config["sharey_axes"],
            dump_config["sharex_axes"],
        )


def main():
    list_of_dump_configs = [dump_config11]
    multipanelise_from_dumps(list_of_dump_configs)


if __name__ == "__main__":
    main()
