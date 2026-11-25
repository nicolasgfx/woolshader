# CPU cases carry one of two labels. retained_cpu: the CPU code of a feature that stays CPU-only
# (the along-yarn fibre generator behind the Fibres-step reference preview, the budget scheduling
# with its support release, and the shared yarn sampling).
# cpu_reference: a CPU implementation with a GPU counterpart on the workflow path (CPU settle
# and gate (d), the CPU tiled knit, legacy fibre growth). Both are always registered; a test is
# removed only together with the behaviour it guards. A case that SKIPs exits with 77
# (SKIP_RETURN_CODE below), so CTest never counts a skip as a pass.
# Each test program belongs to one component and sees that component's private headers and dependencies
# (knit_component_private_access, cmake/components.cmake); knit_tests_workflow, the integration program,
# names the components it inspects. A case keeps its CTest name when it moves to another program:
# knit_register_cases(<program> <prefix> <label or ""> cases...) registers <prefix><case>.
function(knit_register_cases program prefix label)
    foreach(case ${ARGN})
        add_test(NAME ${prefix}${case} COMMAND ${program} ${case})
        if(label)
            set_tests_properties(${prefix}${case} PROPERTIES LABELS ${label})
        endif()
    endforeach()
endfunction()

# knitting: curves, checks, layouts, tubes, plies, the CPU fibre generator, the explicit problem.
add_executable(knit_tests tests.cpp cpu_numerics.cpp ncp_reference.cpp)
knit_component_private_access(knit_tests knitting)
target_link_libraries(knit_tests PRIVATE knit_fibres knit_io)
knit_warnings(knit_tests)
target_compile_definitions(knit_tests PRIVATE KNIT_TEST_DATA="${CMAKE_SOURCE_DIR}/tests/data")
# The renderer, for the strand load test (ppt_load_tiny is skipped without it).
# The separate renderer for the cases that load a package (photon_path_tracer ppt.exe): KNIT_PPT_EXE from the
# environment at configure time; without it those cases skip.
set(KNIT_PPT_EXE "$ENV{KNIT_PPT_EXE}" CACHE FILEPATH "photon_path_tracer ppt.exe")
# One CTest entry per case; the list is printed by `knit_tests --list`.
set(KNIT_TEST_CASES
    chart_skew_measures_sheared_columns spiral_placement_pinned
    export_tube_accepts_matches_check
    bspline_reference bspline_derivatives bspline_continuity bezier_split
    arclength_line arclength_helix arclength_roundtrip arclength_bounds
    sampling_bound
    extrema_speed_interior extrema_curvature_interior
    cross_two_close cross_tangent cross_at_boundary cross_collinear_overlap
    cross_single_span_loop cross_adjacent_span cross_smooth_join
    cross_open_end_extension cross_root_box_spans_knot cross_depth_contraction
    clear_parallel_states clear_intersecting clear_open_end_self clear_hairpin_ok
    clear_hairpin_tight clear_distant_self clear_window_straddle
    stockinette_signature stockinette_periodicity stockinette_scaling_invariance
    presets_expected_verdicts
    rmf_frames tube_uv_seamless compression_two_sided compression_only_where_overlapping
    surface_fold_detected surface_new_collision
    compression_fold_limit compression_parallel_066 compression_third_yarn
    presets_compressed_surface_valid
    motion_certificate
    chart_seam_continuity chart_orientation topology_flat_matches_spans
    hash_matches_bruteforce schur_local_tight_bend length_bounds_prefix
    stitch_queue_balance count_schedule_feasible dec_turn_boundary planar_mapping_equals_flat
    local_signature_matches_flat spiral_cylinder_valid spiral_turn_transition same_yarn_interlock_signature
    flat_inc_signature flat_dec_signature spiral_cone_shaping_signature
    capped_sphere_poles spiral_end_gathered scene_piece_gathers_only_a_visible_end layout_distortion_bounds settle_keeps_topology_on_sphere settle_is_one_relax_round
    stitch_spread_measures_one_stray_control hand_tensions_are_pure hand_columns_slide_apart hand_rate_per_cm2 hand_raised_loop_is_loose hand_thickness_follows_tension_and_spots hand_rest_has_no_jump hand_rest_counts_the_tension_once settle_hand_frees_the_whole_fabric
    invalid_init_rejected repair_with_failing_signature_validates_nothing diameter_monotone_fixed_centreline repair_one_certified_move repair_start_diameter_from_the_report final_surface_gate
    pptstrands_roundtrip
    inside_mesh_matches_winding
    strand_budget_thinning
    ply_interpolation curvature_degenerate_hodograph
    ply_profile_circles ply_profile_bridge tube_uv_wrap tube_uv_jitter
    hash_canonical_cell_order hash_concurrent_queries topology_stitch_near
    mesh_closest_skip_revisits executor_basics executor_budget_and_nesting executor_exceptions
    checks_parallel_identical
    validate_parallel_identical
    ncp_jacobian_fd_edge ncp_jacobian_fd_curv ncp_jacobian_fd_contact ncp_jacobian_fd_band ncp_jacobian_fd_cap
    ncp_contact_sign_support ncp_fixed_eliminated ncp_curv_degenerate_flagged ncp_band_owned_only
    ncp_refresh_cumulative_motion ncp_model_vs_score_divergence

    ncp_same_yarn_arc_gap ncp_hard_target_split ncp_pass_through_detected
    attribution_crossing_rules attribution_close_rules local_signature_attribution_moved_boundary
    # Not registered on purpose (a tool, not a test): fibers_write_snapshots.
    # Audited 2026-10-25 against the runtime (not CPU-reference: shared or without a GPU counterpart):
    fibre_islands_cover_what_they_ask fibre_islands_zero_weight_unchanged
    fibers_deterministic fibers_front_only fibers_occlusion fibres_cut_before_deep_contact
    # - the per-candidate draws: the specification both generators implement (gpu_fibres_* compares against it)
    fibers_target_length_distribution fiber_draw_orientation fiber_draw_thickness
    fibers_length_truncated_lognormal fibers_mix_invariance fibers_curled_strays fiber_random_draws fiber_stray_spread
    # - CPU plies and tube mesh: used with the GPU too (tiled ply tubes, ply maps, the export's and the display's tubes)
    ply_containment ply_single_equals_tube ply_curvature_fold plies_parallel_identical
    tube_mesh_parallel_identical
    # - plies: the parent yarn and the shared candidate allocation (FiberFrame::candidates)
    fibers_parent_yarn_on_plies fibers_density_per_parent
)
foreach(case ${KNIT_TEST_CASES})
    add_test(NAME ${case} COMMAND knit_tests ${case})
endforeach()
set(KNIT_TEST_RETAINED_CPU_CASES
    fibers_root_bridge_serialized
    fibers_surface_all_sides fibers_segment_interior_surface
    fibers_root_attached_exported fibers_seed_per_piece
    fiber_curve_derivative fiber_curve_yarn_coords fiber_polyline_fit fibers_bent_all_accepted
    fiber_emergence_actual fibers_yarn_end_rejected
    fibers_on_plies
    fibers_mix_endpoints fibers_mix_on_plies fibers_random_shape fibers_random_truncation
    fibers_snapshots fibers_parallel_identical fibers_serial_policy_single_thread
    sampling_parallel_identical finalize_parallel_identical
)
foreach(case ${KNIT_TEST_RETAINED_CPU_CASES})
    add_test(NAME ${case} COMMAND knit_tests ${case})
    set_tests_properties(${case} PROPERTIES LABELS retained_cpu)
endforeach()
set(KNIT_TEST_CPU_REFERENCE_CASES
    settle_repair_parallel_identical gate_d_parallel_identical
)
foreach(case ${KNIT_TEST_CPU_REFERENCE_CASES})
    add_test(NAME ${case} COMMAND knit_tests ${case})
    set_tests_properties(${case} PROPERTIES LABELS cpu_reference)
endforeach()
# S6 cases (tests_s6.cpp): a separate executable, listed by `knit_tests_s6 --list`.
add_executable(knit_tests_s6 tests_s6.cpp cpu_numerics.cpp)
knit_component_private_access(knit_tests_s6 scene)
target_link_libraries(knit_tests_s6 PRIVATE knit_jobs knit_export)
knit_warnings(knit_tests_s6)
target_compile_definitions(knit_tests_s6 PRIVATE KNIT_TEST_DATA="${CMAKE_SOURCE_DIR}/tests/data"
                           KNIT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
set(KNIT_S6_TEST_CASES
    fur_switch_refuses_the_fur_workflow scene_settings_roundtrip scene_settings_guided_rows scene_settings_retired_keys scene_settings_retired_knit_keys scene_settings_retired_hair_keys legacy_knit_settings_keep_their_keys scene_settings_refusals scene_settings_stage_keys scene_settings_stage_inputs fnv1a_reference
    knit_params_text_pinned knit_stage_key_pinned ply_map_key_pinned hand_off_leaves_the_keys_alone hand_presets_load_with_an_off stray_presets_load_and_link fibre_stage_text_pinned fibre_hug_keys_only_while_set fibre_drape_preset_transitions tiled_classify_failure_ends_the_part
    texture_profile_parse fibre_presets_assign_all_values yarn_fibres_presets_assign_all_values fibres_stay_in_front object_budgets_are_independent
    pattern_stockinette_bit_identical pattern_reverse_is_mirror purl_signature_expectation pattern_flat_references pattern_rest_scale pattern_rib_layers pattern_garter_wave signature_rib_layer_overlap
    glb_load_transforms glb_required_extension_is_named glb_weld_keeps_corner_uv glb_texture_transform_reaches_the_colour_uv triangle_intersection_cases solid_eligibility_self_intersecting
    scene_cache_roundtrip scene_cache_key viewer_defaults_compute builder_cache_hit_skips_job builder_progressive_publish
    parallel_for_propagates_cancel_token
    filling_open_sheet_side filling_offset_inversion_reported cross_component_clearance_reported
    final_surface_e0_bit_identical inflate_flattening_known_gap reference_check_unresolved_when_compressed tiled_rendered_part_plies
    plies_on_compressed_yarn overlap_coverage_measures_plied_yarn gap_coverage_ignores_tube_settings scene_solids_made_at_load swatch_preset_loads_on_its_uv
    builder_rapid_edits_publish_latest_only builder_scene_switch_cancels builder_stale_completion_ignored builder_shutdown_while_running_joins builder_stage_exception_reported builder_export_not_cancelled_by_edit
    fibre_random_overrides_inherit_identical fibre_random_overrides_applied fibre_thickness_minmax_roundtrip fibre_angle_nominal_equals_physical_straight fibre_angle_measured_on_curved_bowed
    fibre_hug_levels_off fibre_follows_ply_helix
    export_cancel_keeps_previous_output export_failure_keeps_previous_output plymap_cancel_never_cache_hit plymap_corrupt_marker_not_hit publish_recovery_from_old_dir staging_cleanup_on_startup
    two_sheets_fold settle_keeps_the_relaxed_state_whole repair_rejected_attempt_keeps_projection_state
    uv_tiling_cylinder_valid uv_tiling_seam_clearance uv_tiling_mirrored_chart uv_chart_degenerate_folded_overlap lattice_runs_split_merge lattice_program_check mesh_surface_normal_continuous_across_seam boundary_check_certified_vs_sampled
    uv_tiling_plane_equals_flat uv_metric_affine_exact uv_scale_whole_chart_invariant uv_local_stretch_flagged
    texel_bilinear_linear_light palette_reduction
    knit_cache_pipeline tiled_validation_cancel export_closeup_cameras teddy_options_from_settings
    scene_settings_ncp_key_names_the_solver teddy_cache_identity teddy_twin_pieces_chart_alike teddy_defaults_one_path
    coverage_flat_matches_mesh_open_area coverage_domain_fixed coverage_regions_pinned coverage_missing_course coverage_thin_double_sheet coverage_chart_seam
    deadline_marks_part_unresolved_and_rest_not_started knit_parts_in_scene_order
    compress_two_yarns compress_gate_d_clamped_reported max_flattening_setting fold_radius_cpu_gpu_agree
    relax_thresholds_follow_validation capside_fixed_bit_identical
    # Not registered: overlap_seam_covered (over 900 s, stopped unfinished after 14 min on 2026-10-28).
)
foreach(case ${KNIT_S6_TEST_CASES})
    add_test(NAME s6_${case} COMMAND knit_tests_s6 ${case})
endforeach()
if(WIN32)
    # The external hard time cap (tests/capped_run.ps1).
    add_test(NAME capped_run_kills_and_reports
             COMMAND powershell -NoProfile -ExecutionPolicy Bypass -File ${CMAKE_SOURCE_DIR}/tests/capped_run_test.ps1
                     -Script ${CMAKE_SOURCE_DIR}/tests/capped_run.ps1)
endif()
set(KNIT_S6_TEST_CPU_REFERENCE_CASES
    tiled_parallel_identical
)
foreach(case ${KNIT_S6_TEST_CPU_REFERENCE_CASES})
    add_test(NAME s6_${case} COMMAND knit_tests_s6 ${case})
    set_tests_properties(s6_${case} PROPERTIES LABELS cpu_reference)
endforeach()
knit_register_cases(knit_tests_s6 "" "" knit_cache_round_trip knit_cache_attempt_selection ply_maps_periodic)  # from knit_tests

# export: the neutral scene, the package and its materials, the export's tubes (tests_export.cpp), and the
# frozen bytes of every stored format (tests_io_pins.cpp). `knit_tests_export --list`.
add_executable(knit_tests_export tests_export.cpp tests_io_pins.cpp harness_main.cpp cpu_numerics.cpp)
knit_component_private_access(knit_tests_export export)
target_link_libraries(knit_tests_export PRIVATE knit_import)
knit_warnings(knit_tests_export)
target_compile_definitions(knit_tests_export PRIVATE KNIT_TEST_DATA="${CMAKE_SOURCE_DIR}/tests/data" KNIT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
if(EXISTS "${KNIT_PPT_EXE}")
    target_compile_definitions(knit_tests_export PRIVATE KNIT_PPT_EXE="${KNIT_PPT_EXE}")  # ppt_load_tiny skips without it
endif()
knit_register_cases(knit_tests_export "" "" pptstrands_refusals late_taper_keeps_the_body_and_the_path late_taper_segment_count_exact late_taper_modes_keep_the_path late_taper_existing_is_b0 late_taper_onset_keeps_the_body late_taper_onset_relative_to_emergence late_taper_refined_error_bound late_taper_float_ratio_margin late_taper_sparse_rules_scoped late_taper_fallback_profile renderer_memory_matches_the_teddy closeup_part_default_is_the_largest_knitted ppt_export_wool_look_keys ppt_tube_winding_outward ppt_export_bounds_ground ppt_load_tiny ppt_load_grouped_strict export_mix_partition ppt_fibre_material_is_yarn ppt_fibre_material_resolution ppt_export_fibre_groups ppt_unplaced_bytes_unchanged ppt_material_schema12 ppt_material_table_lists_only_referenced ppt_fibre_parts_untextured)  # from knit_tests
knit_register_cases(knit_tests_export s6_ "" tube_split_partition glb_units_roundtrip glb_object_identity io_byte_reader_refuses io_pin_settings_text io_pin_stage_keys io_pin_import_identity io_pin_mesh_read io_pin_prepared_import io_pin_knit_formats io_old_fibre_files_refused io_pin_validation_snapshot io_pin_scene_json io_pin_obj_files io_pin_report_json export_report_escapes_the_settings export_every_format_writes ply_map_set_cache baked_map_set_cache)  # from knit_tests_s6

# import: the wrap through the import session (tests_wrap.cpp) and the properties of each built-in stage
# (tests_wrap_properties.cpp, the flatten in tests_wrap_flatten.cpp). `knit_tests_import --list`.
add_executable(knit_tests_import tests_wrap.cpp tests_wrap_properties.cpp tests_wrap_flatten.cpp harness_main.cpp)
knit_component_private_access(knit_tests_import import)
target_link_libraries(knit_tests_import PRIVATE knit_jobs knit_scene)
knit_warnings(knit_tests_import)
target_compile_definitions(knit_tests_import PRIVATE
    KNIT_TEST_DATA="${CMAKE_SOURCE_DIR}/tests/data" KNIT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
knit_register_cases(knit_tests_import workflow_ "" import_split_by_nodes import_selection_keep_file import_selection_built_in import_selection_of_nothing_is_refused import_selection_identity_and_keys import_identity_unchanged import_wrap_identity_differs import_settings_old_file import_stage_keys_unchanged import_wrap_roundtrip import_wrap_one_primitive import_wrap_session_cancel import_wrap_session_reuses_a_stored_wrap import_retired_engine_keys_load import_keep_file_fur_key_unchanged colourwork_past_the_chart_border wrap_quad_split_planar wrap_quad_split_folded wrap_quad_split_uv_follow_corners wrap_source_soup_y_up_turn wrap_source_soup_duplicate_dropped wrap_source_soup_quads_resplit wrap_remesh_cube wrap_remesh_inside_out wrap_chart_rules_split_a_mostly_stretched_chart wrap_sorted_edges wrap_pack_rectangles_stand_up wrap_pack_islands_do_not_overlap wrap_pack_uv_in_unit_square wrap_pack_gutter_holds wrap_pack_unflagged_untouched wrap_closest_point_regions wrap_closest_point_tie_lower_index wrap_closest_point_brute_force wrap_closest_point_empty_and_not_finite wrap_vertex_normals_cube wrap_surface_file_roundtrip wrap_surface_file_keeps_the_flattened_triangles wrap_whole_wrap_box wrap_whole_wrap_boxes_apart wrap_whole_wrap_touching_boxes wrap_whole_wrap_torus flatten_planar_disc_is_isometric flatten_hemisphere_has_no_flip flatten_boundary_out_of_order flatten_zero_length_edge_mollified flatten_annulus_fails flatten_forced_solver_failure_keeps_start_map flatten_same_bytes_at_any_thread_count flatten_start_map_scaled_to_the_chart flatten_sphere_with_hole_within_bounds flatten_capped_tube_converges flatten_iteration_cap_reported)  # from knit_tests_workflow

# stitchmesh: the stockinette stitch-mesh experiment (tests_stitchmesh.cpp): the reference surface, Method A's time
# field and rings, Method B's orientation field; the construction fixtures, negative controls, virtual shift, guided
# rows and coverage. `knit_tests_stitchmesh --list`.
add_executable(knit_tests_stitchmesh tests_stitchmesh.cpp tests_stitchmesh_markings.cpp tests_stitchmesh_yarn.cpp tests_assisted.cpp harness_main.cpp)
target_link_libraries(knit_tests_stitchmesh PRIVATE knit_stitchmesh knit_model)
knit_warnings(knit_tests_stitchmesh)
target_compile_definitions(knit_tests_stitchmesh PRIVATE KNIT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")  # the frozen settle settings file
knit_register_cases(knit_tests_stitchmesh "" "" ssm_revolution_mesh_valid ssm_revolution_variant_j ssm_time_field_planar_front_exact
    ssm_time_field_tube_matches_arc_length ssm_rings_one_closed_contour_per_level ssm_rings_report_open_contours_at_a_hole
    ssm_orientation_guided_equals_guide ssm_fixtures_as_expected ssm_fixture_open_tube_virtual_helix ssm_fixture_kp_and_d12k
    ssm_fixture_special_quad_strands_stay_closed ssm_negative_controls_name_their_defect ssm_guided_rows_tube_case_d ssm_guided_rows_knit_over_a_bump_too_small_for_a_round
    ssm_relax_forces_balance ssm_relax_rest_and_perturbation ssm_relax_virtual_helix ssm_relax_stall_rule_stops ssm_helix_boundary_contract
    ssm_trim_tube_plane_curve ssm_path_guide_equals_half_plane ssm_path_guide_double_crossing ssm_marked_capsule_build ssm_knit_components_overlap_on_one_piece
    ssm_marked_pole_route_finish_level ssm_marked_pole_route_yarn_is_clean ssm_marked_offered_preview_rejected ssm_prepared_reference_maps_exactly ssm_walk_residual_on_a_crease
    ssm_markings_round_trip_and_refusals ssm_markings_rebind_onto_a_refinement ssm_arm_fixtures_as_expected ssm_marked_rows_follow_the_path ssm_relax_guides_hold ssm_border_ring_is_placed_on_the_bends ssm_relax_border_ring_keeps_its_order ssm_relax_takes_back_steps_that_raise_the_forces
    ssm_border_ring_first_corner_may_move ssm_marked_rows_turned_for_upright_wales ssm_marked_bind_off_strip_is_plain
    ssm_loop_reference_limit_points ssm_marked_loop_preparation ssm_shaping_assessment_model_and_defect
    ssm_settle_failed_raw_construction_not_solved ssm_settle_clean_construction_enters_the_solver ssm_settle_broken_strand_refused_always ssm_settle_failed_converted_construction_not_solved ssm_settle_stop_rule_after_a_target_call ssm_settle_without_progress_ends_stalled ssm_yarn_regular_patch ssm_yarn_adjacent_increases_need_the_face_placement ssm_yarn_shaping_boundary_pole ssm_yarn_checks_fail_with_the_defect ssm_yarn_curves_conversion ssm_yarn_own_plane_frame ssm_yarn_calibrated_templates
    guided_records_from_stitch_mesh guided_settle_defaults_frozen
    ssm_far_finish_named_and_not_offered ssm_proposals_open_tube ssm_proposals_bent_inner_outer ssm_proposals_transform_invariant ssm_proposals_capsule_partial ssm_proposals_branch_refused
    ssm_cup_pole_start_opening_finish ssm_cup_opening_is_a_level_of_the_rows ssm_marked_open_tube_border_loops ssm_row_level_through_vertices ssm_aligned_regular_tube ssm_quick_review_row_crossings ssm_quick_review_matches_layout ssm_quick_review_incomplete_and_self_crossing
    ssm_quick_review_stale_field ssm_quick_review_measures_known)

# markings: surface markings below the mesher (tests_markings.cpp): paths on a mesh, the pieces of an imported file,
# the markings file. `knit_tests_markings --list`.
add_executable(knit_tests_markings tests_markings.cpp harness_main.cpp)
target_link_libraries(knit_tests_markings PRIVATE knit_io knit_model)
knit_warnings(knit_tests_markings)
knit_register_cases(knit_tests_markings "" "" markings_trace_flat_grid markings_trace_cylinder_wrap markings_trace_vertex_points_take_the_short_way markings_pieces_weld_uv_seams
    io_markings_roundtrip io_markings_refuses_unknown_key io_pin_markings_text io_markings_reads_v1_writes_v2 io_markings_reads_v2_writes_v3
    io_markings_reads_v3_writes_v4 io_markings_reads_v4_writes_v5 markings_knit_components_apart markings_setup_accept_reject markings_setup_role_conflicts markings_preview_shadows_roles markings_merge_keeps_user
    markings_border_loop_detection)

# The architecture: what the compiler cannot hold, checked on the sources and on build/components.json.
add_executable(knit_tests_architecture tests_architecture.cpp harness_main.cpp)
target_link_libraries(knit_tests_architecture PRIVATE knit_model)  # harness.h reads the environment (model/env.h)
knit_warnings(knit_tests_architecture)
target_compile_definitions(knit_tests_architecture PRIVATE KNIT_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
    KNIT_COMPONENTS_JSON="${CMAKE_BINARY_DIR}/components.json")
knit_register_cases(knit_tests_architecture s6_ "" architecture_import_doc_complete architecture_import_layers architecture_jobs_doc_complete architecture_jobs_layers architecture_every_test_registered architecture_export_door architecture_fibres_doc_complete architecture_fibres_layers architecture_gpu_cuda_calls_checked architecture_gpu_doc_complete architecture_gpu_layers architecture_knitting_doc_complete architecture_gpu_no_templates architecture_knitting_layers architecture_knitting_no_file_io architecture_model_doc_complete architecture_io_doc_complete architecture_settings_doc_complete architecture_stitchmesh_doc_complete architecture_export_doc_complete architecture_graph_table_current architecture_model_is_plain architecture_model_layers architecture_module_dependencies architecture_no_new_generic_lambdas architecture_no_new_templates architecture_no_stale_paths architecture_scene_doc_complete architecture_scene_layers architecture_viewer_doc_complete architecture_viewer_layers)  # from knit_tests_s6

# Publication safety (tests_publish.cpp): `knit_tests_publish --list`.
add_executable(knit_tests_publish tests_publish.cpp)
target_link_libraries(knit_tests_publish PRIVATE knit_io knit_export)
knit_warnings(knit_tests_publish)
set(KNIT_PUBLISH_TEST_CASES
    publish_recovery_keeps_foreign_siblings publish_recovery_owned_leftovers publish_recovery_token_of_other_target
    publish_refuses_untracked_files publish_refuses_foreign_target publish_keeps_parent_files publish_export_name_rules
    publish_refuses_junction_target publish_into_folder_keeps_other_files publish_into_folder_creates_and_aborts
    publish_into_folder_failure_restores publish_into_folder_open_file publish_into_folder_recovery publish_into_folder_refuses_a_file
)
foreach(case ${KNIT_PUBLISH_TEST_CASES})
    add_test(NAME publish_${case} COMMAND knit_tests_publish ${case})
endforeach()
# The viewer workflow's backend (tests_workflow.cpp): import, results, eligibility,
# export formats, cancellation. `knit_tests_workflow --list`.
add_executable(knit_tests_workflow tests_workflow.cpp tests_guided.cpp tests_marking_session.cpp cpu_numerics.cpp)
knit_integration_test_access(knit_tests_workflow jobs scene export)
knit_warnings(knit_tests_workflow)
# cli_hair_export_has_no_knit runs the command-line tool itself.
add_dependencies(knit_tests_workflow knit_check)
target_compile_definitions(knit_tests_workflow PRIVATE KNIT_CHECK_EXE="$<TARGET_FILE:knit_check>" KNIT_TEST_DATA="${CMAKE_SOURCE_DIR}/tests/data")
if(EXISTS "${KNIT_PPT_EXE}")
    target_compile_definitions(knit_tests_workflow PRIVATE KNIT_PPT_EXE="${KNIT_PPT_EXE}")  # cli_yarn_body_export renders with it
endif()
set(KNIT_WORKFLOW_TEST_CASES
    import_nodes_from_names import_nodes_from_loose_parts import_session_refuses_too_many_nodes import_keep_file_knits_on_charts export_selection_keeps_the_plain_node
    import_obj_material_texture_colours import_ply_three_encodings import_units_and_up_axis import_oversize_rejected_before_parse
    import_header_counts_checked_before_allocation import_malformed_fatal_messages import_repairs_exact_counts
    import_missing_texture_falls_back import_identity_dependencies import_store_roundtrip
    workflow_hand_solver_adds_height workflow_hand_thickness_survives_a_reload workflow_plane_knit_fibres_export export_obj_contents_and_coordinates export_eligibility_rows export_job_enforces_eligibility
    tiled_knit_cache_roundtrip tiled_knit_cache_previous_version_and_knit_again viewer_load_knit_is_used_or_fails tiled_knit_cache_gap_patches fibre_progress_aggregator fresh_mesh_knit_then_fibres_same_session strand_mesh_rebuild_contract builder_cancel_any_job_drops_pending thickness_presets_flat_matrix custom_profile_is_not_final thick_yarn_warning stitch_size_from_yarn gauge_choices_round_trip
    obj_export_mixed_attributes import_session_reset_race import_session_cancel_queued_prepare import_session_stored_match
    yarn_look_roundtrip yarn_look_fibres_construction fibre_progress_and_profile_lines fibre_options_do_not_change_strands
    fibre_budget_matches_full_per_backend pile_keys_and_fur_key fur_field_properties fur_irregularities fur_loops_properties
    workflow_plane_fur_export fibre_preview_reference fur_coverage_rules
    workflow_key_and_migration cli_fur_export_has_no_knit cli_tiled_knit_save_load cli_yarn_body_export cli_yarn_body_maps_export strand_mesh_yarn_maps blender_export_baked_yarn saved_fibre_key_check cli_pieces_manifest_and_saved_fibres cli_retired_research_flags preset_fur_export_like_an_import every_preset_names_its_model blue_swatch_preset yarn_body_refusal_blocks_export pieces_knit_fibres_export fibre_budget_per_object
    import_default_yarn_colour_source presets_open_with_their_texture_colours fur_takes_the_texture_colour_at_its_root stitch_colours_average_each_cell stitch_colours_across_a_chart_seam stitch_colours_do_not_depend_on_the_source
    stitch_colours_along_tubes_and_plies export_texture_colourwork stitch_colours_without_a_texture
    export_unknitted_colour_sources export_tube_triangles_follow_the_tubed_length preset_export_writes_its_details everything_unknitted_exports_plain_meshes export_default_format blender_export_writes_a_blend blender_export_keeps_the_credit blender_export_translucent_yarn_and_cameras blender_export_job_publishes_the_blend blender_export_cancel_stops_blender blender_export_reports_the_script_message partial_base_inside_fabric knit_stand_ins_survive_the_second_part_list teddy_spiral_pieces_result_path_synthetic_final teddy_spiral_pieces_time_cap teddy_pieces_ply_tubes
    # Not registered, run by hand: the evidence tools viewer_chain_evidence and diag_preset_sweep.
    export_writes_into_the_chosen_folder teddy_pieces_drawn_diameter
    guided_source_refuses_other_file guided_ignores_offered_and_rejected guided_setup_fields_not_input guided_stage_boundaries guided_refused_construction guided_cancel_publishes_nothing guided_knit_key_follows_the_guided_version guided_draft_without_a_solve
    guided_settled_yarn_with_failing_mesh_checks_is_not_called_unsettled guided_inputs_keep_the_whole_solver_identity guided_knit_component_enters_the_knit_key
    marking_quick_not_behind_build marking_cancel_superseded marking_quick_latest_only
    stitch_mesh_use_gate stitch_mesh_failed_rebuild_replaces stitch_mesh_obj_matches_build
    # Not registered on purpose (tools): diag_preset_sweep, diag_fur_baseline_measures, diag_fur_preview_latency.
)
foreach(case ${KNIT_WORKFLOW_TEST_CASES})
    add_test(NAME workflow_${case} COMMAND knit_tests_workflow ${case})
endforeach()
# The default with no Blender on the machine: the package, and the Blender entry says why it is off.
add_test(NAME workflow_export_default_format_without_blender COMMAND knit_tests_workflow export_default_format)
set_tests_properties(workflow_export_default_format_without_blender PROPERTIES ENVIRONMENT "KNIT_BLENDER=C:/no/such/blender.exe")
set_tests_properties(workflow_teddy_spiral_pieces_result_path_synthetic_final PROPERTIES TIMEOUT 2400)
set_tests_properties(workflow_pieces_knit_fibres_export PROPERTIES TIMEOUT 1500)
# knit_check --time-cap: the hard cap stops a long run (exit code 4) within seconds of the cap.
add_test(NAME cli_time_cap COMMAND knit_check --scene wool-swatch --time-cap 1)
set_tests_properties(cli_time_cap PROPERTIES PASS_REGULAR_EXPRESSION "time cap of 1 s reached: run stopped"
                     SKIP_REGULAR_EXPRESSION "\\.glb not found" TIMEOUT 30)
# The yarn body of structural fibres (tests_yarn_body.cpp, src/fibres/yarn_body.h): `knit_tests_yarn_body --list`.
add_executable(knit_tests_yarn_body tests_yarn_body.cpp harness_main.cpp)
knit_component_private_access(knit_tests_yarn_body fibres)
knit_warnings(knit_tests_yarn_body)
target_compile_definitions(knit_tests_yarn_body PRIVATE KNIT_TEST_DATA="${CMAKE_SOURCE_DIR}/tests/data")
knit_register_cases(knit_tests_yarn_body "" ""
    yarn_body_radii_and_rules yarn_body_fibre_angle_straight yarn_body_continuous_through_bends yarn_body_deterministic_unique
    yarn_body_follows_envelope yarn_body_regular_definition yarn_body_interior_keeps_shell yarn_body_validation
    yarn_body_overlap_measure yarn_body_sampling_independent yarn_body_offsets_64bit yarn_body_cancel
    yarn_bake_straight_regular yarn_bake_gaps_shadowed yarn_bake_atlas_layout yarn_bake_uv_matches_tube yarn_bake_material_maps)
# Fibre CPU plan cases (tests_fibre.cpp): a separate executable, `knit_tests_fp --list`.
add_executable(knit_tests_fp tests_fibre.cpp cpu_numerics.cpp)
knit_component_private_access(knit_tests_fp fibres)
knit_warnings(knit_tests_fp)
target_compile_definitions(knit_tests_fp PRIVATE KNIT_TEST_DATA="${CMAKE_SOURCE_DIR}/tests/data")
set(KNIT_FP_TEST_CASES
    yarn_position_equals_coordinates
    fibre_params_text_pinned fur_strands_pinned
)
foreach(case ${KNIT_FP_TEST_CASES})
    add_test(NAME fp_${case} COMMAND knit_tests_fp ${case})
endforeach()
set(KNIT_FP_TEST_RETAINED_CPU_CASES
    budget_sched_matches_full_synthetic budget_sched_exact_exhaustion
    budget_sched_ties_across_pieces budget_sched_key_rounds_to_one
    budget_sched_empty_and_rejection_heavy budget_sched_protocol_enforced budget_sched_windows_pinned
    support_inputs_equal_reuse
)
foreach(case ${KNIT_FP_TEST_RETAINED_CPU_CASES})
    add_test(NAME fp_${case} COMMAND knit_tests_fp ${case})
    set_tests_properties(fp_${case} PROPERTIES LABELS retained_cpu)
endforeach()
# GPU tests (tests/gpu/*.cpp): cases needing a device SKIP without one.
file(GLOB KNIT_GPU_TEST_SOURCES CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/tests/gpu/*.cpp)
add_executable(knit_tests_gpu ${KNIT_GPU_TEST_SOURCES} harness_main.cpp cpu_numerics.cpp ncp_reference.cpp)
knit_component_private_access(knit_tests_gpu gpu)
target_link_libraries(knit_tests_gpu PRIVATE knit_scene knit_jobs knit_settings CUDA::cudart_static)  # test_resident.cpp runs the job queue
knit_warnings(knit_tests_gpu)
target_compile_definitions(knit_tests_gpu PRIVATE KNIT_TEST_DATA="${CMAKE_SOURCE_DIR}/tests/data"
                           KNIT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
# One CTest entry per case (`knit_tests_gpu --list`), named gpu_<case> (gpu_ is not doubled),
# read from the TEST(...) lines of the compiled sources, so a known failing case cannot hide a
# new one.
foreach(src ${KNIT_GPU_TEST_SOURCES})
    file(STRINGS ${src} gpu_test_lines REGEX "^TEST\\([a-z0-9_]+\\)")
    foreach(line ${gpu_test_lines})
        string(REGEX REPLACE "^TEST\\(([a-z0-9_]+)\\).*" "\\1" case "${line}")
        if(case MATCHES "^gpu_")
            set(name ${case})
        else()
            set(name gpu_${case})
        endif()
        add_test(NAME ${name} COMMAND knit_tests_gpu ${case})
    endforeach()
endforeach()

if(KNIT_BUILD_VIEWER)
    # UI tests (tests_ui.cpp): the scene panel without a window or GL, on an ImGui built
    # with the test-engine hooks and a throwing IM_ASSERT (tests/imgui_test_config.h).
    add_library(imgui_test_lib STATIC
        ${imgui_SOURCE_DIR}/imgui.cpp
        ${imgui_SOURCE_DIR}/imgui_draw.cpp
        ${imgui_SOURCE_DIR}/imgui_tables.cpp
        ${imgui_SOURCE_DIR}/imgui_widgets.cpp
        ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
        ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl2.cpp)
    target_include_directories(imgui_test_lib PUBLIC ${imgui_SOURCE_DIR} ${imgui_SOURCE_DIR}/backends)
    target_compile_definitions(imgui_test_lib PUBLIC IMGUI_USER_CONFIG="${CMAKE_SOURCE_DIR}/tests/imgui_test_config.h")
    target_link_libraries(imgui_test_lib PUBLIC glfw)
    if(MSVC)
        target_compile_options(imgui_test_lib PRIVATE /W0)
    endif()
    add_executable(knit_tests_ui
        tests_ui.cpp
        ${CMAKE_SOURCE_DIR}/src/viewer/interface.cpp
        ${CMAKE_SOURCE_DIR}/src/viewer/about_window.cpp
        ${CMAKE_SOURCE_DIR}/src/viewer/widgets.cpp
        ${CMAKE_SOURCE_DIR}/src/viewer/panel_import.cpp
        ${CMAKE_SOURCE_DIR}/src/viewer/panel_mark.cpp
        ${CMAKE_SOURCE_DIR}/src/viewer/panel_knit.cpp
        ${CMAKE_SOURCE_DIR}/src/viewer/panel_fibres.cpp
        ${CMAKE_SOURCE_DIR}/src/viewer/panel_fur.cpp
        ${CMAKE_SOURCE_DIR}/src/viewer/panel_yarn_fibres.cpp
        ${CMAKE_SOURCE_DIR}/src/viewer/panel_export.cpp
        ${CMAKE_SOURCE_DIR}/src/viewer/fibre_lines.cpp
        ${CMAKE_SOURCE_DIR}/src/viewer/mesh_lighting.cpp)
    knit_component_private_access(knit_tests_ui viewer)
    target_link_libraries(knit_tests_ui PRIVATE imgui_test_lib glfw OpenGL::GL)
    target_compile_definitions(knit_tests_ui PRIVATE KNIT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
    knit_warnings(knit_tests_ui)
    foreach(case ui_panel_unique_ids_teddy ui_panel_unique_ids_wool_swatch ui_scene_settings_defaults_intended ui_unified_workflow_steps ui_fur_switched_off ui_export_blocker_visible ui_overlay_counts
            ui_progress_classification ui_knit_progress_single_home ui_knit_footer_reachable ui_knit_controls_inventory ui_knit_fibres_disables ui_fibre_yarn_feedback ui_knit_fibres_rendered_as ui_display_yarn_maps_toggle
            ui_construction_tiles ui_fibres_yarn_section ui_yarn_fibres_reveal ui_yarn_fibres_cost_counts_fuzz ui_yarn_map_bake_progress ui_yarn_map_material
            ui_knit_drag_dimension ui_knit_draft_final_actions ui_knit_mesh_detail_mesh_only ui_knit_gauge_from_yarn
            ui_knit_action_notes_per_profile ui_knit_custom_range_ends ui_knit_review_import_round_trip ui_knit_detail_summaries
            ui_export_controls_inventory ui_preset_with_knitting_skipped_exports_its_mesh ui_export_bindings ui_export_fibres_auto_include ui_export_target_matches ui_export_blender_format ui_export_footer_states
            ui_export_footer_reachable ui_export_progress_single_home ui_export_workflow_fixes
            ui_export_detail ui_same_workflow_all_scenes ui_teddy_pieces_scene_switch
            ui_fibre_controls_inventory ui_fibre_footer_reachable ui_fibre_footer_states ui_fibre_preset_custom ui_hand_preset ui_stray_preset ui_island_preset ui_fibre_density_preset
            ui_fibre_preview_target ui_fibre_preview_follows_edits ui_fibre_progress_single_home ui_fibre_share_percent ui_fibre_tabs_preserve
            ui_import_default_rebuilds ui_import_two_stages ui_import_pick_focuses ui_import_drag_is_no_click ui_view_ray_matches_the_camera ui_import_inventory ui_import_footer_reachable ui_import_revalidation ui_import_lacks_charts ui_import_workflow_switch_selects_source ui_fur_disc_preview ui_import_skip_forces_keep ui_import_failed_wrap_keeps_the_choice ui_import_failed_wrap_keep_file_on_request ui_import_opens_with_colourwork ui_start_without_object
            ui_import_long_list_footer_reachable ui_import_pick_reveals_row ui_import_size_readouts
            ui_import_focus_and_tick_are_separate ui_import_second_file_resets_draft ui_import_invalid_size ui_import_custom_survives_height
            ui_import_suggestion_names_unit ui_import_mixed_keeps_plain_component
            ui_export_writes_into_the_chosen_folder ui_text_snapshot ui_text_snapshot_local_models ui_flat_swatch_is_a_note ui_knit_patterns ui_imgui_ini_migration ui_ssao_kernel ui_camera_turntable ui_fibre_populations_checkbox ui_camera_drag_lock ui_shown_meshes_cover_view ui_viewer_options ui_script_steps ui_fibre_colour_sources ui_platform_requests ui_fibre_mask_image_changed ui_import_unit_choices ui_import_skip_rules ui_import_cancel_paths ui_import_rules_without_a_frame ui_import_scripted_skip_keeps_charts ui_knit_final_button ui_tubes_follow_level ui_stitch_size_on_change ui_fibre_chunks ui_layout_at_rest
            ui_camera_presets ui_view_toolbar_inventory ui_view_mode_availability ui_view_mode_open_file ui_view_mode_partial_rebuild ui_display_rows ui_compare_with_source ui_view_popovers ui_toolbar_narrow ui_render_preview_window
            ui_about_window ui_about_licence_texts ui_about_public_paths ui_about_credits
            ui_mark_stage_listed ui_mark_pick_adds_point ui_mark_undo_redo ui_mark_proposal_accept_reject ui_mark_save_load ui_mark_move_point_edits_proposal ui_mark_preparation_saved
            ui_mark_knit_this_piece ui_mark_build_outdated_after_edit ui_knit_route_section
            ui_mark_setup_preview_accept_undo ui_mark_preview_after_accept ui_mark_proposals_keep_knit_route ui_mark_setup_outdated_basis
            ui_mark_quick_feedback_states ui_mark_save_load_setups ui_mark_highlight_on_review_traces ui_mark_quick_display_matches_request
            ui_mark_stitch_view_in_toolbar ui_mark_stitch_view_after_build ui_mark_stitch_view_not_stolen ui_mark_stitch_view_switch_back
            ui_mark_stitch_view_empty_rebuild ui_mark_export_stitch_mesh ui_mark_stitch_view_click_selects_only)
        add_test(NAME ui_${case} COMMAND knit_tests_ui ${case})
    endforeach()
endif()

# The author's scene and prototype tools (tools/, not part of the published source): registered only where the
# folder and a Python interpreter are present. make_brown_manifest.py on a synthetic two-bear manifest (the brown
# teddy preset's pieces).
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND AND EXISTS ${CMAKE_SOURCE_DIR}/tools/scenes/make_brown_manifest.py)
    add_test(NAME tools_make_brown_manifest COMMAND ${Python3_EXECUTABLE} -I ${CMAKE_SOURCE_DIR}/tools/scenes/make_brown_manifest.py --self-test)
    set_tests_properties(tools_make_brown_manifest PROPERTIES PASS_REGULAR_EXPRESSION "self-test: PASS" TIMEOUT 60)
endif()
# tools/fur_noise: the spawn-point prototype's unittest suite (needs numpy, scipy, Pillow; the window test skips
# without a display). discover puts the absolute start directory on sys.path, which -I does not.
if(Python3_Interpreter_FOUND AND EXISTS ${CMAKE_SOURCE_DIR}/tools/fur_noise)
    add_test(NAME tools_fur_noise COMMAND ${Python3_EXECUTABLE} -I -m unittest discover -s ${CMAKE_SOURCE_DIR}/tools/fur_noise -p "test_*.py")
    set_tests_properties(tools_fur_noise PROPERTIES TIMEOUT 300)
endif()

# Every test runner exits with 77 when the case it ran was skipped (a missing asset or device):
# CTest reports it as skipped, never as passed.
get_property(knit_all_tests DIRECTORY PROPERTY TESTS)
set_tests_properties(${knit_all_tests} PROPERTIES SKIP_RETURN_CODE 77)
# The Fur workflow is off by default (kFurByDefault, not released yet); the tests cover both workflows, so it is
# switched on for every case. The cases about the switch set it themselves (setFurEnabled).
set_property(TEST ${knit_all_tests} APPEND PROPERTY ENVIRONMENT "KNIT_FUR=1")
